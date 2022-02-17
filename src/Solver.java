import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.ParseException;
import java.util.*;

/**
 * This class contains all methods for executing a multi-echelon FRBS with an automatic KB generation.
 * It provides methods for:
 * - Applying 3-step procedure (WM, SA, LA-Tuning) for training a global FRBS (exact model) and a local FRBS (heuristic).
 * - Testing models on entire SC
 */
public class Solver {
	private SC sc;
	private List<Inventory> inventories; // List of active inventories
	private Map<Inventory,List<Material>> materials;
	private KnowledgeBaseSet initKbSet;
	private KnowledgeBaseSet bestKbSet; // Currently best knowledge base set
	private KnowledgeBaseSet bestKbSetWM; // Best knowledge base after applying Wang-Mendel method
	private KnowledgeBaseSet bestKbSetSA; // Best knowledge base after applying Simulate Annealing
	private KnowledgeBaseSet bestKbSetLA; // Best knowledge base after applying LA tuning
	private Map<String,KnowledgeBase> bestKnowledgeBasesHeuristic = new HashMap<String,KnowledgeBase>();
	private int numLabels;
	// Parameters for Simulated Annealing
	private double alpha;
	private double temp;
	private double percentageSwitches;
	private int stoppingCount;
	private List<WMGeneration> wmObjects;
	private KnowledgeBaseSet kbSet;
	private Map<String,Double> totalCosts;
	private double totalScCosts;
	private double scFillRate;
	// Parameters for Genetic Tuning
	private int N;
	private int bitsGene;
	private double phi;
	private int numCycles;
	private int maxIterations;
	private double gammaFit;
	private double phiFit;
	
	private List<Inventory> plants;
	private List<Inventory> warehouses;
	private List<Inventory> distributionCenters;
	
	private Writer writer;
	private int scenario;
	// Lists for storing run times
	private List<Long> runtimesExact;
	private List<Long> runtimesHeuristic;
	
	/**
	 * Default constructor.
	 * @param sc SC object
	 * @param numLabels Number of fuzzy sets (labels) per fuzzy variable
	 * @param alpha Parameter for reducing temperature --> T_new = T_old * alpha
	 * @param temp Initial temperature
	 * @param percentageSwitches Percentage of switches that are applied in an iteration of the Simulated algorithm
	 * @param stoppingCount Number of iterations without improvement after which the Simulated algorithm stops
	 * @param N Number of individuals in population of LA-Tuning algorithm
	 * @param bitsGene Number of bits used for encoding gray codes in LA-Tuning algorithm
	 * @param phi Coefficient for reducing threshold of LA-Tuning algorithm
	 * @param numCycles Maximal number of cycles of LA-Tuning algorithm
	 * @param maxIterations Maximal number of iterations of LA-Tuning algorithm
	 * @param gammaFit Parameter for prioritizing SC costs
	 * @param phiFit Parameter for prioritizing fill rate 
	 * @param scenario Integer that identifies the current scenario of the sensitivity analysis
	 */
	public Solver(SC sc, int numLabels, double alpha, double temp, double percentageSwitches, int stoppingCount, int N, 
			int bitsGene, double phi, int numCycles, int maxIterations, double gammaFit, double phiFit, int scenario) {
		this.sc = sc;
		List<Inventory> allInventories = sc.getInventories();
		this.inventories = new ArrayList<Inventory>();
		// Only plants, warehouses and distribution centers are within the own SC --> are part of the FRBS
		for(Inventory i : allInventories) {
			if(i.getStageType().equals("PLANT") || i.getStageType().equals("WH") || i.getStageType().equals("DIST")) {
				this.inventories.add(i);
			}
		}
		this.materials = new HashMap<Inventory,List<Material>>();
		// Retrieve materials relevant for each inventory
		for(Inventory i : inventories) {
			List<Material> invMaterials = i.getMaterials();
			materials.put(i, invMaterials);
		}
		this.initKbSet = null;
		this.bestKbSet = null;
		this.bestKbSetWM = null;
		this.bestKbSetSA = null;
		this.bestKbSetLA = null;
		this.bestKnowledgeBasesHeuristic = new HashMap<>();
		this.numLabels = numLabels;
		this.alpha = alpha;
		this.temp = temp;
		this.percentageSwitches = percentageSwitches;
		this.stoppingCount = stoppingCount;
		this.N = N;
		this.bitsGene = bitsGene;
		this.phi = phi;
		this.numCycles = numCycles;
		this.maxIterations = maxIterations;
		this.gammaFit = gammaFit;
		this.phiFit = phiFit;
		this.wmObjects = new ArrayList<WMGeneration>();
		this.totalCosts = new HashMap<>();
		this.totalScCosts = 0.0;
		this.scFillRate = 0.0;
		
		this.plants = new ArrayList<Inventory>();
		this.warehouses = new ArrayList<Inventory>();
		this.distributionCenters = new ArrayList<Inventory>();
		for(Inventory i : inventories) {
			if(i.getStageType().equals("PLANT")) {
				plants.add(i);
			}
			if(i.getStageType().equals("WH")) {
				warehouses.add(i);
			}
			if(i.getStageType().equals("DIST")) {
				distributionCenters.add(i);
			}
		}
		this.writer = new Writer(sc);
		this.scenario = scenario;
		this.runtimesExact = new ArrayList<Long>();
		this.runtimesHeuristic = new ArrayList<Long>();
	}
	
	/**
	 * Creates an initial fuzzy knowledge base for each inventory-material by applying the Wang-Mendel method.
	 * @throws Exception 
	 */
	public void createInitialKnowledgeBases(boolean multiEchelon) throws Exception {
		long startTime = System.nanoTime();
		Map<String,KnowledgeBase> knowledgeBases = new HashMap<String,KnowledgeBase>();
	
		// Apply the Wang-Mendel method to each inventory and material
		for(Inventory i : inventories) {
			// Only consider PLANT-, WH- and DIST-type inventories
			if(i.isActiveInv()) {
				List<Material> materialList = materials.get(i);
				for(Material m : materialList) {
					WMGeneration wmObject = new WMGeneration(sc, i, m, numLabels);
					wmObject.createDataBase();
					wmObject.createRuleBase();
					wmObject.createKnowledgeBase(multiEchelon);
					KnowledgeBase initKb = wmObject.getWmKnowledgeBase();
					String tuple = "" + i.getInventoryId() + "," + m.getMaterialCode();
					knowledgeBases.put(tuple, initKb);
				}
			}
		}
		// Create initial knowledge base set
		this.initKbSet = new KnowledgeBaseSet(sc, knowledgeBases, gammaFit, phiFit);
		// This knowledge base set is currently the best
		this.bestKbSet = initKbSet;
		this.bestKbSetWM = initKbSet;
		writer.writeFinalResultsTraining(bestKbSet, null, 3, scenario);
		long runTime = System.nanoTime() - startTime;
		// Print the runtime of the WM method
		System.out.println("Runtime initialization exact method: " + runTime / (Math.pow(10, 9)) + " seconds");	
	}
	
	/**
	 * Applies Simulated Annealing in the global FRBS.
	 * @throws Exception 
	 */
	public void applySimulatedAnnealing() throws Exception {
		// Compute number of switches
		int totalNumRules = 0;
		for (var entry : bestKbSet.getKnowledgeBases().entrySet()) {
		    totalNumRules += entry.getValue().getRuleBase().size();
		}
		int numSwitchesPerRun = (int) (totalNumRules * percentageSwitches);
		
		// Create, initialize and run the Simulated Annealing Method
		GlobalSA gsa = new GlobalSA(sc, numLabels, bestKbSet, alpha, temp, numSwitchesPerRun, stoppingCount, gammaFit, phiFit);
		gsa.initialize();
		gsa.runSimulatedAnnealing();
		
		// Update optimal global knowledge base
		if(gsa.getBestKnowledgeBaseSet().getFitnessValue() > bestKbSet.getFitnessValue()) {
			this.bestKbSet = gsa.getBestKnowledgeBaseSet();
		}
		this.bestKbSetSA = gsa.getBestKnowledgeBaseSet();
	}
	
	/**
	 * Applies LA-Tuning (Genetic) algorithm in the global FRBS.
	 * @param simple Indicates whether rule reduction is included (false) in the algorithm or not (true)
	 * @throws Exception
	 */
	public void applyGeneticTuning(boolean simple) throws Exception {
		// Create, initialize and run the Genetic Tuning algorithm
		GlobalGTSimple ggt;
		if(simple) {
			ggt = new GlobalGTSimple(sc, bestKbSet, numLabels, N, bitsGene, phi, gammaFit, phiFit);
			ggt.runGeneticTuning(numCycles,maxIterations);
		} else {
			ggt = new GlobalGTSimple(sc, bestKbSet, numLabels, N, bitsGene, phi, gammaFit, phiFit);
			ggt.runGeneticTuning(numCycles,maxIterations);
		}
		
		// Update optimal global knowledge base
		this.bestKbSet = ggt.getBestKnowledgeBaseSet();
		this.bestKbSetLA = ggt.getBestKnowledgeBaseSet();
		writer.writeLAIterationValuesGlobal(ggt.getFitnessValues(), ggt.getObjValues(), ggt.getFillRates(), scenario);
	}
	
	/** Trains global multi-echelon FRBS.
	 * @param multiEchelon Must be true for this method
	 * @throws Exception
	 */
	public void trainFisModel(boolean multiEchelon) throws Exception {
		System.out.println("");
		System.out.println("Exact model:");
		System.out.println("");
		long startTime = System.nanoTime();
		// 1) Apply Wang-Mendel model
		this.createInitialKnowledgeBases(multiEchelon);
		// Print current results
		System.out.println("");
		System.out.println("Results after Wagner-Mendel: ");
		System.out.println("Fitness value: " + bestKbSet.getFitnessValue());
		System.out.println("Total SC costs: " + Math.round(bestKbSet.getObjVal() * Math.pow(10, -6) * 100.0) / 100.0 + " M");
		System.out.println("SC fill rate: " + Math.round(bestKbSet.getFillRate() * 100.0) / 100.0);
		System.out.println("Holding costs: " + bestKbSet.getSimulation().getTotalHoldingCosts());
		
		// 2) Apply Simulated Annealing algorithm
		this.applySimulatedAnnealing();
		// Print current results
		System.out.println("");
		System.out.println("Results after Simulated Annealing: ");
		System.out.println("Fitness value: " + bestKbSet.getFitnessValue());
		System.out.println("Total SC costs: " + Math.round(bestKbSet.getObjVal() * Math.pow(10, -6) * 100.0) / 100.0 + " M");
		System.out.println("SC fill rate: " + Math.round(bestKbSet.getFillRate() * 100.0) / 100.0);
		System.out.println("Holding costs: " + bestKbSet.getSimulation().getTotalHoldingCosts());
		System.out.println("");
		
		// 3) Apply LA-Tuning (Genetic tuning)
		this.applyGeneticTuning(true);
		System.out.println("Results after Genetic Tuning: ");
		System.out.println("Fitness value: " + bestKbSet.getFitnessValue());
		System.out.println("Total SC costs: " + Math.round(bestKbSet.getObjVal() * Math.pow(10, -6) * 100.0) / 100.0 + " M");
		System.out.println("SC fill rate: " + Math.round(bestKbSet.getFillRate() * 100.0) / 100.0);
		System.out.println("Holding costs: " + bestKbSet.getSimulation().getTotalHoldingCosts());
		
		// Write .txt files which contain info about the data bases/rule bases, final overall results and the interim results after applying the separate steps of the algorithm
		writer.writeDBsAndRBs(false, bestKbSet, scenario);
		writer.writeFinalResultsTraining(bestKbSet, null, 1, scenario);
		writer.writeInterimTrainingResultsExact(bestKbSetWM, bestKbSetSA, bestKbSetLA, scenario);
		long runTime = System.nanoTime() - startTime;
		runtimesExact.add((long) (runTime / (Math.pow(10, 9))));
		// Write .txt file with the run time
		writer.writeRuntime(runTime, false, scenario);
		// Print run time
		System.out.println("Runtime Exact model: " + runTime / (Math.pow(10, 9)) + " seconds");	
	}
	
	/** Applies the 3-step training procedure (WM, SA, LA-Tuning) to a particular inventory and material.
	 * @param inventory 
	 * @param material
	 * @param numLabels Number of fuzzy sets used for representing a fuzzy variable
	 * @param simple Indicates whether rule reduction is included (false) in the algorithm or not (true)
	 * @throws Exception
	 */
	public void trainFisModelLocal(Inventory inventory, Material material, int numLabels) throws Exception {
		long startTime = System.nanoTime();
		String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();	
		// Local Wang-Mendel algorithm
		WMGeneration wm = new WMGeneration(sc, inventory, material, numLabels);
		wm.createDataBase();
		wm.createRuleBase();
		wm.createKnowledgeBase(false);
		KnowledgeBase initKb = wm.getWmKnowledgeBase();
		// Print current results
		System.out.println("Results after Wagner-Mendel: ");
		System.out.println("Fitness value: " + initKb.getFitnessValue());
		System.out.println("Total SC costs: " + Math.round(initKb.getObjVal() * Math.pow(10, -6) * 100.0) / 100.0 + " M");
		System.out.println("SC fill rate: " + Math.round(initKb.getFillRate() * 100.0) / 100.0);
		System.out.println("Holding costs: " + initKb.getSimulation().getTotalHoldingCosts());
		
		// Local Simulated Annealing
		int numSwitches = (int) wm.getRuleBase().size() / 4;
		SimulatedAnnealing sa = new SimulatedAnnealing(sc, 3, inventory, material, initKb, alpha, temp, numSwitches, stoppingCount, gammaFit, phiFit);
		sa.initialize();
		sa.runSimulatedAnnealing();
		// Print current results
		System.out.println("");
		System.out.println("Results after Simulated Annealing: ");
		System.out.println("Fitness value: " + sa.getBestKnowledgeBase().getFitnessValue());
		System.out.println("Total SC costs: " + Math.round(sa.getBestKnowledgeBase().getObjVal() * Math.pow(10, -6) * 100.0) / 100.0 + " M");
		System.out.println("SC fill rate: " + Math.round(sa.getBestKnowledgeBase().getFillRate() * 100.0) / 100.0);
		System.out.println("Holding costs: " + sa.getBestKnowledgeBase().getSimulation().getTotalHoldingCosts());
		
		List<String> labels = sc.getLabels3();
		if(numLabels == 5) {
			labels = sc.getLabels5();
		}
		// Local Genetic Tuning
		GeneticTuningSimple gt = new GeneticTuningSimple(sc, sa.getBestKnowledgeBase(), 3, N, bitsGene, phi, gammaFit, phiFit);
		gt.runGeneticTuning(numCycles,maxIterations);
		System.out.println("");
		System.out.println("Results after Genetic Tuning: ");
		System.out.println("Fitness value: " + gt.getBestKnowledgeBase().getFitnessValue());
		System.out.println("Total SC costs: " + Math.round(gt.getBestKnowledgeBase().getObjVal() * Math.pow(10, -6) * 1000.0) / 1000.0 + " M");
		System.out.println("SC fill rate: " + Math.round(gt.getBestKnowledgeBase().getFillRate() * 100.0) / 100.0);
		System.out.println("Holding costs: " + gt.getBestKnowledgeBase().getSimulation().getTotalHoldingCosts());
		bestKnowledgeBasesHeuristic.put(tuple, gt.getBestKnowledgeBase());
		// Write interim results into .txt file
		writer.writeInterimTrainingResultsLocal(tuple, wm.getWmKnowledgeBase(), sa.getBestKnowledgeBase(), gt.getBestKnowledgeBase(), scenario);
		// Write LA-Tuning fitness per iteration into .txt file
		writer.writeLAIterationValuesLocal(tuple, gt.getFitnessValues(), gt.getObjValues(), gt.getFillRates(), scenario);
		long runTime = System.nanoTime() - startTime;
		// Print run time
		System.out.println("Runtime tuple " + tuple + ": " + runTime / (Math.pow(10, 9)) + " seconds");
	}
	
	/**
	 * Trains global multi-echelon FRBS Heuristic (applies 3-step procedure to each inventory-material separately).
	 * @throws Exception 
	 */
	public void trainFISModelHeuristic() throws Exception {
		System.out.println("");
		System.out.println("Heuristic:");
		System.out.println("");
		long startTime = System.nanoTime();
		// Create one new quantity vector per inventory-material of previous echelon (Warehouses)
		HashMap<String,TreeMap<Date,Integer>> newQntVectors = new HashMap<String,TreeMap<Date,Integer>>();
		// Initialize new quantity vectors of previous echelon
		for(Inventory warehouse : warehouses) {
			for(Material material : warehouse.getMaterials()) {
				String tuple = "" + warehouse.getInventoryId() + "," + material.getMaterialCode();
				TreeMap<Date,Integer> emptyMap = new TreeMap<Date,Integer>();
				for(Date date : sc.getDates()) {
					if(!date.after(sc.getSplitDate())) {
						emptyMap.put(date, 0);
					} else {
						emptyMap.put(date, sc.getDemandsMap().get(tuple).get(date));
					}
				}
				newQntVectors.put(tuple, emptyMap);
			}
		}
		
		// 1) Simulate the inventories of distribution centers
		for(Inventory distributionCenter : distributionCenters) {
			for(Material material : distributionCenter.getMaterials()) {
				String tuple = "" + distributionCenter.getInventoryId() + "," + material.getMaterialCode();
				// Train local FRBS
				this.trainFisModelLocal(distributionCenter, material, numLabels);
				System.out.println("");
				System.out.println("Inventory/Material " + distributionCenter.getInventoryId() + "," + material.getMaterialCode() + " completed.");
				System.out.println("");
				KnowledgeBase kb = bestKnowledgeBasesHeuristic.get(tuple);
				Simulation simulation = kb.getSimulation();
				// Retrieve order quantities of whole training horizon
				TreeMap<Date,Integer> orderQuantities = simulation.getOrderQuantitiesFIS().get(tuple);
				// Retrieve predecessor
				String predecessor = "" + distributionCenter.getPredecessor() + "," + material.getMaterialCode();
				Inventory predecessorInv = this.getInventoryByCode(distributionCenter.getPredecessor());
				int indexSuccessor = predecessorInv.getSuccessors().indexOf(distributionCenter.getInventoryId());
				String predecessorTuple = "" + distributionCenter.getPredecessor() + "," + material.getMaterialCode();
				TreeMap<Date,Integer> orderQuantitiesCopy = orderQuantities;
				// Replace sub-demands of predecessor (to this inventory) by new order quantities
				sc.getScSubDemands().get(predecessorTuple).set(indexSuccessor, orderQuantitiesCopy);
				// Add new order quantities to new quantity vector of predecessor
				for(Date trainingDate : sc.getTrainingDates()) {
					newQntVectors.get(predecessor).put(trainingDate, newQntVectors.get(predecessor).get(trainingDate) + orderQuantities.get(trainingDate));
				}
			}
		}
		// Replace total demands of predecessors by new quantity vectors
		for (var entry : newQntVectors.entrySet()) {
			String tuple = entry.getKey();
			TreeMap<Date,Integer> mapToAdd = new TreeMap<Date,Integer>(entry.getValue());
		    sc.getDemandsMap().put(tuple, mapToAdd);
		}
		// Empty new quantity vectors
		newQntVectors.clear();
		
		// 2) Simulate the inventories of warehouses 
		// Initialize new quantity vectors of previous echelon
		for(Inventory plant : plants) {
			for(Material material : plant.getMaterials()) {
				String tuple = "" + plant.getInventoryId() + "," + material.getMaterialCode();
				TreeMap<Date,Integer> emptyMap = new TreeMap<Date,Integer>();
				for(Date date : sc.getDates()) {
					if(!date.after(sc.getSplitDate())) {
						emptyMap.put(date, 0);
					} else {
						emptyMap.put(date, sc.getDemandsMap().get(tuple).get(date));
					}
				}
				newQntVectors.put(tuple, emptyMap);
			}
		}
		
		for(Inventory warehouse : warehouses) {
			for(Material material : warehouse.getMaterials()) {
				String tuple = "" + warehouse.getInventoryId() + "," + material.getMaterialCode();
				// Train local FRBS
				this.trainFisModelLocal(warehouse, material, numLabels);
				System.out.println("");
				System.out.println("Inventory/Material " + warehouse.getInventoryId() + "," + material.getMaterialCode() + " completed.");
				System.out.println("");
				KnowledgeBase kb = bestKnowledgeBasesHeuristic.get(tuple);
				Simulation simulation = kb.getSimulation();
				// Retrieve order quantities of whole training horizon
				TreeMap<Date,Integer> orderQuantities = simulation.getOrderQuantitiesFIS().get(tuple);
				// NOTE: demands of plants are expressed in terms of raw materials
				for(Material matROH : material.getIngredients().keySet()) {
					// Retrieve raw material predecessor
					String predecessor = "" + warehouse.getPredecessor() + "," + matROH.getMaterialCode();
					for(Date trainingDate : sc.getTrainingDates()) {
						// Retrieve the order quantity of this particular raw material
						double amount;
						if(matROH.getMaterialCode().equals("CC-P01") || matROH.getMaterialCode().equals("CC-P02") || matROH.getMaterialCode().equals("CC-P03") || matROH.getMaterialCode().equals("CC-P04")) {
							amount = orderQuantities.get(trainingDate) * material.getIngredients().get(matROH);
						} else {
							amount = orderQuantities.get(trainingDate) * material.getIngredients().get(matROH) * material.getUnitSize();
						}	
						int qntToAdd = (int) Math.ceil(amount);
						// Add new order quantities to new quantity vector of raw material predecessor
						newQntVectors.get(predecessor).put(trainingDate, newQntVectors.get(predecessor).get(trainingDate) + qntToAdd);
					}
				}
			}
		}
		
		// Replace sub-demands and total demands of raw material predecessors by new quantity vectors (sub-demand vector is equal to total demand vector)
		// NOTE: demands of plants are expressed in terms of raw materials
		for (var entry : newQntVectors.entrySet()) {
			String tuple = entry.getKey();
			TreeMap<Date,Integer> mapToAdd = new TreeMap<Date,Integer>(entry.getValue());
			TreeMap<Date,Integer> subDemandsToAdd = new TreeMap<Date,Integer>(entry.getValue());
		    sc.getDemandsMap().put(tuple, mapToAdd);
		    sc.getScSubDemands().get(tuple).set(0, subDemandsToAdd);
		}
		
		// 3) Simulate the inventories of plants
		for(Inventory plant : plants) {
			for(Material material : plant.getMaterials()) {
				this.trainFisModelLocal(plant, material, numLabels);
				System.out.println("");
				System.out.println("Inventory/Material " + plant.getInventoryId() + "," + material.getMaterialCode() + " completed.");
				System.out.println("");
			}
		}
		
		// Reset input data sets before applying FRBS to entire SC
		sc.getDataGenerator().addDataToSc();
		
		// Apply local fuzzy models to entire supply chain
		// Create knowledge base set  
		bestKbSet = new KnowledgeBaseSet(sc, this.bestKnowledgeBasesHeuristic, gammaFit, phiFit);
		// Write result .txt files
		writer.writeDBsAndRBs(true, bestKbSet, scenario);
		writer.writeFinalResultsTraining(bestKbSet, null, 2, scenario);
		writer.writeInventoryLevelsTraining(bestKbSet.getSimulation(), null, 2);
		writer.writeInventoryPositionsTraining(bestKbSet.getSimulation(), null, 2);
		// Print training results of heuristic
		System.out.println("Results after HEURISTIC: ");
		System.out.println("Fitness value: " + bestKbSet.getFitnessValue());
		System.out.println("Total SC costs: " + Math.round(bestKbSet.getObjVal() * Math.pow(10, -6) * 100.0) / 100.0 + " M");
		System.out.println("SC fill rate: " + Math.round(bestKbSet.getFillRate() * 100.0) / 100.0);
		System.out.println("Holding costs: " + bestKbSet.getSimulation().getTotalHoldingCosts());
		long runTime = System.nanoTime() - startTime;
		runtimesHeuristic.add((long) (runTime / (Math.pow(10, 9))));
		// Write run time into .txt file
		writer.writeRuntime(runTime, true, scenario);
		// Print run time
		System.out.println("Runtime Heuristic: " + runTime / (Math.pow(10, 9)) + " seconds");	
	}
	
	/** Tests FRBS. A FRBS (exact or heuristic) is applied to test data.
	 * @param heuristic Is true if the heuristic is tested and false if the exact model is tested
	 * @throws Exception
	 */
	public void testModel(boolean heuristic) throws Exception {
		// Create .fcl file of each inventory-material --> contains info about fuzzy model (DB and RB)
		for(Inventory i : sc.getSortedInventories()) {
			for(Material m : i.getMaterials()) {
				String tuple = "" + i.getInventoryId() + "," + m.getMaterialCode();
				KnowledgeBase kb =  bestKbSet.getKnowledgeBases().get(tuple);
				Writer kbWriter = new Writer(sc);
				double[] dataBase = kb.getDataBase();
				HashSet<String[]> ruleBase = kb.getRuleBase();
				if(i.getStageType().equals("PLANT")) {
					kbWriter.writeFCL(i, m, dataBase, ruleBase, sc.getVarNamesL(), sc.getLabels3(), 0);
				} else if(i.getStageType().equals("WH")) {
					kbWriter.writeFCL(i, m, dataBase, ruleBase, sc.getVarNamesM(), sc.getLabels3(), 0);
				} else {
					kbWriter.writeFCL(i, m, dataBase, ruleBase, sc.getVarNamesS(), sc.getLabels3(), 0);
				}
			}
		}
		// Simulate entire SC during test phase by applying the fuzzy models obtained from the KB generation method
		Simulation simulationFIS = new Simulation(sc,0);
		simulationFIS.simulateSupplyChain(false);
		int modelType;
		if(heuristic) {
			modelType = 2;
		} else {
			modelType = 1;
		}
		// Write final results into .txt file
		writer.writeFinalResultsTest(simulationFIS, null, modelType, scenario);
		if(modelType == 2) {
			writer.writeInventoryLevelsTest(simulationFIS, null, modelType);
			writer.writeInventoryPositionsTest(simulationFIS, null, modelType);
		}
		this.totalCosts = simulationFIS.getTotalCosts();
		this.totalScCosts = simulationFIS.getTotalScCosts();
		this.scFillRate = simulationFIS.getSCFillRate();
		// Print results for each inventory and material
		System.out.println("");
		// Print SC objective value and fill rate
		System.out.println("Results of testing Fuzzy Model: ");
		System.out.println("Total SC costs: " + Math.round(totalScCosts * Math.pow(10, -6) * 100.0) / 100.0 + " M");
		System.out.println("SC fill rate: " + Math.round(scFillRate * 100.0) / 100.0);
		System.out.println("Holding costs: " + simulationFIS.getTotalHoldingCosts());
		System.out.println("");
		// Delete .fcl files
		this.deleteFclFiles();
	}

	/**
	 * @return the total costs of the SC
	 */
	public double getTotalScCosts() {
		for (var entry : bestKbSet.getKnowledgeBases().entrySet()) {
			totalScCosts += entry.getValue().getObjVal();
		}
		return totalScCosts;
	}
	
	/**
	 * @return the sc
	 */
	public SC getSc() {
		return sc;
	}

	/**
	 * @return the inventories
	 */
	public List<Inventory> getInventories() {
		return inventories;
	}

	/**
	 * @return the materials
	 */
	public Map<Inventory, List<Material>> getMaterials() {
		return materials;
	}

	/**
	 * @return the numLabels
	 */
	public int getNumLabels() {
		return numLabels;
	}

	/**
	 * @return the wmObjects
	 */
	public List<WMGeneration> getWmObjects() {
		return wmObjects;
	}
		
	/**
	 * @return the runtimesExact
	 */
	public List<Long> getRuntimesExact() {
		return runtimesExact;
	}

	/**
	 * @return the runtimesHeuristic
	 */
	public List<Long> getRuntimesHeuristic() {
		return runtimesHeuristic;
	}

	/**
	 * Deletes all .fcl files in the directory.
	 */
	public void deleteFclFiles() {
		File directory = new File("C:\\Users\\iliad\\Desktop\\Eclipse_workspace\\KB_Generation");
		File[] files = directory.listFiles();
		for (File f : files)
		{
		    if (f.getName().startsWith("kbGeneration"))
		    {
		      f.delete();
		    }
		}
	}
	
    /** Returns the inventory object corresponding to an inventoryId.
     * @param inventoryId
     * @return
     */
    public Inventory getInventoryByCode(String inventoryId) {
    	Inventory inventory = null;
    	for(Inventory i : inventories) {
    		if(i.getInventoryId().equals(inventoryId)) {
    			inventory = i;
    			break;
    		}
    	}
    	return inventory;
    }
}