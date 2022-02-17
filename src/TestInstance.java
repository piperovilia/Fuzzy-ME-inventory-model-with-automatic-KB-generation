import java.io.File;
import java.io.FileNotFoundException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 * A class to compute the average test performance of the three models: EOQ, exact FRBS, heuristic FRBS
 */
public class TestInstance {
	private SC sc;
	private List<Map<String,KnowledgeBase>> knowledgeBasesExact; // List containing KB sets resulting from exact FRBS for each scenario
	private List<Map<String,KnowledgeBase>> knowledgeBasesHeuristic; // List containing KB sets resulting from heuristic FRBS for each scenario
	private List<Map<String,TreeMap<Date,Integer>>> trainingDemand; // List containing training demands of each scenario (relevant for EOQ)
	private List<Map<String,TreeMap<Date,Integer>>> trainingLeadTime; // List containing training lead times of each scenario (relevant for EOQ)
	private boolean multiEchelon;
	private double gammaFit; 
	private double phiFit;
	private List<List<Integer>> demandLists;
	private List<List<Double>> cdfMapsDemand;
	private List<List<Double>> demandParameters;
	private List<List<Double>> leadTimeParameters;
	private Map<Material,Double> meanPrices;
	private List<Map<Material,Double>> marketPriceStdvs;
	
	// Lists used for tracking results (per method)
	// Total SC costs
	private List<List<Double>> totalScCostsEOQ;
	private List<List<Double>> totalScCostsExact;
	private List<List<Double>> totalScCostsHeuristic;
	// Fill rates
	private List<List<Double>> scFillRateEOQ;
	private List<List<Double>> scFillRateExact;
	private List<List<Double>> scFillRateHeuristic;
	// Fitness values
	private List<List<Double>> fitnessValueEOQ;
	private List<List<Double>> fitnessValueExact;
	private List<List<Double>> fitnessValueHeuristic;
	// Holding costs
	private List<List<Double>> holdingCostsEOQ;
	private List<List<Double>> holdingCostsExact;
	private List<List<Double>> holdingCostsHeuristic;
	// Penalty costs
	private List<List<Double>> penaltyCostsEOQ;
	private List<List<Double>> penaltyCostsExact;
	private List<List<Double>> penaltyCostsHeuristic;
	// Order costs
	private List<List<Double>> orderCostsEOQ;
	private List<List<Double>> orderCostsExact;
	private List<List<Double>> orderCostsHeuristic;
	// Setup costs
	private List<List<Double>> setupCostsEOQ;
	private List<List<Double>> setupCostsExact;
	private List<List<Double>> setupCostsHeuristic;
	// Production costs
	private List<List<Double>> productionCostsEOQ;
	private List<List<Double>> productionCostsExact;
	private List<List<Double>> productionCostsHeuristic;
	// Purchase costs
	private List<List<Double>> purchaseCostsEOQ;
	private List<List<Double>> purchaseCostsExact;
	private List<List<Double>> purchaseCostsHeuristic;
	// Transport costs
	private List<List<Double>> transportCostsEOQ;
	private List<List<Double>> transportCostsExact;
	private List<List<Double>> transportCostsHeuristic;
		
	/**
	 * Default constructor.
	 * @param sc SC object
	 * @param gammaFit Parameter to prioritize SC costs
	 * @param phiFit Parameter to prioritize fill rate
	 * @param demandLists Contains demand quantities if a custom discrete demand distribution is applied
	 * @param cdfMapsDemand Contains CDFs of demands if custom discrete demand distribution is applied
	 * @param leadTimeParameters Parameters of lead time distributions
	 * @param meanPrices Mean raw material prices for Normal distributions
	 * @param marketPriceStdvs Standard deviations of raw material prices for Normal distribution
	 */
	public TestInstance(SC sc, double gammaFit, double phiFit, List<List<Integer>> demandLists, List<List<Double>> cdfMapsDemand, List<List<Double>> leadTimeParameters,
			List<List<Double>> demandParameters, Map<Material,Double> meanPrices, List<Map<Material,Double>> marketPriceStdvs) {
		this.sc = sc;
		this.knowledgeBasesExact = new ArrayList<Map<String,KnowledgeBase>>();
		this.knowledgeBasesHeuristic = new ArrayList<Map<String,KnowledgeBase>>();
		this.trainingDemand = new ArrayList<Map<String,TreeMap<Date,Integer>>>();
		this.trainingLeadTime = new ArrayList<Map<String,TreeMap<Date,Integer>>>();
		for(int i=0; i<8; i++) {
			Map<String,KnowledgeBase> emptyMapExact = new HashMap<>();
			Map<String,KnowledgeBase> emptyMapHeuristic = new HashMap<>();
			Map<String,TreeMap<Date,Integer>> emptyMapDemands = new HashMap<>();
			Map<String,TreeMap<Date,Integer>> emptyMapLeadTimes = new HashMap<>();
			this.knowledgeBasesExact.add(emptyMapExact);
			this.knowledgeBasesHeuristic.add(emptyMapHeuristic);
			this.trainingDemand.add(emptyMapDemands);
			this.trainingLeadTime.add(emptyMapLeadTimes);
		}
		this.multiEchelon = true;
		this.gammaFit = gammaFit;
		this.phiFit = phiFit;
		this.demandLists = demandLists;
		this.cdfMapsDemand = cdfMapsDemand;
		this.leadTimeParameters = leadTimeParameters;
		this.demandParameters = demandParameters;
		this.meanPrices = meanPrices;
		this.marketPriceStdvs = marketPriceStdvs;
		
		// Initialize tracking lists
		this.totalScCostsEOQ = new ArrayList<List<Double>>();
		this.initializeEmptyList(totalScCostsEOQ);
		this.totalScCostsExact = new ArrayList<List<Double>>();
		this.initializeEmptyList(totalScCostsExact);
		this.totalScCostsHeuristic = new ArrayList<List<Double>>();
		this.initializeEmptyList(totalScCostsHeuristic);
		
		this.scFillRateEOQ = new ArrayList<List<Double>>();
		this.initializeEmptyList(scFillRateEOQ);
		this.scFillRateExact = new ArrayList<List<Double>>();
		this.initializeEmptyList(scFillRateExact);
		this.scFillRateHeuristic = new ArrayList<List<Double>>();
		this.initializeEmptyList(scFillRateHeuristic);
		
		this.fitnessValueEOQ = new ArrayList<List<Double>>();
		this.initializeEmptyList(fitnessValueEOQ);
		this.fitnessValueExact = new ArrayList<List<Double>>();
		this.initializeEmptyList(fitnessValueExact);
		this.fitnessValueHeuristic = new ArrayList<List<Double>>();
		this.initializeEmptyList(fitnessValueHeuristic);
		
		this.holdingCostsEOQ = new ArrayList<List<Double>>();
		this.initializeEmptyList(holdingCostsEOQ);
		this.holdingCostsExact = new ArrayList<List<Double>>();
		this.initializeEmptyList(holdingCostsExact);
		this.holdingCostsHeuristic = new ArrayList<List<Double>>();
		this.initializeEmptyList(holdingCostsHeuristic);
		
		this.penaltyCostsEOQ = new ArrayList<List<Double>>();
		this.initializeEmptyList(penaltyCostsEOQ);
		this.penaltyCostsExact = new ArrayList<List<Double>>();
		this.initializeEmptyList(penaltyCostsExact);
		this.penaltyCostsHeuristic = new ArrayList<List<Double>>();
		this.initializeEmptyList(penaltyCostsHeuristic);
		
		this.orderCostsEOQ = new ArrayList<List<Double>>();
		this.initializeEmptyList(orderCostsEOQ);
		this.orderCostsExact = new ArrayList<List<Double>>();
		this.initializeEmptyList(orderCostsExact);
		this.orderCostsHeuristic = new ArrayList<List<Double>>();
		this.initializeEmptyList(orderCostsHeuristic);
		
		this.setupCostsEOQ = new ArrayList<List<Double>>();
		this.initializeEmptyList(setupCostsEOQ);
		this.setupCostsExact = new ArrayList<List<Double>>();
		this.initializeEmptyList(setupCostsExact);
		this.setupCostsHeuristic = new ArrayList<List<Double>>();
		this.initializeEmptyList(setupCostsHeuristic);
		
		this.productionCostsEOQ = new ArrayList<List<Double>>();
		this.initializeEmptyList(productionCostsEOQ);
		this.productionCostsExact = new ArrayList<List<Double>>();
		this.initializeEmptyList(productionCostsExact);
		this.productionCostsHeuristic = new ArrayList<List<Double>>();
		this.initializeEmptyList(productionCostsHeuristic);
		
		this.purchaseCostsEOQ = new ArrayList<List<Double>>();
		this.initializeEmptyList(purchaseCostsEOQ);
		this.purchaseCostsExact = new ArrayList<List<Double>>();
		this.initializeEmptyList(purchaseCostsExact);
		this.purchaseCostsHeuristic = new ArrayList<List<Double>>();
		this.initializeEmptyList(purchaseCostsHeuristic);
		
		this.transportCostsEOQ = new ArrayList<List<Double>>();
		this.initializeEmptyList(transportCostsEOQ);
		this.transportCostsExact = new ArrayList<List<Double>>();
		this.initializeEmptyList(transportCostsExact);
		this.transportCostsHeuristic = new ArrayList<List<Double>>();
		this.initializeEmptyList(transportCostsHeuristic);	
	}
	
	/**
	 * Initialize a given list with 8 (number of scenarios) empty lists of doubles.
	 * @param list List to be initialized
	 */
	public void initializeEmptyList(List<List<Double>> list) {
		for(int i=0; i<8; i++) {
			List<Double> emptyList = new ArrayList<Double>();
			list.add(emptyList);
		}
	}

	/**
	 * Recreates KBs (DB + RB), demands and lead times from training (for each scenario). The method extracts the relevant information from .txt files. 
	 * @param heuristic Is true if method is applied to the heuristic and false if applied to the exact model
	 * @throws Exception
	 */
	public void createModelInputs(boolean heuristic) throws Exception {
		String modelType;
		if(heuristic) {
			modelType = "Heuristic";
		} else {
			modelType = "Exact";
		}
		for(int i=1; i<=8; i++) {
			for(Inventory inventory : sc.getSortedInventories()) {
				for(Material material : inventory.getMaterials()) {
					// Create DB
					String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
					List<Double> dataBaseList = new ArrayList<Double>();
					String fileStrDb = "dataBase" + tuple + modelType + "Scenario" + i + ".txt";
					File fileDataBase = new File(fileStrDb);					
			        try (Scanner s = new Scanner(fileDataBase)) {
			        	s.nextLine(); // Skip line with field names
			        	while(s.hasNext()) {
			            	String valueStr = s.next();
			            	double value = Double.parseDouble(valueStr);
			            	dataBaseList.add(value);
			        	}   
			        }
					double[] dataBase = new double[inventory.getNumVars() * 3 * 3];
					for(int j=0; j<dataBaseList.size(); j++) {
						dataBase[j] = dataBaseList.get(j);
					}
					
					// Create RB
					String fileStrRb = "ruleBase" + tuple + modelType + "Scenario" + i + ".txt";
					File fileRuleBase = new File(fileStrRb);
					int numVars = inventory.getNumVars();
					HashSet<String[]> ruleBase = new HashSet<String[]>();
					
			        try (Scanner s = new Scanner(fileRuleBase)) {
			        	int counter = 0;
			        	String[] rule = new String[numVars];
			        	s.nextLine(); // Skip line with field names
			        	while(s.hasNext()) {
			            	String label = s.next();
			            	rule[counter] = label;
			            	if(counter == numVars - 1) {
			            		counter = 0;
			            		ruleBase.add(rule);
			            		rule = new String[numVars];
			            	} else {
			            		counter++;
			            	}
			        	}   
			        }
			        
			        // Create KB by using DB and RB from above
			        KnowledgeBase kb = new KnowledgeBase(sc, inventory, material, dataBase, ruleBase, multiEchelon, gammaFit, phiFit);
			        if(heuristic) {
			        	knowledgeBasesHeuristic.get(i-1).put(tuple, kb);
			        } else {
			        	knowledgeBasesExact.get(i-1).put(tuple, kb);
			        }
			        
			        // Create demand and lead time lists based on input data files 
			        String fileStrInput = "inputData" + tuple + "," + i + ".txt";
			        File inputData = new File(fileStrInput);
			        TreeMap<Date,Integer> demands = new TreeMap<Date,Integer>();
			        TreeMap<Date,Integer> leadTimes = new TreeMap<Date,Integer>();
			        try (Scanner s = new Scanner(inputData)) {
			        	s.nextLine(); // Skip line with field names
			        	while(s.hasNext()) {
			        		String dateStr = s.next();
			        		Date date = new SimpleDateFormat("yyyy-MM-dd").parse(dateStr);
			            	int demand = s.nextInt();
			            	s.next();
			            	s.next();
			            	int leadTime;
			            	if(inventory.getStageType().equals("DIST")) {
			            		leadTime = 24;
			            	} else {
			            		leadTime = s.nextInt();
			            	}
			            	demands.put(date, demand);
			            	leadTimes.put(date, leadTime);
			            	s.nextLine();
			        	}   
			        }
			        trainingDemand.get(i-1).put(tuple, demands);
			        trainingLeadTime.get(i-1).put(tuple, leadTimes);
				}
			}
		}
	}
	
	/**
	 * Runs 'iterations' tests for each uncertainty scenario.
	 * @param iterations Number of test runs
	 * @param numDates Number of dates (training + test)
	 * @param startDate Start date of time horizon
	 * @param splitRatio Percentage of dates used for training --> [0, 1]
	 * @param isDemand Is true if demand is used and false if demand changes
	 * @throws Exception
	 */
	public void runTest(int iterations, int numDates, Date startDate, double splitRatio, boolean isDemand, boolean continuousDemand) throws Exception {
		// Create all model inputs for the exact FRBS and the heuristic FRBS
		this.createModelInputs(false);
		this.createModelInputs(true);
		int counter = 0; // Counter to keep track of scenario
		// Run all scenarios(2^3 scenarios)
		for(int i=0; i<2; i++) {
			for(int j=0; j<2; j++) {
				for(int k=0; k<2; k++) {
					for(int l=0; l<iterations; l++) {
						DataReader readerNew = new DataReader(splitRatio, isDemand);
						readerNew.readDataSimulation(startDate, numDates);
						SC scObject = readerNew.getSupplyChain();
						// DATA GENERATION
						List<Integer> demandsList = demandLists.get(i);
						List<Double> cdfMapDemand = cdfMapsDemand.get(i);
						double meanDemand = demandParameters.get(i).get(0);
						double stdvDemand = demandParameters.get(i).get(1);
						double shapeWh = leadTimeParameters.get(j).get(0);
						double scaleWh = leadTimeParameters.get(j).get(1);
						double shapeP = leadTimeParameters.get(j).get(0);
						double scaleP = leadTimeParameters.get(j).get(1);
						Map<Material,Double> priceStdvs  = marketPriceStdvs.get(k);
						boolean isSimulation = true;
						scObject.setIsSimulation(isSimulation);
						DataGeneration dataGenerator;
						if(continuousDemand) {
							dataGenerator = new DataGeneration(scObject, readerNew, scObject.getDates(), continuousDemand, meanDemand, stdvDemand, null, null, shapeWh, scaleWh, shapeP, scaleP, meanPrices, priceStdvs);
						} else {
							dataGenerator = new DataGeneration(scObject, readerNew, scObject.getDates(), continuousDemand, 0.0, 0.0, demandsList, cdfMapDemand, shapeWh, scaleWh, shapeP, scaleP, meanPrices, priceStdvs);
						}
						dataGenerator.generateSCData();
						
						// Simulate EOQ model
						Simulation simulationEOQ = new Simulation(scObject, 0);
						simulationEOQ.simulateSupplyChainEOQ(trainingDemand.get(counter),trainingLeadTime.get(counter));
						// Compute maximal SC costs
						double maxScCosts = 5 * simulationEOQ.getTotalScCosts();
						// Compute various KPIs
						this.totalScCostsEOQ.get(counter).add(simulationEOQ.getTotalScCosts());
						this.scFillRateEOQ.get(counter).add(simulationEOQ.getSCFillRate());
						this.fitnessValueEOQ.get(counter).add((double) ((1 - simulationEOQ.getTotalScCosts() / maxScCosts) * simulationEOQ.getSCFillRate()));
						this.holdingCostsEOQ.get(counter).add(simulationEOQ.getTotalHoldingCosts());
						this.penaltyCostsEOQ.get(counter).add(simulationEOQ.getTotalPenaltyCosts());
						this.orderCostsEOQ.get(counter).add(simulationEOQ.getTotalOrderCosts());
						this.setupCostsEOQ.get(counter).add(simulationEOQ.getTotalSetupCosts());
						this.productionCostsEOQ.get(counter).add(simulationEOQ.getTotalProductionCosts());
						this.purchaseCostsEOQ.get(counter).add(simulationEOQ.getTotalUnitCosts());
						this.transportCostsEOQ.get(counter).add(simulationEOQ.getTotalTransportationCosts());
						
						// Reset data in SC
						dataGenerator.addDataToSc();
						
						// Simulate exact FRBS
						// Write .fcl files
						for(Inventory inv : scObject.getSortedInventories()) {
							for(Material mat : inv.getMaterials()) {
								String tuple = "" + inv.getInventoryId() + "," + mat.getMaterialCode();
								KnowledgeBase kb =  knowledgeBasesExact.get(counter).get(tuple);
								
								// Apply FIS model obtained by training
								Writer kbWriter = new Writer(scObject);
								double[] dataBase = kb.getDataBase();
								HashSet<String[]> ruleBase = kb.getRuleBase();
								if(inv.getStageType().equals("PLANT")) {
									kbWriter.writeFCL(inv, mat, dataBase, ruleBase, scObject.getVarNamesL(), scObject.getLabels3(), 0);
								} else if(inv.getStageType().equals("WH")) {
									kbWriter.writeFCL(inv, mat, dataBase, ruleBase, scObject.getVarNamesM(), scObject.getLabels3(), 0);
								} else {
									kbWriter.writeFCL(inv, mat, dataBase, ruleBase, scObject.getVarNamesS(), scObject.getLabels3(), 0);
								}
							}
						}	
						// Run simulation
						Simulation simulationExact = new Simulation(scObject, 0);
						simulationExact.simulateSupplyChain(false);
						// Store KPIs
						this.totalScCostsExact.get(counter).add(simulationExact.getTotalScCosts());
						this.scFillRateExact.get(counter).add(simulationExact.getSCFillRate());
						this.fitnessValueExact.get(counter).add((double) ((1 - simulationExact.getTotalScCosts() / maxScCosts) * simulationExact.getSCFillRate()));
						this.holdingCostsExact.get(counter).add(simulationExact.getTotalHoldingCosts());
						this.penaltyCostsExact.get(counter).add(simulationExact.getTotalPenaltyCosts());
						this.orderCostsExact.get(counter).add(simulationExact.getTotalOrderCosts());
						this.setupCostsExact.get(counter).add(simulationExact.getTotalSetupCosts());
						this.productionCostsExact.get(counter).add(simulationExact.getTotalProductionCosts());
						this.purchaseCostsExact.get(counter).add(simulationExact.getTotalUnitCosts());
						this.transportCostsExact.get(counter).add(simulationExact.getTotalTransportationCosts());
						
						// Delete .fcl files
						this.deleteFclFiles();
						
						// Reset data in SC
						dataGenerator.addDataToSc();
						
						// Simulate heuristic FRBS
						for(Inventory inv : scObject.getSortedInventories()) {
							for(Material mat : inv.getMaterials()) {
								String tuple = "" + inv.getInventoryId() + "," + mat.getMaterialCode();
								KnowledgeBase kb =  knowledgeBasesHeuristic.get(counter).get(tuple);
								
								// Apply FIS model obtained by training
								Writer kbWriter = new Writer(scObject);
								double[] dataBase = kb.getDataBase();
								HashSet<String[]> ruleBase = kb.getRuleBase();
								if(inv.getStageType().equals("PLANT")) {
									kbWriter.writeFCL(inv, mat, dataBase, ruleBase, scObject.getVarNamesL(), scObject.getLabels3(), 0);
								} else if(inv.getStageType().equals("WH")) {
									kbWriter.writeFCL(inv, mat, dataBase, ruleBase, scObject.getVarNamesM(), scObject.getLabels3(), 0);
								} else {
									kbWriter.writeFCL(inv, mat, dataBase, ruleBase, scObject.getVarNamesS(), scObject.getLabels3(), 0);
								}
							}
						}
						// Run simulation
						Simulation simulationHeuristic = new Simulation(scObject,0);
						simulationHeuristic.simulateSupplyChain(false);
						// Store KPIs
						this.totalScCostsHeuristic.get(counter).add(simulationHeuristic.getTotalScCosts());
						this.scFillRateHeuristic.get(counter).add(simulationHeuristic.getSCFillRate());
						this.fitnessValueHeuristic.get(counter).add((double) ((1 - simulationHeuristic.getTotalScCosts() / maxScCosts) * simulationHeuristic.getSCFillRate()));
						this.holdingCostsHeuristic.get(counter).add(simulationHeuristic.getTotalHoldingCosts());
						this.penaltyCostsHeuristic.get(counter).add(simulationHeuristic.getTotalPenaltyCosts());
						this.orderCostsHeuristic.get(counter).add(simulationHeuristic.getTotalOrderCosts());
						this.setupCostsHeuristic.get(counter).add(simulationHeuristic.getTotalSetupCosts());
						this.productionCostsHeuristic.get(counter).add(simulationHeuristic.getTotalProductionCosts());
						this.purchaseCostsHeuristic.get(counter).add(simulationHeuristic.getTotalUnitCosts());
						this.transportCostsHeuristic.get(counter).add(simulationHeuristic.getTotalTransportationCosts());
						// Delete .fcl files
						this.deleteFclFiles();
					}
					
					System.out.println("Secnario: " + (counter + 1));
					
					// Average values of EOQ model
					System.out.println("EOQ: ");
					System.out.println("Average total costs: " + this.getAverageOfList(totalScCostsEOQ.get(counter)) + " | MIN: " + Collections.min(totalScCostsEOQ.get(counter)) + " | MAX: " + Collections.max(totalScCostsEOQ.get(counter)));
					System.out.println("Average fill rate: " + this.getAverageOfList(scFillRateEOQ.get(counter)) + " | MIN: " + Collections.min(scFillRateEOQ.get(counter)) + " | MAX: " + Collections.max(scFillRateEOQ.get(counter)));
					System.out.println("Average fitness value: " + this.getAverageOfList(fitnessValueEOQ.get(counter)) + " | MIN: " + Collections.min(fitnessValueEOQ.get(counter)) + " | MAX: " + Collections.max(fitnessValueEOQ.get(counter)));
					System.out.println("Average holding costs: " + this.getAverageOfList(holdingCostsEOQ.get(counter)));
					System.out.println("Average penalty costs: " + this.getAverageOfList(penaltyCostsEOQ.get(counter)));
					System.out.println("Average order costs: " + this.getAverageOfList(orderCostsEOQ.get(counter)));
					System.out.println("Average setup costs: " + this.getAverageOfList(setupCostsEOQ.get(counter)));
					System.out.println("Average production costs: " + this.getAverageOfList(productionCostsEOQ.get(counter)));
					System.out.println("Average purchase costs: " + this.getAverageOfList(purchaseCostsEOQ.get(counter)));
					System.out.println("Average holding costs: " + this.getAverageOfList(transportCostsEOQ.get(counter)));
					System.out.println("");
					
					// Average values of Exact model
					System.out.println("EXACT: ");
					System.out.println("Average total costs: " + this.getAverageOfList(totalScCostsExact.get(counter)) + " | MIN: " + Collections.min(totalScCostsExact.get(counter)) + " | MAX: " + Collections.max(totalScCostsExact.get(counter)));
					System.out.println("Average fill rate: " + this.getAverageOfList(scFillRateExact.get(counter)) + " | MIN: " + Collections.min(scFillRateExact.get(counter)) + " | MAX: " + Collections.max(scFillRateExact.get(counter)));
					System.out.println("Average fitness value: " + this.getAverageOfList(fitnessValueExact.get(counter)) + " | MIN: " + Collections.min(fitnessValueExact.get(counter)) + " | MAX: " + Collections.max(fitnessValueExact.get(counter)));
					System.out.println("Average holding costs: " + this.getAverageOfList(holdingCostsExact.get(counter)));
					System.out.println("Average penalty costs: " + this.getAverageOfList(penaltyCostsExact.get(counter)));
					System.out.println("Average order costs: " + this.getAverageOfList(orderCostsExact.get(counter)));
					System.out.println("Average setup costs: " + this.getAverageOfList(setupCostsExact.get(counter)));
					System.out.println("Average production costs: " + this.getAverageOfList(productionCostsExact.get(counter)));
					System.out.println("Average purchase costs: " + this.getAverageOfList(purchaseCostsExact.get(counter)));
					System.out.println("Average holding costs: " + this.getAverageOfList(transportCostsExact.get(counter)));
					System.out.println("");
					
					// Average values of Heuristic
					System.out.println("HEURISTIC: ");
					System.out.println("Average total costs: " + this.getAverageOfList(totalScCostsHeuristic.get(counter)) + " | MIN: " + Collections.min(totalScCostsHeuristic.get(counter)) + " | MAX: " + Collections.max(totalScCostsHeuristic.get(counter)));
					System.out.println("Average fill rate: " + this.getAverageOfList(scFillRateHeuristic.get(counter)) + " | MIN: " + Collections.min(scFillRateHeuristic.get(counter)) + " | MAX: " + Collections.max(scFillRateHeuristic.get(counter)));
					System.out.println("Average fitness value: " + this.getAverageOfList(fitnessValueHeuristic.get(counter)) + " | MIN: " + Collections.min(fitnessValueHeuristic.get(counter)) + " | MAX: " + Collections.max(fitnessValueHeuristic.get(counter)));
					System.out.println("Average holding costs: " + this.getAverageOfList(holdingCostsHeuristic.get(counter)));
					System.out.println("Average penalty costs: " + this.getAverageOfList(penaltyCostsHeuristic.get(counter)));
					System.out.println("Average order costs: " + this.getAverageOfList(orderCostsHeuristic.get(counter)));
					System.out.println("Average setup costs: " + this.getAverageOfList(setupCostsHeuristic.get(counter)));
					System.out.println("Average production costs: " + this.getAverageOfList(productionCostsHeuristic.get(counter)));
					System.out.println("Average purchase costs: " + this.getAverageOfList(purchaseCostsHeuristic.get(counter)));
					System.out.println("Average holding costs: " + this.getAverageOfList(transportCostsHeuristic.get(counter)));
					System.out.println("");
					
					// Increment counter by 1
					counter += 1;
				}
			}
		}	
	}
	
	/** Returns average of a list of double
	 * @param list
	 * @return average of list
	 */
	public double getAverageOfList(List<Double> list) {
		double sum = 0.0;
		for(int i=0; i<list.size(); i++) {
			sum += list.get(i);
		}
		return (double) sum / list.size();
	}

	/**
	 * @return the knowledgeBasesExact
	 */
	public List<Map<String, KnowledgeBase>> getKnowledgeBasesExact() {
		return knowledgeBasesExact;
	}

	/**
	 * @return the knowledgeBasesHeuristic
	 */
	public List<Map<String, KnowledgeBase>> getKnowledgeBasesHeuristic() {
		return knowledgeBasesHeuristic;
	}

	/**
	 * @return the trainingDemand
	 */
	public List<Map<String, TreeMap<Date, Integer>>> getTrainingDemand() {
		return trainingDemand;
	}

	/**
	 * @return the trainingLeadTime
	 */
	public List<Map<String, TreeMap<Date, Integer>>> getTrainingLeadTime() {
		return trainingLeadTime;
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
}