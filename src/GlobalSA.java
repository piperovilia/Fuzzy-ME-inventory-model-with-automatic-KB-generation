import java.io.File;
import java.util.*;

/**
 * SIMULATED ANNEALING (GLOBAL):
 * A class to optimize the rule bases of SC by applying a Simulated Annealing algorithm.
 * This algorithm is applied to the entire SC --> considers all the information along the chain.
 */
public class GlobalSA {
	private SC sc;
	private int numLabels;
	private KnowledgeBaseSet initKbSet; // Initial knowledge base set
	private List<String[]> initGlobalRuleBase; // Contains the rules of all initial RBs in the SC
	private List<String[]> currGlobalRuleBase; // Contains the rules of all current RBs in the SC
	private double initFitVal; // Initial fitness value
	private double bestFitVal; // Stores currently best fitness value
	private double lastFitVal; // Stores last fitness value
	private double alpha;
	private double temp;
	private int numSwitchesPerRun;
	private int stoppingCount;
	private double gammaFit;
	private double phiFit;
	private int iteration;
	private List<File> files;
	private int counterNoProgress; // Counts how many iterations have been applied without having progress
    private long runTime;
	private Random random;
	private TreeSet<KnowledgeBaseSet> knowledgeBaseSets; // Collection of knowledge base sets that were created in the course of the algorithm
	
	/**
	 * @param sc
	 * @param numLabels Number of labels per fuzzy variable
	 * @param inventory
	 * @param material
	 * @param initKb Initial knowledge base
	 * @param alpha Coefficient with which current temperature is multiplied
	 * @param temp Initial temperature
	 * @param numSwitchesPerRun Number of switched rule outputs per iteration
	 * @param stoppingCount Number of iterations without improvement to terminate algorithm
	 * @param gammaFit Parameter for prioritizing total SC costs
	 * @param phiFit Parameter for prioritizing fill rate
	 */
	public GlobalSA(SC sc, int numLabels, KnowledgeBaseSet initKbSet, double alpha, double temp, int numSwitchesPerRun, int stoppingCount, double gammaFit, double phiFit) {
		this.sc = sc;
		this.initKbSet = initKbSet;
		this.numLabels = numLabels;
		// Initialize global rule base set --> contains the rules of all KBs in the SC
		this.initGlobalRuleBase = new ArrayList<String[]>();
		for (var entry : initKbSet.getKnowledgeBases().entrySet()) {
			HashSet<String[]> ruleBase = entry.getValue().getRuleBase();
		    initGlobalRuleBase.addAll(ruleBase);
		}
		// Tuning parameters
		this.alpha = alpha;
		this.temp = temp;
		this.numSwitchesPerRun = numSwitchesPerRun;
		this.stoppingCount = stoppingCount;
		this.gammaFit = gammaFit;
		this.phiFit = phiFit;
		this.iteration = 1;
		this.files = new ArrayList<File>();
		this.counterNoProgress = 0;
		this.random = new Random();
		this.knowledgeBaseSets = new TreeSet<KnowledgeBaseSet>();
	}
	
	/** Computes the initial objective value given the initial data and rule bases for a particular inventory and material.
	 * @throws Exception
	 */
	public void initialize() throws Exception {
		knowledgeBaseSets.add(initKbSet);
		initFitVal = initKbSet.getFitnessValue();
		this.currGlobalRuleBase = new ArrayList<String[]>(initGlobalRuleBase);
		this.bestFitVal = initFitVal;
		this.lastFitVal = initFitVal;
	}

	/**
	 * Runs the Simulated Annealing algorithm until a stopping criterion is met.
	 * @throws Exception 
	 */
	public void runSimulatedAnnealing() throws Exception {
		this.initialize();
		while(counterNoProgress < stoppingCount) {
			// Create candidate solution
			List<String[]> candidateSolution = this.getRandomCandidateSolution(numSwitchesPerRun);
			// Accept or reject new candidate solution
			this.acceptOrReject(candidateSolution,true);
			// Update temperature
			temp = temp * alpha;
		}
		// Delete .fcl files in the directory
		this.deleteFclFiles();
	}
	
	/** Accepts or rejects a given candidate solution according to the total SC inventory costs incurred by the new solution.
	 * @param candidateSolution Candidate solution
	 * @param multiEchelon Indicates whether the method is applied in the context of a multi-echelon optimization (true) or not (false)
	 * @throws Exception 
	 */
	public void acceptOrReject(List<String[]> candidateSolution, boolean multiEchelon) throws Exception {
		// Retrieve objective value of new rule base set by applying the new Fuzzy Knowledge Bases
		KnowledgeBaseSet newKbSet = this.getNewKnowledgeBaseSet(candidateSolution);
		// Add new KnowledgeBase object to list
		knowledgeBaseSets.add(newKbSet);
		double newFitVal = newKbSet.getFitnessValue();
		double newObjVal = newKbSet.getObjVal();
		double newFillRate = newKbSet.getFillRate();
		if(newFitVal > lastFitVal) {
			// Directly accept new rule base
			currGlobalRuleBase = candidateSolution;
			lastFitVal = newFitVal;
		} else {
			// New objective value is not better than last one --> Decide if you accept it anyways
			double rand = random.nextDouble();
			double deltaX = Math.abs(newFitVal-lastFitVal);
			if(rand < Math.exp(-deltaX/temp)) {
				// Accept new rule base
				currGlobalRuleBase = candidateSolution;
				lastFitVal = newFitVal;
			}
		}
		// Check if new objective value is best so far
		if(newFitVal > bestFitVal) {
			bestFitVal = newFitVal;
			// Set the counter for not making progress to 0
			counterNoProgress = 0;
		} else {
			// No improvement was made --> increase the counter for not making progress by 1
			counterNoProgress += 1;
		}
	}
	
	/** Creates and returns a new knowledge base set using the candidate solution
	 * @param candidateRbSet New candidate solution
	 * @throws Exception
	 */
	public KnowledgeBaseSet getNewKnowledgeBaseSet(List<String[]> candidateRbSet) throws Exception {
		int pos = 0; // Position in candidate rule base set
		
		Map<String,KnowledgeBase> newKnowledgeBases = new HashMap<String,KnowledgeBase>();
		
		for (var entry : initKbSet.getKnowledgeBases().entrySet()) {
			KnowledgeBase oldKb = entry.getValue();
			// Extract information of old knowledge base to create new knowledge base object
			Inventory inventory = oldKb.getInventory();
			Material material = oldKb.getMaterial();
			double[] dataBase = oldKb.getDataBase();
			// Obtain new rule base
			HashSet<String[]> ruleBase = new HashSet<String[]>();
			int numRules = oldKb.getRuleBase().size();
			// Extract corresponding rule base
			for(int i=pos; i<pos+numRules; i++) {
				ruleBase.add(candidateRbSet.get(i));
			}
			// Create new knowledge base
			KnowledgeBase newKb = new KnowledgeBase(sc, inventory, material, dataBase, ruleBase, true, gammaFit, phiFit);
			// Add new knowledge base to knowledge base set
			newKnowledgeBases.put(entry.getKey(), newKb);
			// Go to next rule base in candidate rule base set
			pos += numRules;
		}
		
		// Create new knowledge base set object
		KnowledgeBaseSet newKbSet = new KnowledgeBaseSet(sc, newKnowledgeBases, gammaFit, phiFit);
		
		return newKbSet;
	}
	
	/** Selects a candidate for the neighborhood of solutions.
	 *  A neighborhood is defined as: randomly pick 'numChanges' rules and for each rule randomly change the output of the order quantity to a neighboring label.
	 * @param numChanges
	 */
	public List<String[]> getRandomCandidateSolution(int numChanges) {
		// Retrieve the current rule base
		List<String[]> candidateSolution = new ArrayList<String[]>(currGlobalRuleBase);
		// Retrieve the total number of rules
		int numRules = candidateSolution.size();
		List<Integer> listIndices = new ArrayList<Integer>();
		for(int i=0; i<numRules; i++) {
			listIndices.add(i);
		}
		Collections.shuffle(listIndices);
		// Conduct changes of the output
		for(int i=0; i<numChanges; i++) {
			// Pick a random rule index
			int randInt = listIndices.get(i);
			// Retrieve rule at random index
			String[] rule = candidateSolution.get(randInt);
			// Retrieve output
			String output = rule[rule.length - 1];
			// Take correct label names
			List<String> labelNames;
			if(numLabels == 3) {
				labelNames = sc.getLabels3();
			} else {
				labelNames = sc.getLabels5();
			}
			if(output.equals(labelNames.get(0))) {
				// If output has lowest label, take one label higher
				rule[rule.length - 1] = labelNames.get(1);
			} else if(output.equals(labelNames.get(labelNames.size() - 1))) {
				// If output has highest label, take one label lower
				rule[rule.length - 1] = labelNames.get(labelNames.size() - 2);
			} else {
				// If output is in the middle, take one label lower
				int index = labelNames.indexOf(output);
				rule[rule.length - 1] = labelNames.get(index - 1);
			}
		}
		return candidateSolution;
	}
	
	/**
	 * @return Knowledge base set with lowest total SC costs
	 * @throws Exception
	 */
	public KnowledgeBaseSet getBestKnowledgeBaseSet() throws Exception {
		return knowledgeBaseSets.first();
	}
	
	/**
	 * @return the sc
	 */
	public SC getSc() {
		return sc;
	}	
	
	/**
	 * @return the initKb
	 */
	public KnowledgeBaseSet getinitKbSet() {
		return initKbSet;
	}

	/**
	 * @return the objective value of the initial knowledge base set
	 */
	public double getInitObjValue() {
		return initKbSet.getObjVal();
	}
	
	/**
	 * @return the fill rate of the initial knowledge base set
	 */
	public double getInitFillRate() {
		return initKbSet.getObjVal();
	}
	
	/**
	 * @return the fitness value of the initial knowledge base set
	 */
	public double getInitFitVal() {
		return initKbSet.getFitnessValue();
	}

	/**
	 * @return the objective value of the best knowledge base set
	 */
	public double getBestObj() {
		return knowledgeBaseSets.first().getObjVal();
	}

	/**
	 * @return the fill rate of the best knowledge base set
	 */
	public double getFillRateAtBest() {
		return knowledgeBaseSets.first().getFillRate();
	}
	
	/**
	 * @return the fitness value of the best knowledge base set
	 */
	public double getBestFitVal() {
		return knowledgeBaseSets.first().getFitnessValue();
	}
	
	/**
	 * @return the run time
	 */
	public long getRunTime() {
		return runTime;
	}

	/** Prints rule bases of best knowledgeBaseSet.
	 * @param inventory
	 * @param material
	 */
	public void printBestRb() {
		for (var entry : knowledgeBaseSets.first().getKnowledgeBases().entrySet()) {
			KnowledgeBase kb = entry.getValue();
		    System.out.println("Inventory/Material: " + kb.getInventory().getInventoryId() + "/" + kb.getMaterial().getMaterialCode());
		    kb.printRb();
		    System.out.println("");
		}
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