import java.io.File;
import java.util.*;

/**
 * SIMULATED ANNEALING (LOCAL):
 * A class to optimize a given rule base by applying a Simulated Annealing algorithm.
 * This algorithm is applied to a single inventory and material.
 */
public class SimulatedAnnealing {
	private SC sc;
	private Inventory inventory;
	private Material material;
	private KnowledgeBase initKb; // Initial knowledge base
	private int numLabels;
	private double[] dataBase; // Initial data base
	private HashSet<String[]> initRuleBase; // Initial rule base
	private HashSet<String[]> currRuleBase; // Stores current rule base
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
    private long startTime;
    private long runTime;
	private Random random;
	private TreeSet<KnowledgeBase> knowledgeBases; // Set of knowledge bases that were created in the course of the algorithm
	
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
	public SimulatedAnnealing(SC sc, int numLabels, Inventory inventory, Material material, KnowledgeBase initKb, double alpha, 
			double temp, int numSwitchesPerRun, int stoppingCount, double gammaFit, double phiFit) {
		super();
		this.sc = sc;
		this.inventory = inventory;
		this.material = material;
		this.initKb = initKb;
		this.numLabels = numLabels;
		this.dataBase = initKb.getDataBase();
		// Initialize rule bases
		this.initRuleBase = initKb.getRuleBase();
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
		this.knowledgeBases = new TreeSet<KnowledgeBase>();
	}
	
	/** Computes the initial fitness value given the initial data and rule bases for a particular inventory and material.
	 * @throws Exception
	 */
	public void initialize() throws Exception {
		knowledgeBases.add(initKb);
		this.initFitVal = initKb.getFitnessValue();
		this.currRuleBase = new HashSet<String[]>(initKb.getRuleBase());
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
			HashSet<String[]> candidateSolution = this.getRandomCandidateSolution(numSwitchesPerRun);
			// Accept or reject new candidate solution
			this.acceptOrReject(candidateSolution, false);
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
	public void acceptOrReject(HashSet<String[]> candidateSolution, boolean multiEchelon) throws Exception {
		// Retrieve objective value of new rule base by applying the new Fuzzy Knowledge Base
		KnowledgeBase newKb = new KnowledgeBase(sc, inventory, material, dataBase, candidateSolution, multiEchelon, gammaFit, phiFit);
		knowledgeBases.add(newKb);
		double newObjVal = newKb.getObjVal();
		double newFillRate = newKb.getFillRate();
		double newFitVal = newKb.getFitnessValue();
		if(newFitVal > lastFitVal) {
			// Directly accept new rule base
			currRuleBase = candidateSolution;
			lastFitVal = newFitVal;
		} else {
			// New objective value is not better than last one --> Decide if you accept it anyways
			double rand = random.nextDouble();
			double deltaX = Math.abs(newFitVal-lastFitVal);
			if(rand < Math.exp(-deltaX/temp)) {
				// Accept new rule base
				currRuleBase = candidateSolution;
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
	
	/** Selects a candidate for the neighborhood of solutions.
	 *  A neighborhood is defined as: randomly pick 'numChanges' rules and for each rule randomly change the output of the order quantity to a neighboring label.
	 * @param numChanges Number of rule outputs to be changed
	 */
	public HashSet<String[]> getRandomCandidateSolution(int numChanges) {
		// Retrieve the current rule base
		ArrayList<String[]> candidateSolution = new ArrayList<String[]>(currRuleBase);
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
				// Choose random direction: New label can be one lower or one higher
				boolean randBool = random.nextBoolean();
				// Retrieve index of current output label
				int index = labelNames.indexOf(output);
				if(randBool) {
					// Take left label
					rule[rule.length - 1] = labelNames.get(index - 1);
				} else {
					// Take right label
					rule[rule.length - 1] = labelNames.get(index + 1);
				}
			}
		}
		HashSet<String[]> candidateSolutionFinal = new HashSet<String[]>(candidateSolution);
		return candidateSolutionFinal;
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
	
	/**
	 * @return the sc
	 */
	public SC getSc() {
		return sc;
	}	
	
	/**
	 * @return the initKb
	 */
	public KnowledgeBase getInitKb() {
		return initKb;
	}

	/**
	 * @return the dataBase
	 */
	public double[] getDataBase() {
		return dataBase;
	}

	/**
	 * @return the bestRuleBase
	 */
	public HashSet<String[]> getBestRuleBase() {
		return knowledgeBases.first().getRuleBase();
	}

	/**
	 * @return the objective value of the initial knowledge base
	 */
	public double getInitObjValue() {
		return initKb.getObjVal();
	}
	
	/**
	 * @return the fill rate of the initial knowledge base
	 */
	public double getInitFillRate() {
		return initKb.getFillRate();
	}
	
	/**
	 * @return the fitness value of the initial knowledge base
	 */
	public double getInitFitVal() {
		return initFitVal;
	}

	/**
	 * @return the objective value of the best knowledge base
	 */
	public double getBestObj() {
		return knowledgeBases.first().getObjVal();
	}

	/**
	 * @return the fill rate of the best knowledge base
	 */
	public double getFillRateAtBest() {
		return knowledgeBases.first().getFillRate();
	}
	
	/**
	 * @return the fill rate of the best knowledge base
	 */
	public double getBestFitVal() {
		return knowledgeBases.first().getFitnessValue();
	}
	
	/**
	 * @return the run time
	 */
	public long getRunTime() {
		return runTime;
	}

	/** Prints best rule base obtained.
	 * @param inventory
	 * @param material
	 */
	public void printBestRb() {
		for(String[] a : knowledgeBases.first().getRuleBase()) {
			if(inventory.getNumVars() == 5) {
    			if(sc.isDemand()) {
    				System.out.println("IF " + "demand IS " + a[0] + " AND inventory position IS " + a[1] + " AND lead time IS " + a[2] +  " AND market price IS " + a[3] + /*" AND supplierReliability IS " + rule[3] + */ " THEN orderQuantity IS " + a[4] + ";");
    			} else {
    				System.out.println("IF " + "demand change IS " + a[0] + " AND inventory position IS " + a[1] + " AND lead time IS " + a[2] +  " AND market price IS " + a[3] + /*" AND supplierReliability IS " + rule[3] + */ " THEN orderQuantity IS " + a[4] + ";");
    			}
			} else if(inventory.getNumVars() == 4){
    			if(sc.isDemand()) {
    				System.out.println("IF " + "demand IS " + a[0] + " AND inventory position IS " + a[1] + " AND lead time IS " + a[2] + /*" AND supplierReliability IS " + rule[3] + */ " THEN order quantity IS " + a[3] + ";");
    			} else {
    				System.out.println("IF " + "demand change IS " + a[0] + " AND inventory position IS " + a[1] + " AND lead time IS " + a[2] + /*" AND supplierReliability IS " + rule[3] + */ " THEN order quantity IS " + a[3] + ";");
    			}
			} else {
    			if(sc.isDemand()) {
    				System.out.println("IF " + "demand IS " + a[0] + " AND inventory position IS " + a[1] + /*" AND supplierReliability IS " + rule[3] + */ " THEN order quantity IS " + a[2] + ";");
    			} else {
    				System.out.println("IF " + "demand change IS " + a[0] + " AND inventory position IS " + a[1] + /*" AND supplierReliability IS " + rule[3] + */ " THEN order quantity IS " + a[2] + ";");
    			}
			}
		}
	}
	
	/**
	 * @return Knowledge base with lowest total SC costs
	 * @throws Exception
	 */
	public KnowledgeBase getBestKnowledgeBase() throws Exception {
		return knowledgeBases.first();
	}
}