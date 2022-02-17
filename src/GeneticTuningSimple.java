import java.util.*;

/**
 * LOCAL LA-TUNING algorithm (WITHOUT RB reduction):
 * This class tunes the membership function parameters and the number of rules in a RB by applying a Genetic Algorithm.
 * This algorithm is applied to a single inventory and material.
 */
public class GeneticTuningSimple {
	private SC sc;
	private KnowledgeBase initKb;
	private Inventory inventory;
	private Material material;
	private int numLabels;
	private int N;
	private Random rand;
	private int numGenes; // Total number of genes in chromosome
	private TreeMap<KnowledgeBase,double[]> population;
	private ArrayList<double[]> currentOffspring; // Contains the offspring of last iteration
	private String tuple;
	private int bitsGene;
	private Map<Integer,String> grayCodes; // Map that contains a gray code for each integer
	private Map<double[],Integer> grayCodeRanges; // Ranges used to transform a double to an integer
	private double initThreshold; // Initial threshold for a restart
	private double threshold; // Current threshold for a restart
	private double phi;
	private double gammaFit;
	private double phiFit;
	private int counterCycles;
	private double[] initIndividual;
	private List<Double> fitnessValues;
	private List<Double> objValues;
	private List<Double> fillRates;
	
	/** Default constructor.
	 * @param sc
	 * @param initKb Initial knowledge base
	 * @param numLabels Number of labels per system variable
	 * @param N Number of individuals in population
	 * @param bitsGene Bits per gene for transforming real-coded array to gray coded array
	 * @param phi Percentage for decreasing ltThreshold after every iteration -> Value between 0 and 1
	 * @param gammaFit Parameter for prioritizing SC costs
	 * @param phiFit Parameter for prioritizing fill rate
	 * @throws Exception 
	 */
	public GeneticTuningSimple(SC sc, KnowledgeBase initKb, int numLabels, int N, int bitsGene, double phi, double gammaFit, double phiFit) throws Exception {
		this.sc = sc;
		this.N = N;
		this.initKb = initKb;
		this.inventory = initKb.getInventory();
		this.material = initKb.getMaterial();
		this.numLabels = numLabels;
		this.population = new TreeMap<KnowledgeBase,double[]>();
		this.currentOffspring = new ArrayList<double[]>();
		// Number of genes in a chromosome/individual (Amplitude variation also applied to output variable)
		this.numGenes = 2 * inventory.getNumVars()*numLabels;
		double[] individual = new double[numGenes];
		// Initialize population with an individual having all MF genes with value 0 (no displacement or amplitude variation)
		for(int j=0; j<numGenes; j++) {
			individual[j] = 0;
		}
		initIndividual = individual;
		// Add first individual to the population
		population.put(initKb, individual);
		// Random generator
		this.rand = new Random();
		// "Inventory,Material" String combination 
		this.tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
		this.grayCodes = new HashMap<Integer,String>();
		this.bitsGene = bitsGene;
		this.grayCodeRanges = new HashMap<double[],Integer>();
		this.generateGrayCodes(bitsGene);
		this.initThreshold = (numGenes * bitsGene) / 4;
		this.threshold = initThreshold;
		this.phi = phi;
		this.gammaFit = gammaFit;
		this.phiFit = phiFit;
		this.counterCycles = 0;
		this.fitnessValues = new ArrayList<Double>();
		fitnessValues.add(initKb.getFitnessValue());
		this.objValues = new ArrayList<Double>();
		objValues.add(initKb.getObjVal());
		this.fillRates = new ArrayList<Double>();
		fillRates.add(initKb.getFillRate());
	}
	
	/** Adds N individuals to the population by taking random values [-0.5,0.5) for MF genes.
	 * @throws Exception 
	 */
	public void initialize() throws Exception {
		for(int i=0; i<N-1; i++) {
			double[] newIndividual = new double[numGenes];
			// Add random values to new individual
			for(int j=0; j<numGenes; j++) {
				// Generate real value between -0.5 and 0.5
				double randVal = rand.nextDouble() - 0.5;
				newIndividual[j] = randVal;		
			}
			// Create knowledge base
			KnowledgeBase newKb = this.getKbFromIndividual(newIndividual);
			// Add new individual to population
			population.put(newKb, newIndividual);
		}
	}
	
	/** Runs genetic algorithm 'numIterations' iterations or or 'numCycles' cycles.
	 * @param numIterations Number of iterations
	 * @throws Exception
	 */
	public void runGeneticTuning(int numCycles, int maxIterations) throws Exception {
		int counterIterations = 0;
		// Initialize population
		this.initialize();
		// Store currently best fitness value, objective value and fill rate
		fitnessValues.add(this.getBestKnowledgeBase().getFitnessValue());
		objValues.add(this.getBestKnowledgeBase().getObjVal());
		fillRates.add(this.getBestKnowledgeBase().getFillRate());
		while(counterCycles < numCycles && counterIterations < maxIterations) {
			// Create offspring
			this.createOffspring();
			// Update population: Pick N best individuals
			this.updatePopulation();
			fitnessValues.add(this.getBestKnowledgeBase().getFitnessValue());
			objValues.add(this.getBestKnowledgeBase().getObjVal());
			fillRates.add(this.getBestKnowledgeBase().getFillRate());
			// Empty lists of current offsprings
			currentOffspring.clear();
			// Print results
			/*
			System.out.println("Iteration " + (counterIterations+1) + " completed.");
			System.out.println("Current best fitness value: " + this.getBestFitVal());
			System.out.println("Current objective value: " + this.getBestObjVal());
			System.out.println("L_t threshold: " + threshold);
			System.out.println("");
			*/
			counterIterations += 1;
		}
	}
	
	/**
	 * Leaves the fittest N individuals after offspring was created. Also, updates threshold and restarts population if necessary.
	 * @throws Exception 
	 */
	public void updatePopulation() throws Exception {
		// Population is sorted in terms of total SC costs in ascending order
		int count = 0;
		TreeMap<KnowledgeBase,double[]> newPopulation = new TreeMap<KnowledgeBase,double[]>();
		for (Map.Entry<KnowledgeBase,double[]> entry : population.entrySet()) {
			if (count >= N) {
				break;
			}
			newPopulation.put(entry.getKey(), entry.getValue());
			count++;
		}
		this.population = newPopulation;
		
		// Check whether new offspring is included in the updated population
		int numNewOffspring = 0;
		for (double[] offspring : currentOffspring) {
			for (var entry : population.entrySet()) {
				double[] individual = entry.getValue();
				if(Arrays.equals(offspring, individual)) {
					numNewOffspring += 1;
					break;
				}
			}
		}
		// Decrement threshold by 1 if no offspring is included in the new population
		if(numNewOffspring == 0) {
			threshold -= 1;
		}
		// Update threshold:  
		threshold = threshold - threshold * phi;
		// Restart population if threshold < 0.
		if(threshold < 0) {
			this.restart();
		}
	}
	
	/**
	 * Restarts population.
	 * @throws Exception
	 */
	public void restart() throws Exception {
		// Retrieve best individual of current best population
		KnowledgeBase bestKb = population.firstKey();
		double[] bestIndividual = population.firstEntry().getValue();
		// Create new population
		TreeMap<KnowledgeBase,double[]> newPopulation = new TreeMap<KnowledgeBase,double[]>();
		// Add best individual of old population to new population
		newPopulation.put(bestKb,bestIndividual);
		// Create the remaining N-1 individuals around the best individual
		for(int i=0; i<N-1; i++) {
			// New individual is randomly picked around the currently best individual
			double[] newIndividual = new double[numGenes];
			for(int j=0; j<numGenes; j++) {
				double randShift = -0.125 + (0.125 + 0.125) * rand.nextDouble();
				// If new gene < -0.5, take -0.5. Else
				if(bestIndividual[j] + randShift < -0.5) {
					newIndividual[j] = -0.5;
				}
				// If new gene > 0.5, take 0.5
				if(bestIndividual[j] + randShift > 0.5) {
					newIndividual[j] = 0.5;
				}
				// Else, adapt new gene (-0.5 <= new gene <= 0.5)
				if(bestIndividual[j] + randShift >= -0.5 && bestIndividual[j] + randShift <= 0.5) {
					newIndividual[j] = bestIndividual[j] + randShift;;
				} 
			}
			// Add new individual to new population
			KnowledgeBase newKb = this.getKbFromIndividual(newIndividual);
			newPopulation.put(newKb, newIndividual);
		}
		// Update population
		population = newPopulation;
		// Update threshold
		threshold = initThreshold;
		// Increment cycle counter by 1
		counterCycles += 1;
	}
	
	/**
	 * Returns list of all feasible pairing partners of each individual.
	 * @return
	 */
	public Map<double[],List<double[]>> getPairsMap() {
		// RETRIEVE PAIRS OF PARENTS WHICH CAN BE PAIRED
		List<double[]> firstParents = new ArrayList<double[]>();
		List<double[]> secondParents = new ArrayList<double[]>();
		List<double[]> individualsList = new ArrayList<double[]>(population.values());
		
		for (int i=0; i<individualsList.size(); i++) {
			for (int j=i+1; j<individualsList.size(); j++) {
				// Check whether hammingDistance/2 > currThreshold
				int hammingDistCt = this.getHammingDistance(individualsList.get(i), individualsList.get(j));
				if((double) hammingDistCt / 2 > threshold) {
					firstParents.add(individualsList.get(i));
					secondParents.add(individualsList.get(j));
				}
			}
		}

		// Map all possible pairing partners to one particular individual
		Map<double[],List<double[]>> pairsMap = new HashMap<>();
		// Initialize lists
		for(int i=0; i<firstParents.size(); i++) {
			List<double[]> emptyList = new ArrayList<double[]>();
			pairsMap.put(firstParents.get(i), emptyList);
		}
		// Add pairing partners to lists
		for(int i=0; i<firstParents.size(); i++) {
			pairsMap.get(firstParents.get(i)).add(secondParents.get(i));
		}	
		return pairsMap;
	}
	
	/**
	 * Creates new offspring by using PCBLX crossover.
	 * @throws Exception 
	 */
	public void createOffspring() throws Exception {	
		// Map all possible pairing partners to one particular individual
		Map<double[],List<double[]>> pairsMap = this.getPairsMap();
		
		// CREATE OFFSPRING
		double a = -0.5;
		double b = 0.5;
		
		// PCBLX CROSSOVER
		for (var entry : pairsMap.entrySet()) {
			double[] parent1 = entry.getKey();	
			// Retrieve pairing partner
		    List<double[]> potentialPartners = entry.getValue();
		    double[] parent2;
		    if(potentialPartners.size() == 1) {
				parent2 = entry.getValue().get(0);
		    } else {
				// Pick random partner
		    	parent2 = entry.getValue().get(rand.nextInt(entry.getValue().size()));
		    }
		    List<double[]> newOffspring = this.createPCBLXCrossover(parent1, parent2, a, b);
		    // Create new knowledge base sets with new offsprings
		    KnowledgeBase newKb1 = this.getKbFromIndividual(newOffspring.get(0));
		    KnowledgeBase newKb2 = this.getKbFromIndividual(newOffspring.get(1));
		    // Add 2 new offsprings to population
		    population.put(newKb1, newOffspring.get(0));
		    population.put(newKb2, newOffspring.get(1));
		}
	}
	
	/**
	 * Creates new offspring by applying PCBLX crossover operator. 
	 * @param parent1 First parent to pair
	 * @param parent2 Second parent to pair
	 * @param a
	 * @param b
	 * @return List of new offspring
	 */
	public List<double[]> createPCBLXCrossover(double[] parent1, double[] parent2, double a, double b) {
		List<double[]> newOffsprings = new ArrayList<double[]>();
		double[] offspring1 = new double[numGenes];
		double[] offspring2 = new double[numGenes];
		
		for(int i=0; i<numGenes; i++) {
			// 1st offspring
			double x1 = parent1[i];
			double y1 = parent2[i];
			double I1 = Math.abs(x1-y1);
			double l1 = Math.max(a, x1-I1);
			double u1 = Math.min(b, x1+I1);
			// Randomly (uniformly) choose number from interval [l_i,u_i]
			double z1 = l1 + (u1 - l1) * rand.nextDouble();
			offspring1[i] = z1;
			// 2nd offspring
			double x2 = parent2[i];
			double y2 = parent1[i];
			double I2 = Math.abs(x2-y2);
			double l2 = Math.max(a, x2-I2);
			double u2 = Math.min(b, x2+I2);
			// Randomly (uniformly) choose number from interval [l_i,u_i]
			double z2 = l2 + (u2 - l2) * rand.nextDouble();
			offspring2[i] = z2;
		}
		// Add offspring to list
		newOffsprings.add(offspring1);
		newOffsprings.add(offspring2);
		
		// Add new offsprings to list of current offsprings
		currentOffspring.add(offspring1);
		currentOffspring.add(offspring2);
		return newOffsprings;
	}
	
	/** Transforms a given individual to a knowledge base object and returns it.
	 * @param individual
	 * @return new knowledge base object
	 * @throws Exception
	 */
	public KnowledgeBase getKbFromIndividual(double[] individual) throws Exception {
		double[] dataBase = initKb.getDataBase().clone();
		HashSet<String[]> ruleBase = initKb.getRuleBase();
		// MF Tuning
		for(int i=0; i<inventory.getNumVars(); i++) {
			// Related to C_L
			int posDb = i * numLabels * 3;
			int posPop = i * numLabels;
			// Related to C_A
			int posPop2 = i * numLabels + inventory.getNumVars() * numLabels;		
			double maxShift = dataBase[posDb+2] - dataBase[posDb];
			for(int j=0; j<numLabels; j++) {
				// Related to C_L
				double shift = (maxShift/2) * individual[posPop+j];
				if(j == 0 && dataBase[posDb + 3*j + 1] + shift < 0) {
					shift = shift - (dataBase[posDb + 3*j + 1] + shift);
				}
				dataBase[posDb + 3*j] = dataBase[posDb + 3*j] + shift;
				dataBase[posDb + 3*j + 1] = dataBase[posDb + 3*j + 1] + shift;
				dataBase[posDb + 3*j + 2] = dataBase[posDb + 3*j + 2] + shift;
				// Related to C_A
				double shift2 = maxShift * individual[posPop2+j];
				// Left point of MF
				dataBase[posDb + 3*j] = dataBase[posDb + 3*j] - shift2/2;
				// Right point of MF
				dataBase[posDb + 3*j + 2] = dataBase[posDb + 3*j + 2] + shift2/2;
			}
		}
		// Create knowledge base object
		KnowledgeBase newKb = new KnowledgeBase(sc, inventory, material, dataBase, ruleBase, false, gammaFit, phiFit);
		return newKb;
	}
	
	/**
	 * @return The hamming distance between two 'real-coded' parents (vectors) represented as gray code vectors. 
	 * A hamming distance is the number of digits that are different in two vectors of binary digits.
	 */
	public int getHammingDistance(double[] parent1, double[] parent2) {
		long startTime = System.nanoTime();
		int hammingDistance = 0;
		// Transform parent vectors to gray code vectors (every entry of the parent vector is expressed as a gray code consisting of 'bitsGene' binary values)
		// Only compute hamming distance if parent vectors are of same size, otherwise it is 0.
		if(parent1.length == parent2.length) {
			int[] grayCodeArr1 = new int[numGenes * bitsGene];
			int[] grayCodeArr2 = new int[numGenes * bitsGene];
			for(int i=0; i<numGenes; i++) {
				// Get grayCodeArr1
				double val1 = parent1[i];
				for (var entry : grayCodeRanges.entrySet()) {
					if(val1 > entry.getKey()[0] && val1 <= entry.getKey()[1]) {
						int val1Int = entry.getValue();
						String grayCode1 = grayCodes.get(val1Int);
						for(int j=0; j<grayCode1.length(); j++) {
							String charAsString = String.valueOf(grayCode1.charAt(j)); 
							int charAsInt = Integer.parseInt(charAsString);
							grayCodeArr1[i*grayCode1.length() + j] = charAsInt;
						}
						break;
					}
				}
				// Get grayCodeArr2
				double val2 = parent2[i];
				for (var entry : grayCodeRanges.entrySet()) {
					if(val2 > entry.getKey()[0] && val2 <= entry.getKey()[1]) {
						int val2Int = entry.getValue();
						String grayCode2 = grayCodes.get(val2Int);
						for(int j=0; j<grayCode2.length(); j++) {
							String charAsString = String.valueOf(grayCode2.charAt(j)); 
							int charAsInt = Integer.parseInt(charAsString);
							grayCodeArr2[i*grayCode2.length() + j] = charAsInt;
						}
						break;
					}
				}
			}
			// Compute Hamming distance
			for(int k=0; k<grayCodeArr1.length; k++) {
				if(grayCodeArr1[k] != grayCodeArr2[k]) {
					hammingDistance += 1;
				}
			}
		}
		return hammingDistance;
	}
	
	/** This function generates all n bit Gray codes and adds the generated codes to grayCodes map. 
	 * @param n Number of bits
	 */
	public void generateGrayCodes(int n) {
		// Base case
	    if (n <= 0) {
	    	return;
	    }     
	    // 'arr' will store all generated codes
	    ArrayList<String> arr = new ArrayList<String> ();
	    // Start with one-bit pattern
	    arr.add("0");
	    arr.add("1");
	 
	    // Every iteration of this loop generates 2*i codes from previously generated i codes.
	    int i, j;
	    for (i = 2; i < (1<<n); i = i<<1) {
	        // Enter the previously generated codes again in arr[] in reverse order. Nor arr[] has double number of codes.
	        for (j = i-1 ; j >= 0 ; j--) {
	        	arr.add(arr.get(j));
	        }   
	        // Append 0 to the first half
	        for (j = 0 ; j < i ; j++) {
	        	arr.set(j, "0" + arr.get(j));
	        }
	            
	        // Append 1 to the second half
	        for (j = i ; j < 2*i ; j++) {
	        	arr.set(j, "1" + arr.get(j));
	        }            
	    }
	    // Print contents of arr[]
	    for (i = 0 ; i < arr.size() ; i++ ) {
	    	grayCodes.put(i, arr.get(i));
	    }
	    // Create mapping ranges for transforming doubles to decimals
	    double rangeDist = (double) 1/arr.size();
	    double position = -0.5;
	    for (i = 0 ; i < arr.size() ; i++ ) {
	    	double[] range = { position, position + rangeDist };
	    	grayCodeRanges.put(range, i);
	    	position += rangeDist;
	    }
	}
	
	/**
	 * Prints objective values and individuals of the current population.
	 */
	public void printPopulation() {
		for (var entry : population.entrySet()) {
			System.out.println("");
		    System.out.println("Objective value: " + entry.getKey().getObjVal());
		    for(double d : entry.getValue()) {
		    	System.out.print(Math.round(d * 100.0) / 100.0 + " | ");
		    }
		    System.out.println("");
		}
	}
	
	/**
	 * Prints a given individual.
	 */
	public void printIndividual(double[] individual) {
		System.out.println("Individual: ");
		for (double d : individual) {
		    System.out.print(d + " | ");
		}
		System.out.println("");
	}

	/**
	 * @return the population
	 */
	public TreeMap<KnowledgeBase, double[]> getPopulation() {
		return population;
	}
	
	/**
	 * @return the knowledge base with the lowest total SC costs
	 */
	public KnowledgeBase getBestKnowledgeBase() {
		return population.firstKey();
	}

	/**
	 * @return the bestObjVal
	 */
	public double getBestObjVal() {
		return population.firstKey().getObjVal();
	}
	
	/**
	 * @return the bestObjVal
	 */
	public double getBestFitVal() {
		return population.firstKey().getFitnessValue();
	}
	
	/**
	 * @return the fitnessValues
	 */
	public List<Double> getFitnessValues() {
		return fitnessValues;
	}

	/**
	 * @return the objValues
	 */
	public List<Double> getObjValues() {
		return objValues;
	}

	/**
	 * @return the fillRates
	 */
	public List<Double> getFillRates() {
		return fillRates;
	}
}