import java.util.*;

/**
 * GLOBAL LA-TUNING algorithm (WITH RB reduction):
 * This class tunes the membership function parameters and the number of rules in a RB by applying a Genetic Algorithm.
 * This algorithm is applied to the entire SC --> considers all the information along the chain.
 */
public class GlobalGT {
	private SC sc;
	private KnowledgeBaseSet initKbSet;
	private int numLabels;
	private int N;
	private Random rand;
	private int numGenesCt; // Number of genes in C_T part
	private int numGenesCs; // Number of genes in C_S part
	private int numGenes; // Total number of genes in chromosome
	private TreeMap<KnowledgeBaseSet,double[]> population;
	private ArrayList<double[]> currentCtOffspring; // Contains the offspring related to the C_T part of last iteration
	private ArrayList<double[]> currentCsOffspring; // Contains the offspring related to the C_S part of last iteration
	private String tuple;
	private int bitsGene;
	private Map<Integer,String> grayCodes; // Map that contains a gray code for each integer
	private Map<double[],Integer> grayCodeRanges; // Ranges used to transform a double to an integer
	private double initLtThreshold; // Initial threshold for restarting C_T part
	private double ltThreshold; // Current threshold for restarting C_T part
	private double lsThreshold; // Current threshold for restarting C_T part
	private double phi;
	private double gammaFit;
	private double phiFit;
	private boolean firstRestart;
	private double[] bestCsIndividual;
	private int counterCycles;
	
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
	public GlobalGT(SC sc, KnowledgeBaseSet initKbSet, int numLabels, int N, int bitsGene, double phi, double gammaFit, double phiFit) throws Exception {
		this.sc = sc;
		this.N = N;
		this.initKbSet = initKbSet;
		this.numLabels = numLabels;
		this.population = new TreeMap<KnowledgeBaseSet,double[]>();
		this.currentCtOffspring = new ArrayList<double[]>();
		this.currentCsOffspring = new ArrayList<double[]>();
		// Number of genes in a chromosome/individual
		// NOTE: in this version NO amplitude variation is applied to output variable 
		// ADD UP NUMGENES OF ALL KBs (GLOBAL)
		this.numGenesCt = 0;
		for(Inventory i : sc.getInventories()) {
			if(i.isActiveInv()) {
				int numMaterials = i.getMaterials().size();
				numGenesCt += i.getNumVars()*numMaterials*numLabels + (i.getNumVars()-1)*numMaterials*numLabels;
			}
		}
		this.numGenesCs = 0;
		// Add up the number of rules of each KB in the SC
		for (var entry : initKbSet.getKnowledgeBases().entrySet()) {
		    numGenesCs += entry.getValue().getRuleBase().size();
		}
		this.numGenes = numGenesCt + numGenesCs;
		double[] individual = new double[numGenes];
		// Initialize population with an individual having all MF genes with value 0 (no displacement or amplitude variation), 
		// and all rule genes with value 1 (all rules included)
		for(int j=0; j<numGenes; j++) {
			if(j < numGenesCt) {
				individual[j] = 0;
			} else {
				individual[j] = 1;
			}
		}
		// Add first individual to the population
		population.put(initKbSet, individual);
		// Random generator
		this.rand = new Random();
		this.grayCodes = new HashMap<Integer,String>();
		this.bitsGene = bitsGene;
		this.grayCodeRanges = new HashMap<double[],Integer>();
		this.generateGrayCodes(bitsGene);
		this.initLtThreshold = (numGenesCt * bitsGene) / 4;
		this.ltThreshold = initLtThreshold;
		this.lsThreshold = numGenesCs / 4;
		this.phi = phi;
		this.gammaFit = gammaFit;
		this.phiFit = phiFit;
		this.firstRestart = false;
		this.bestCsIndividual = new double[numGenesCs];
		this.counterCycles = 0;
	}
	
	/** Adds N individuals to the population by taking random values [-0.5,0.5) for MF genes and random values {0,1} for rule genes.
	 * @throws Exception 
	 */
	public void initialize() throws Exception {
		for(int i=0; i<N-1; i++) {
			double[] newIndividual = new double[numGenes];
			// Add random values to new individual
			for(int j=0; j<numGenes; j++) {
				if(j < numGenesCt) {
					// Generate real value between -0.5 and 0.5
					double randVal = rand.nextDouble() - 0.5;
					newIndividual[j] = randVal;			
				} else {
					// Generate 0 or 1
					int randVal = rand.nextDouble() >= 0.5? 1 : 0;
					newIndividual[j] = randVal;
				}
			}
			// Create knowledge base
			KnowledgeBaseSet newKbSet = this.getKbSetFromIndividual(newIndividual);
			// Add new individual to population
			population.put(newKbSet, newIndividual);
		}
	}
	
	/**  Runs genetic algorithm 'numIterations' iterations or or 'numCycles' cycles.
	 * @param numIterations Number of iterations
	 * @throws Exception
	 */
	public void runGeneticTuning(int numCycles, int maxIterations) throws Exception {
		int counterIterations = 0;
		// Initialize populations
		this.initialize();
		while(counterCycles < numCycles && counterIterations < maxIterations) {
			// Create offspring
			if(firstRestart) {
				this.createOffspringAfterRestart();
			} else {
				this.createOffspringBeforeRestart();
			}
			// Update population: Pick N best individuals
			this.updatePopulation();
			// Empty lists of current offsprings
			currentCtOffspring.clear();
			currentCsOffspring.clear();
			// Print results
			/*
			System.out.println("Iteration " + (counterIterations+1) + " completed.");
			System.out.println("Current best objective value: " + this.getBestObjVal());
			System.out.println("L_t threshold: " + ltThreshold + " | L_s threshold: " + lsThreshold);
			System.out.println("");
			*/
			counterIterations += 1;
		}
	}
	
	/**
	 * Leaves the fittest N individuals after offspring was created. Also, updates ltThreshold and restarts population if necessary.
	 * @throws Exception 
	 */
	public void updatePopulation() throws Exception {
		// Population is sorted in terms of total SC costs in ascending order
		int count = 0;
		TreeMap<KnowledgeBaseSet,double[]> newPopulation = new TreeMap<KnowledgeBaseSet,double[]>();
		for (Map.Entry<KnowledgeBaseSet,double[]> entry : population.entrySet()) {
			if (count >= N) {
				break;
			}
			newPopulation.put(entry.getKey(), entry.getValue());
			count++;
		}
		this.population = newPopulation;
		
		// Check whether new offspring is included in the updated population
		ArrayList<double[]> ctParts = new ArrayList<double[]>();
		ArrayList<double[]> csParts = new ArrayList<double[]>();
		
		for (var entry : population.entrySet()) {
			double[] ct = Arrays.copyOfRange(entry.getValue(), 0, numGenesCt);
			ctParts.add(ct);
			double[] cs = Arrays.copyOfRange(entry.getValue(), numGenesCt, numGenes);
			csParts.add(cs);
		}
		
		int numNewCtOffspring = 0;
		for (double[] offspring : currentCtOffspring) {
			for(double[] ctPart : ctParts) {
				if(Arrays.equals(offspring, ctPart)) {
					numNewCtOffspring += 1;
					break;
				}
			}
		}
		// Decrement ltThreshold by 1 if no C_T offspring is included in the new population
		if(numNewCtOffspring == 0) {
			ltThreshold -= 1;
		}
		
		int numNewCsOffspring = 0;
		for (double[] offspring : currentCsOffspring) {
			for(double[] csPart : csParts) {
				if(Arrays.equals(offspring, csPart)) {
					numNewCsOffspring += 1;
					break;
				}
			}
		}
		// Decrement ltThreshold by 1 if no C_T offspring is included in the new population
		if(numNewCsOffspring == 0) {
			lsThreshold -= 1;
		}
		
		// Update ltThreshold:  
		ltThreshold = ltThreshold - initLtThreshold * phi;
		// Restart population if ltThreshold and lsThreshold are < 0 BEFORE the first restart. After, restart if ltThreshold < 0.
		if(!firstRestart) {
			if(ltThreshold < 0 && lsThreshold < 0) {
				this.restart();
				firstRestart = true;
			}
		} else {
			if(ltThreshold < 0) {
				this.restart();
			}
		}
	}
	
	/**
	 * Restarts population.
	 * @throws Exception
	 */
	public void restart() throws Exception {
		// Retrieve best individual of current best population
		//System.out.println(population.size());
		KnowledgeBaseSet bestKbSet = population.firstKey();
		double[] bestIndividual = population.firstEntry().getValue();
		//this.printIndividual(bestIndividual);
		// Fix the C_S part of the best individual if 1st restart
		if(!firstRestart) {
			for(int i=numGenesCt; i<numGenes; i++) {
				bestCsIndividual[i-numGenesCt] = bestIndividual[i];
			}
		}
		// Create new population
		TreeMap<KnowledgeBaseSet,double[]> newPopulation = new TreeMap<KnowledgeBaseSet,double[]>();
		// Add best individual of old population to new population
		newPopulation.put(bestKbSet,bestIndividual);
		// Create the remaining N-1 individuals around the best individual
		for(int i=0; i<N-1; i++) {
			// New C_T part is randomly picked around best individual's C_T part
			double[] newIndividual = new double[numGenes];
			for(int j=0; j<numGenes; j++) {
				if(j < numGenesCt) {
					// Related to C_T part
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
				} else {
					// Related to C_S part
					// C_T part of best individual is fixed for all new individuals
					newIndividual[j] = bestCsIndividual[j-numGenesCt];
				}
			}
			// Add new individual to new population
			KnowledgeBaseSet newKbSet = this.getKbSetFromIndividual(newIndividual);
			newPopulation.put(newKbSet, newIndividual);
		}
		// Update population
		population = newPopulation;
		// Update ltThreshold
		ltThreshold = initLtThreshold;
		// Increment cycle counter by 1
		counterCycles += 1;
		//System.out.println("Cycle completed.");
	}
	
	/**
	 * Returns list of all feasible pairing partners of each individual.
	 * @param isCt Indicates whether pairs are retrieved for C_T or C_S part
	 * @return
	 */
	public Map<double[],List<double[]>> getPairsMap(boolean isCt) {
		// RETRIEVE PAIRS OF PARENTS WHICH CAN BE PAIRED
		// Related to C_T:
		List<double[]> firstParents = new ArrayList<double[]>();
		List<double[]> secondParents = new ArrayList<double[]>();
		List<double[]> individualsList = new ArrayList<double[]>(population.values());
		for (int i=0; i<individualsList.size(); i++) {
			for (int j=i+1; j<individualsList.size(); j++) {
				// Check whether hammingDistance/2 > ltThreshold or lsThreshold
				if(isCt) {
					int hammingDistCt = this.getHammingDistanceCt(individualsList.get(i), individualsList.get(j));
					if(hammingDistCt / 2 > ltThreshold) {
						firstParents.add(individualsList.get(i));
						secondParents.add(individualsList.get(j));
					}
				} else {
					// Check whether hammingDistance/2 > lsThreshold
					int hammingDistCs = this.getHammingDistanceCs(individualsList.get(i), individualsList.get(j));
					if((double) hammingDistCs / 2 > lsThreshold) {
						firstParents.add(individualsList.get(i));
						secondParents.add(individualsList.get(j));
					}
				}
			}
		}

		// Map all possible pairing partners to one particular C_T individual
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
	 * Creates new offspring by using PCBLX crossover for the C_T part and HUX crossover for the C_S part (Before 1st restart).
	 * @throws Exception 
	 */
	public void createOffspringBeforeRestart() throws Exception {	
		// Map all possible pairing partners to one particular C_T individual
		Map<double[],List<double[]>> pairsMapCt = this.getPairsMap(true);
		// Map all possible pairing partners to one particular C_S individual
		Map<double[],List<double[]>> pairsMapCs = this.getPairsMap(false);
		
		// CREATE OFFSPRING
		double a = -0.5;
		double b = 0.5;
		
		List<double[]> usedParentsCs = new ArrayList<double[]>();
		
		// PCBLX CROSSOVER related to C_T
		for (var entry : pairsMapCt.entrySet()) {
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
		    List<double[]> newOffspringCt = this.createPCBLXCrossover(parent1, parent2, a, b);
		    List<double[]> newOffspringCs;
	
		    // CASE 1: C_T offspring and C_S offspring
		    if(pairsMapCs.containsKey(parent1) && pairsMapCs.get(parent1).contains(parent2)) {
		    	newOffspringCs = this.createHUXCrossover(parent1, parent2);
		    	usedParentsCs.add(parent1);
		    	usedParentsCs.add(parent2);
		    } else if(pairsMapCs.containsKey(parent2) && pairsMapCs.get(parent2).contains(parent1)) {
		    	newOffspringCs = this.createHUXCrossover(parent2, parent1);
		    	usedParentsCs.add(parent1);
		    	usedParentsCs.add(parent2);
		    } else {
		    	// CASE 2: C_T offspring, no C_S offspring -> take old C_S part
				double[] parent1CsPart = new double[numGenesCs];
				double[] parent2CsPart = new double[numGenesCs];
				for(int i=numGenesCt; i<numGenes; i++) {
					parent1CsPart[i-numGenesCt] = parent1[i];
					parent2CsPart[i-numGenesCt] = parent2[i];
				}
				newOffspringCs = new ArrayList<double[]>();
		    	newOffspringCs.add(parent1CsPart);
		    	newOffspringCs.add(parent2CsPart);
		    }
	    	// Create 1st and 2nd offspring and add to population
		    this.concatenateAndAddToPopulation(newOffspringCt.get(0), newOffspringCs.get(0));
		    this.concatenateAndAddToPopulation(newOffspringCt.get(1), newOffspringCs.get(1));	
		}
		
		// CASE 3: No C_T offspring, C_S offspring -> take old C_T part
		for (var entry : pairsMapCs.entrySet()) {
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
			if(!usedParentsCs.contains(parent1) && !usedParentsCs.contains(parent2)) {
				List<double[]> newOffspringCs = this.createHUXCrossover(parent1, parent2);
				// Take old C_T part
				List<double[]> newOffspringCt = new ArrayList<double[]>();
				double[] parent1CtPart = new double[numGenesCt];
				double[] parent2CtPart = new double[numGenesCt];
				for(int i=0; i<numGenesCt; i++) {
					parent1CtPart[i] = parent1[i];
					parent2CtPart[i] = parent2[i];
				}
		    	newOffspringCt.add(parent1CtPart);
		    	newOffspringCt.add(parent2CtPart);
		    	
		    	// Create 1st and 2nd offspring and add to population
			    this.concatenateAndAddToPopulation(newOffspringCt.get(0), newOffspringCs.get(0));
			    this.concatenateAndAddToPopulation(newOffspringCt.get(1), newOffspringCs.get(1));
			}	 
		}
	}
	
	/**
	 * Creates new offspring by using PCBLX crossover for the C_T part and takes the C_S part of the best individual at the first restart (after 1st restart).
	 * @throws Exception 
	 */
	public void createOffspringAfterRestart() throws Exception {	
		int counterOffspringCt = 0;
		// Map all possible pairing partners to one particular C_T individual
		Map<double[],List<double[]>> pairsMapCt = this.getPairsMap(true);
		
		// CREATE OFFSPRING
		double a = -0.5;
		double b = 0.5;
		// Pick random pairs of parents and create offspring
		// PCBLX CROSSOVER related to C_T
		for (var entry : pairsMapCt.entrySet()) {
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
		    
		    List<double[]> newOffspringCt = this.createPCBLXCrossover(parent1, parent2, a, b);
		    // Take old C_S part
		    List<double[]> newOffspringCs = new ArrayList<double[]>();
			double[] parent1CsPart = new double[numGenesCs];
			double[] parent2CsPart = new double[numGenesCs];
			for(int i=numGenesCt; i<numGenes; i++) {
				parent1CsPart[i-numGenesCt] = parent1[i];
				parent2CsPart[i-numGenesCt] = parent2[i];
			}	
			
	    	newOffspringCs.add(parent1CsPart);
	    	newOffspringCs.add(parent2CsPart);

	    	// Create 1st and 2nd offspring and add to population
		    this.concatenateAndAddToPopulation(newOffspringCt.get(0), newOffspringCs.get(0));
		    this.concatenateAndAddToPopulation(newOffspringCt.get(1), newOffspringCs.get(1));
		    counterOffspringCt += 2;
		}
		// If no offspring can be created, decrement ltThreshold by 1
		if(counterOffspringCt == 0) {
			ltThreshold -= 1;
		}
	}
	
	/**
	 * Combines given C_T and C_S part to a new individual and adds individual to population.
	 * @param array1 C_T part
	 * @param array2 C_S part
	 * @throws Exception
	 */
	public void concatenateAndAddToPopulation(double[] array1, double[] array2) throws Exception {
	    double[] newOffspring = Arrays.copyOf(array1, array1.length + array2.length);
	    System.arraycopy(array2, 0, newOffspring, array1.length, array2.length);
	    KnowledgeBaseSet newKbSet = this.getKbSetFromIndividual(newOffspring);
	    population.put(newKbSet, newOffspring);
	}
	
	/**
	 * Creates new offspring of C_T parts by applying PCBLX crossover operator.
	 * @param parent1 First parent to pair
	 * @param parent2 Second parent to pair
	 * @param a
	 * @param b
	 * @return List of new offspring
	 */
	public List<double[]> createPCBLXCrossover(double[] parent1, double[] parent2, double a, double b) {
		List<double[]> newOffsprings = new ArrayList<double[]>();
		double[] offspringCt1 = new double[numGenesCt];
		double[] offspringCt2 = new double[numGenesCt];
		
		for(int i=0; i<numGenesCt; i++) {
			// 1st offspring
			double x1 = parent1[i];
			double y1 = parent2[i];
			double I1 = Math.abs(x1-y1);
			double l1 = Math.max(a, x1-I1);
			double u1 = Math.min(b, x1+I1);
			// Randomly (uniformly) choose number from interval [l_i,u_i]
			double z1 = l1 + (u1 - l1) * rand.nextDouble();
			offspringCt1[i] = z1;
			// 2nd offspring
			double x2 = parent2[i];
			double y2 = parent1[i];
			double I2 = Math.abs(x2-y2);
			double l2 = Math.max(a, x2-I2);
			double u2 = Math.min(b, x2+I2);
			// Randomly (uniformly) choose number from interval [l_i,u_i]
			double z2 = l2 + (u2 - l2) * rand.nextDouble();
			offspringCt2[i] = z2;
		}
		// Add offspring to list
		newOffsprings.add(offspringCt1);
		newOffsprings.add(offspringCt2);
		
		// Add new offsprings to list of current CT offsprings
		currentCtOffspring.add(offspringCt1);
		currentCtOffspring.add(offspringCt2);
		
		return newOffsprings;
	}
	
	/**
	 * Creates new offspring of C_S parts by applying HUX crossover operator.
	 * @param parent1 First parent to pair
	 * @param parent2 Second parent to pair
	 * @return List of new offspring
	 */
	public List<double[]> createHUXCrossover(double[] parent1, double[] parent2) {
		List<double[]> newOffsprings = new ArrayList<double[]>();
		double[] offspringCs1 = new double[numGenesCs];
		double[] offspringCs2 = new double[numGenesCs];
		
		for(int i=numGenesCt; i<numGenes; i++) {
			offspringCs1[i-numGenesCt] = parent1[i];
			offspringCs2[i-numGenesCt] = parent2[i];
		}
		
		List<Integer> indicesUnmatched = new ArrayList<Integer>();
		
		for(int i=numGenesCt; i<numGenes; i++) {
			if(parent1[i] != parent2[i]) {
				indicesUnmatched.add(i);
			}
		}
		
		// Shuffle list with indices
		Collections.shuffle(indicesUnmatched);
		for(int j=0; j<(int) indicesUnmatched.size()/2; j++) {
			offspringCs1[indicesUnmatched.get(j) - numGenesCt] = parent2[indicesUnmatched.get(j)];
			offspringCs2[indicesUnmatched.get(j) - numGenesCt] = parent1[indicesUnmatched.get(j)];
		}
		
		// Add offspring to list
		newOffsprings.add(offspringCs1);
		newOffsprings.add(offspringCs2);
		
		// Add new offsprings to list of current CS offsprings
		currentCsOffspring.add(offspringCs1);
		currentCsOffspring.add(offspringCs2);
		
		return newOffsprings;
	}
	
	/** Transforms a given individual to a knowledge base set object and returns it.
	 * @param individual
	 * @return new knowledge base set
	 * @throws Exception
	 */
	public KnowledgeBaseSet getKbSetFromIndividual(double[] individual) throws Exception {
		Map<String,KnowledgeBase> newKnowledgeBases = new HashMap<String,KnowledgeBase>();
		int globalPos = 0;
		int globalPosRb = numGenesCt;
		int counterRules = 0;
		// Go through all KBs (global)
		for (var entry : initKbSet.getKnowledgeBases().entrySet()) {
			KnowledgeBase kb = entry.getValue();
			Inventory inventory = kb.getInventory();
			Material material = kb.getMaterial();
			double[] dataBase = kb.getDataBase().clone();
			
			// MF Tuning
			for(int i=0; i<inventory.getNumVars(); i++) {
				// Related to C_L
				int posDb = i * numLabels * 3;
				int posPop = globalPos + i * numLabels;
				// Related to C_A
				int posPop2 = globalPos + i * numLabels + inventory.getNumVars() * numLabels;		
				double maxShift = dataBase[posDb+2] - dataBase[posDb];
				for(int j=0; j<numLabels; j++) {
					// Related to C_L
					double shift = (maxShift/2) * individual[posPop+j];
					dataBase[posDb + 3*j] = dataBase[posDb + 3*j] + shift;
					dataBase[posDb + 3*j + 1] = dataBase[posDb + 3*j + 1] + shift;
					dataBase[posDb + 3*j + 2] = dataBase[posDb + 3*j + 2] + shift;
					if(i<inventory.getNumVars()-1) {
						// Related to C_A
						double shift2 = maxShift * individual[posPop2+j];
						// Left point of MF
						dataBase[posDb + 3*j] = dataBase[posDb + 3*j] - shift2/2;
						// Right point of MF
						dataBase[posDb + 3*j + 2] = dataBase[posDb + 3*j + 2] + shift2/2;
					}
				}
			}
			globalPos += inventory.getNumVars()*numLabels + (inventory.getNumVars()-1)*numLabels;
			// Rule tuning
			HashSet<String[]> ruleBaseSet = new HashSet<String[]>(kb.getRuleBase());
			List<String[]> ruleBaseList = new ArrayList<String[]>(ruleBaseSet);
			List<String[]> ruleBase = new ArrayList<String[]>();
			int startPos = globalPosRb + counterRules;
			for(int i=startPos; i<startPos+ruleBaseSet.size(); i++) {
				if(individual[i] == 1) {
					ruleBase.add(ruleBaseList.get(i-startPos));
				}
			}
			counterRules += ruleBaseSet.size();
			HashSet<String[]> newRuleBase = new HashSet<String[]>(ruleBase);
			// Create knowledge base object
			KnowledgeBase newKb = new KnowledgeBase(sc, inventory, material, dataBase, newRuleBase, true, gammaFit, phiFit);
			newKnowledgeBases.put(entry.getKey(), newKb);
		}
		// Create knowledge base set
		KnowledgeBaseSet newKbSet = new KnowledgeBaseSet(sc, newKnowledgeBases, gammaFit, phiFit);
		return newKbSet;
	}
	
	/**
	 * @return The hamming distance between two 'real-coded' parents (vectors) represented as gray code vectors. 
	 * A hamming distance is the number of digits that are different in two vectors of binary digits.
	 */
	public int getHammingDistanceCt(double[] parent1, double[] parent2) {
		int hammingDistance = 0;
		// Transform parent vectors to gray code vectors (every entry of the parent vector is expressed as a gray code consisting of 'bitsGene' binary values)
		// Only compute hamming distance if parent vectors are of same size, otherwise it is 0.
		if(parent1.length == parent2.length) {
			int[] grayCodeArr1 = new int[numGenesCt * bitsGene];
			int[] grayCodeArr2 = new int[numGenesCt * bitsGene];
			for(int i=0; i<numGenesCt; i++) {
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
			// Compute hamming distance
			for(int k=0; k<grayCodeArr1.length; k++) {
				if(grayCodeArr1[k] != grayCodeArr2[k]) {
					hammingDistance += 1;
				}
			}
		}	
		return hammingDistance;
	}
	
	/**
	 * @return The hamming distance between two 'binary-coded' parents (vectors). 
	 * A hamming distance is the number of digits that are different in two vectors (arrays).
	 */
	public int getHammingDistanceCs(double[] parent1, double[] parent2) {
		int hammingDistance = 0;
		if(parent1.length == parent2.length) {
			// Make sure to go to C_S part of individual
			for(int i=numGenesCt; i<parent1.length; i++) {
				if(parent1[i] != parent2[i]) {
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
	 * Prints individual.
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
	public TreeMap<KnowledgeBaseSet, double[]> getPopulation() {
		return population;
	}
	
	/**
	 * @return the knowledge base with the lowest total SC costs
	 */
	public KnowledgeBaseSet getBestKnowledgeBaseSet() {
		return population.firstKey();
	}

	/**
	 * @return the bestObjVal
	 */
	public double getBestObjVal() {
		return population.firstKey().getObjVal();
	}
}
