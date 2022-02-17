import java.io.File;
import java.util.List;
import java.util.Map;

/**
 * A class to model a set of all knowledge bases in the SC (relevant for multi-echelon case).
 * A KnowledgeBaseSet automatically runs a complete simulation of the SC to directly compute the total SC costs, the SC fill rate and the fitness value linked that result from the given KBs.
 * In list and maps, Knowledge Base Sets are sorted according to their corresponding fitness values in a decreasing order.
 */
public class KnowledgeBaseSet implements Comparable<KnowledgeBaseSet> {
	private SC sc;
	private Map<String,KnowledgeBase> knowledgeBases;
	private double gammaFit;
	private double phiFit;
	private double objVal;
	private double fillRate;
	private double maxScCosts;
	private double fitnessValue;
	private Simulation simulation;
	
	/**
	 * @param sc Supply chain object
	 * @param knowledgeBases List of knowledge bases (a KB per inventory-material)
	 * @param gammaFit Parameter for prioritizing the SC costs
	 * @param phiFit Parameter for prioritizing the fill rate
	 * @throws Exception
	 */
	public KnowledgeBaseSet(SC sc, Map<String, KnowledgeBase> knowledgeBases, double gammaFit, double phiFit) throws Exception {
		this.sc = sc;
		this.knowledgeBases = knowledgeBases;
		this.gammaFit = gammaFit;
		this.phiFit = phiFit;
		
		// Write all .fcl files (one per inventory-material) 
		for (var entry : knowledgeBases.entrySet()) {
			Writer writer = new Writer(sc);
			if(entry.getValue().getInventory().getStageType().equals("PLANT")) {
				writer.writeFCL(entry.getValue().getInventory(), entry.getValue().getMaterial(), entry.getValue().getDataBase(), entry.getValue().getRuleBase(), sc.getVarNamesL(), sc.getLabels3(), 0);
			} else if(entry.getValue().getInventory().getStageType().equals("WH")) {
				writer.writeFCL(entry.getValue().getInventory(), entry.getValue().getMaterial(), entry.getValue().getDataBase(), entry.getValue().getRuleBase(), sc.getVarNamesM(), sc.getLabels3(), 0);
			} else {
				writer.writeFCL(entry.getValue().getInventory(), entry.getValue().getMaterial(), entry.getValue().getDataBase(), entry.getValue().getRuleBase(), sc.getVarNamesS(), sc.getLabels3(), 0);
			}
		}
		
		// Retrieve the maximal SC costs by summing up the maximal costs per inventory-material
		this.maxScCosts = 0.0;
		for(String tuple : sc.getTuples()) {
			maxScCosts += sc.getMaxScCost().get(tuple);
		}
		
		// Compute objective value, fill rate and fitness value of this knowledge base set
		this.simulation = new Simulation(sc,0);
		simulation.simulateSupplyChain(true);
		this.objVal = simulation.getTotalScCosts();
		this.fillRate = simulation.getSCFillRate();
		this.computeFitnessValue();
		// Delete all .fcl files
		this.deleteFclFiles();
	}
	
	/**
	 * @return the sc
	 */
	public SC getSc() {
		return sc;
	}
	
	/**
	 * Computes fitness value of this knowledge base set
	 */
	public void computeFitnessValue() {
		double normalizedCost = (1- (double) objVal/maxScCosts);   
		fitnessValue = Math.pow(normalizedCost, gammaFit) * Math.pow(fillRate, phiFit);
	}
	
	/**
	 * @return the fitnessValue
	 */
	public double getFitnessValue() {
		return fitnessValue;
	}
	
	/**
	 * @param fitnessValue the fitnessValue to set
	 */
	public void setFitnessValue(double fitnessValue) {
		this.fitnessValue = fitnessValue;
	}

	/**
	 * @return the knowledgeBases
	 */
	public Map<String, KnowledgeBase> getKnowledgeBases() {
		return knowledgeBases;
	}

	/**
	 * @return the objVal
	 */
	public double getObjVal() {
		return objVal;
	}
	
	/**
	 * @return the fillRate
	 */
	public double getFillRate() {
		return fillRate;
	}
	
	/**
	 * @return the simulation
	 */
	public Simulation getSimulation() {
		return simulation;
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
     * Enables sorting lists of knowledge bases according to their objective value in an ascending order. 
     */
    @Override public int compareTo(KnowledgeBaseSet kbSet) {
        if (this.getFitnessValue() < kbSet.getFitnessValue()) {
            // If current object is greater, then return 1
            return 1;
        }
        else if (this.getFitnessValue() > kbSet.getFitnessValue()) {
            // If current object is greater, then return -1
            return -1;
        }
        else {
            // If current object is equal to o, then return 0
            return 1;
        }
    }
}