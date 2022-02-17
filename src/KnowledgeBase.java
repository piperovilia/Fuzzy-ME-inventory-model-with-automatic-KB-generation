import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

/**
 * A class to model a Fuzzy Knowledge Base consisting of a Data Base and a Rule Base. A Knowledge Base is linked to a particular inventory and material.
 * If the Knowledge Base is locally applied, at creation it automatically runs a simulation of this inventory-material to compute the total SC costs, the SC fill rate and the fitness value.
 * In list and maps, Knowledge Bases are sorted according to their corresponding fitness values in a decreasing order.
 */
public class KnowledgeBase implements Comparable<KnowledgeBase> {
	private SC sc;
	private int numLabels;
	private Inventory inventory;
	private Material material;
	private Simulation simulation;
	private double objVal;
	private double fillRate;
	private double[] dataBase;
	private HashSet<String[]> ruleBase;
	private boolean multiEchelon;
	private double gammaFit;
	private double phiFit;
	private double maxScCosts;
	private double fitnessValue;
	
	/**
	 * @param sc Supply chain object
	 * @param inventory Inventory related to KB
	 * @param material Material related to KB
	 * @param dataBase DB related to KB
	 * @param ruleBase RB related to KB
	 * @param multiEchelon Indicates whether this KB is created in the context of a multi-echelon iventory system
	 * @param gammaFit Parameter for prioritizing SC costs
	 * @param phiFit Parameter for prioritizing fill rate
	 * @throws Exception
	 */
	public KnowledgeBase(SC sc, Inventory inventory, Material material, double[] dataBase,
			HashSet<String[]> ruleBase, boolean multiEchelon, double gammaFit, double phiFit) throws Exception {
		this.sc = sc;
		this.inventory = inventory;
		this.material = material;
		this.dataBase = dataBase;
		this.ruleBase = ruleBase;
		this.numLabels = dataBase.length / (inventory.getNumVars() * 3);
		this.multiEchelon = multiEchelon;
		this.gammaFit = gammaFit;
		this.phiFit = phiFit;
		
		List<String> labels = sc.getLabels3();
		if(numLabels == 5) {
			labels = sc.getLabels5();
		}
		
		// Computation of objective value, fill rate and fitenss value (only if on single echelon applied)
		if(!multiEchelon) {
			Writer writer = new Writer(sc);
			if(inventory.getStageType().equals("PLANT")) {
				writer.writeFCL(inventory, material, dataBase, ruleBase, sc.getVarNamesL(), labels, 0);
			} else if(inventory.getStageType().equals("WH")) {
				writer.writeFCL(inventory, material, dataBase, ruleBase, sc.getVarNamesM(), labels, 0);
			} else {
				writer.writeFCL(inventory, material, dataBase, ruleBase, sc.getVarNamesS(), labels, 0);
			}
			String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
			this.simulation = new Simulation(sc,0);
			simulation.simulateInventory(inventory, material, dataBase, ruleBase, true, multiEchelon);
			this.objVal = simulation.getTotalCosts().get(tuple);
			this.fillRate = simulation.getFillRate(tuple);
			maxScCosts = sc.getMaxScCost().get(tuple);
			this.computeFitnessValue();
			this.deleteFclFiles();
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

	/**
	 * @return the sc
	 */
	public SC getSc() {
		return sc;
	}
	
	/**
	 * @return the inventory
	 */
	public Inventory getInventory() {
		return inventory;
	}

	/**
	 * @return the material
	 */
	public Material getMaterial() {
		return material;
	}

	/**
	 * @return the simulation
	 */
	public Simulation getSimulation() {
		return simulation;
	}

	/**
	 * Computes fitness value of this knowledge base
	 */
	public void computeFitnessValue() {
		double normalizedCost = (1- (double) objVal/maxScCosts);   
		fitnessValue = Math.pow(normalizedCost, gammaFit) * Math.pow(fillRate, phiFit);
	}

	/**
	 * @return the fitness value of this knowledge base
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
	 * @return the dataBase
	 */
	public double[] getDataBase() {
		return dataBase;
	}

	/**
	 * @return the ruleBase
	 */
	public HashSet<String[]> getRuleBase() {
		return ruleBase;
	}

	/**
	 * @return the multiEchelon
	 */
	public boolean isMultiEchelon() {
		return multiEchelon;
	}
	
	/**
	 * @return the maxScCosts
	 */
	public double getMaxScCosts() {
		return maxScCosts;
	}

	/** Prints rule base.
	 * @param inventory
	 * @param material
	 */
	public void printRb() {
		for(String[] a : ruleBase) {
			if(inventory.getStageType().equals("PLANT")) {
				System.out.println("IF " + "demand IS " + a[0] + " AND inventory position IS " + a[1] + " AND lead time IS " + a[2] + " AND market price IS " + a[3] + /*" AND supplier reliability IS " + a[3] +*/ " THEN order quantity IS " + a[4]);
			} else if(inventory.getStageType().equals("WH")) {
				System.out.println("IF " + "demand IS " + a[0] + " AND inventory position IS " + a[1] + " AND lead time IS " + a[2] + /*" AND supplier reliability IS " + a[3] +*/ " THEN order quantity IS " + a[3]);
			} else {
				System.out.println("IF " + "demand IS " + a[0] + " AND inventory position IS " + a[1] + /*" AND supplier reliability IS " + a[3] +*/ " THEN order quantity IS " + a[2]);
			}			
		}
	}
	
    /**
     * Enables sorting lists of knowledge bases according to their fitness value in an ascending order. 
     */
    @Override public int compareTo(KnowledgeBase kb) {
        if (this.getFitnessValue() < kb.getFitnessValue()) {
            // If current object is greater, then return 1
            return 1;
        }
        else if (this.getFitnessValue() > kb.getFitnessValue()) {
            // If current object is greater, then return -1
            return -1;
        }
        else {
            // If current object is equal to o, then return 0
            return 1;
        }
    }
    
    /*
	@Override
	public boolean equals(Object o) {
		// null check
		if (o == null) {
			return false;
		}
		// this instance check
		if (this == o) {
			return true;
		} 
		// instanceof check and actual value check
		if (o instanceof KnowledgeBase) {
			KnowledgeBase otherKb = (KnowledgeBase) o;
			boolean isDbEqual = Arrays.equals(this.dataBase, otherKb.getDataBase());
			boolean isRbEqual = true;
			if(this.ruleBase.size() == otherKb.getRuleBase().size()) {
				for(String[] rule : this.ruleBase) {
					if(!otherKb.getRuleBase().stream().anyMatch(c -> Arrays.equals(c, rule))) {
						isRbEqual = false;
						break;
					}
				}
			} else {
				isRbEqual = false;
			}
			if(isDbEqual == true && isRbEqual == true && this.inventory.getInventoryId().equals(otherKb.getInventory().getInventoryId()) 
					&& this.material.getMaterialCode().equals(otherKb.getMaterial().getMaterialCode())) {
				return true;
			}
		}
		return false;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + Arrays.hashCode(dataBase);
		result = prime * result + ((inventory == null) ? 0 : inventory.hashCode());
		result = prime * result + ((material == null) ? 0 : material.hashCode());
		result = prime * result + ((ruleBase == null) ? 0 : ruleBase.hashCode());
		return result;
	}
	*/
}