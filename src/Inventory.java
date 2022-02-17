import java.util.*;

/**
 * A class to model an inventory.
 *
 */
public class Inventory {
	private String inventoryId;
	private String stageType;
	private String description;
	private HashMap<Material,Double> holdingCosts;
	private HashMap<Material,Double> shortageCosts;
	private HashMap<String,HashMap<Material,Double>> transportCosts;
	private HashMap<Material,Double> setupCosts;
	private HashMap<Material,Double> orderCosts;
	private HashMap<Material,Double> prodUnitCosts;
	// Each material has one demands/demandChanges map per successor
	private String predecessor;
	private List<String> successors;
	private List<Material> materials;
	private Date splitDate;
	private int numVars;
	private boolean isActiveInv;
	private TreeSet<Date> dates;
	
	/**
	 * @param inventoryId The ID
	 * @param stageType "P" (plant), "WH" (warehouse) or "DIST" (distribution center)
	 * @param description
	 * @param holdingCosts Map with holding costs per material
	 * @param shortageCosts Map with penalty costs per material
	 * @param transportCosts Map with transport costs
	 * @param setupCosts Map with setup costs per material
	 * @param orderCosts Map with order costs per material
	 * @param prodUnitCosts Map with production costs per unit per material
	 * @param predecessor Predecessor inventory of this inventory
	 * @param successors List of successors of this inventory
	 * @param materials List of material at this inventory
	 * @param splitDate Last date of training phase
	 * @param isActiveInv Indicates whether this inventory is a P, WH or DIST (active)
	 * @param dates List of all dates
	 */
	public Inventory(String inventoryId, String stageType, String description, HashMap<Material,Double> holdingCosts, HashMap<Material,Double> shortageCosts, 
			HashMap<String,HashMap<Material,Double>> transportCosts, HashMap<Material, Double> setupCosts, HashMap<Material, Double> orderCosts,
			HashMap<Material, Double> prodUnitCosts, String predecessor, List<String> successors, List<Material> materials, Date splitDate, boolean isActiveInv, TreeSet<Date> dates) {
		this.inventoryId = inventoryId;
		this.stageType = stageType;
		this.description = description;
		this.holdingCosts = holdingCosts;
		this.shortageCosts = shortageCosts;
		this.transportCosts = transportCosts;
		this.setupCosts = setupCosts;
		this.orderCosts = orderCosts;
		this.prodUnitCosts = prodUnitCosts;
		this.successors = successors;
		this.predecessor = predecessor;
		this.materials = materials;
		this.splitDate = splitDate;
		// Each type of inventory has a different number of fuzzy variables
		if(this.stageType.equals("PLANT")) {
			this.numVars = 5;
		} else if (stageType.equals("WH")) {
			this.numVars = 4;
		} else {
			this.numVars = 3;
		}
		this.isActiveInv = isActiveInv;
		this.dates = dates;
	}

	/**
	 * @return the inventoryId
	 */
	public String getInventoryId() {
		return inventoryId;
	}

	/**
	 * @return the stage type of this inventory
	 */
	public String getStageType() {
		return stageType;
	}
	
	/**
	 * @return the description of this inventory
	 */
	public String getDescription() {
		return description;
	}

	/**
	 * @return the holdingCosts
	 */
	public HashMap<Material, Double> getHoldingCosts() {
		return holdingCosts;
	}

	/**
	 * @return the shortageCosts
	 */
	public HashMap<Material, Double> getShortageCosts() {
		return shortageCosts;
	}

	/**
	 * @return the transport costs
	 */
	public HashMap<String,HashMap<Material,Double>> getTransportCosts() {
		return transportCosts;
	}
	
	/**
	 * @return the setup costs
	 */
	public Map<Material, Double> getSetupCosts() {
		return setupCosts;
	}
	
	/**
	 * @return the orderCosts
	 */
	public HashMap<Material, Double> getOrderCosts() {
		return orderCosts;
	}

	/**
	 * @return the production unit costs
	 */
	public HashMap<Material, Double> getProdUnitCosts() {
		return prodUnitCosts;
	}
	
	/**
	 * @return the predecessor
	 */
	public String getPredecessor() {
		return predecessor;
	}

	/**
	 * @return the successors
	 */
	public List<String> getSuccessors() {
		return successors;
	}

	/** Returns materials being stored at this inventory.
	 * @return the materials
	 */
	public List<Material> getMaterials() {
		return materials;
	}

	/**
	 * @return the splitDate
	 */
	public Date getSplitDate() {
		return splitDate;
	}

	/**
	 * @return the numVars
	 */
	public int getNumVars() {
		return numVars;
	}

	/**
	 * @return the isActiveInv
	 */
	public boolean isActiveInv() {
		return isActiveInv;
	}

	/**
	 * @return the dates
	 */
	public TreeSet<Date> getDates() {
		return dates;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((inventoryId == null) ? 0 : inventoryId.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Inventory other = (Inventory) obj;
		if (inventoryId == null) {
			if (other.inventoryId != null)
				return false;
		} else if (!inventoryId.equals(other.inventoryId))
			return false;
		return true;
	}
}