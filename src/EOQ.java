import java.util.*;
import java.util.Date;
import java.util.TreeMap;

/**
 * A class to determine the optimal batch quantity and re-oroder level for a particular inventory and material. 
 * The optimal batch quantity is based on the Economic Order Quantity (EOQ). 
 * The optimal re-order level considers a safety stock that assumes Normal demand and a Gamma distributed lead time.
 * 2 methods offered: 1) Normal distribution with uncertainty on the lead time, 2) Average-Max formula
 */
public class EOQ {
	private SC sc;
	private Inventory inventory;
	private Material material;
	private TreeMap<Date,Integer> demands;
	private TreeMap<Date,Integer> leadTimes;
	
	private double serviceLevel;
	private boolean avgMax; // Indicates if method 1 or 2 is applied
	// For method 1
	private double k; // Safety factor (Normal distribution of demand)
	private double stdvLt;
	private double avgDemand;
	private double avgLt;
	// For method 2
	private int maxDemand;
	private int maxLt;
	
	private int orderQnt; // Q
	private int safetyStock; // SS
	private int reorderLevel; // s
	private int orderUpToLevel; // S
	
	private double holdingCosts; // h
	private double orderCosts; // o
	
	
	/** Default constructor.
	 * @param sc SC object
	 * @param inventory Inventory of interest
	 * @param material Material of interest
	 * @param demands Demands observed up to this date
	 * @param leadTimes Lead times observed up to this date
	 * @param serviceLevel Service level S1 (probability of no stock out) - Choose out of {0.9, 0.95, 0.98, 0.99}
	 */
	public EOQ(SC sc, Inventory inventory, Material material, TreeMap<Date,Integer> demands, TreeMap<Date,Integer> leadTimes, double serviceLevel, boolean avgMax) {
		this.inventory = inventory;
		this.material = material;
		String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode(); // Identifier of inventory and material
		this.demands = demands;
		this.leadTimes = leadTimes;
		this.serviceLevel = serviceLevel;
		this.avgMax = avgMax;
		
		// Define safety factor k
		if(serviceLevel == 0.9) {
			this.k = 1.28;
		}
		if(serviceLevel == 0.95) {
			this.k = 1.64;
		}		
		if(serviceLevel == 0.98) {
			this.k = 2.05;
		}		
		if(serviceLevel == 0.99) {
			this.k = 2.33;
		}
		
		// Compute average demand
		if(sc.isSimulation()) {
			this.avgDemand = this.getRealAverageDemand();
		} else {
			this.avgDemand = this.getAverageDemand();
		}
		// Compute average lead time
		this.avgLt = this.getAverageLeadTime();
		// Compute max lead time
		this.maxLt = (Collections.max(leadTimes.values())) / 24;
		// Compute standard deviation of lead time
		this.stdvLt = this.getStdvLeadTime();
		// Compute max demand
		this.maxDemand = (Collections.max(demands.values()));
		
		// Retrieve holding and order costs
		this.holdingCosts = inventory.getHoldingCosts().get(material);
		if(inventory.getStageType().equals("WH")) {
			this.orderCosts = inventory.getSetupCosts().get(material);
		} else {
			this.orderCosts = inventory.getOrderCosts().get(material);
		}
				
		// Compute SS
		if(avgMax) {
			// Method 2
			this.safetyStock = (int) Math.ceil(maxDemand*maxLt - avgDemand * avgLt);
		} else {
			// Method 1
			this.safetyStock = (int) Math.ceil(k * stdvLt * avgDemand);
		}
		// Compute s
		this.reorderLevel = (int) Math.ceil(safetyStock + avgDemand * avgLt);
		// Compute Q
		this.orderQnt = (int) Math.ceil(Math.sqrt((2*orderCosts*avgDemand) / holdingCosts));
		// Compute S
		this.orderUpToLevel = reorderLevel + orderQnt;
	}
	
	/**
	 * @return Standard deviation of the lead time
	 */
	public double getStdvLeadTime() {
		double sum = 0.0;
		for (var entry : leadTimes.entrySet()) {
		    sum += Math.pow(entry.getValue()/24 - avgLt, 2);
		}
		return Math.sqrt(sum / (leadTimes.size() - 1));
	}
	
	/** Returns the average demand per time unit by only considering periods with positive demand.
	 * @param demandsMap
	 */
	public double getAverageDemand() {
		int sum = 0;
		int counter = 0;
		for (var entry : demands.entrySet()) {
		    sum += entry.getValue();
		    if(entry.getValue() > 0) {
		    	counter += 1;
		    }
		}
		if(counter == 0) {
			return 0.0;
		} else {
			return (double) sum / counter;
		}
	}
	
	/** Returns the average demand per time unit.
	 * @param demandsMap
	 */
	public double getRealAverageDemand() {
		int sum = 0;
		for (var entry : demands.entrySet()) {
		    sum += entry.getValue();
		}
		return (double) sum / demands.size();
	}
	
	/** Returns the average lead time per time unit.
	 * @param demandsMap
	 */
	public double getAverageLeadTime() {
		int sum = 0;
		for (var entry : leadTimes.entrySet()) {
		    sum += entry.getValue()/24;
		}
		return (double) sum / leadTimes.size();
	}

	/**
	 * @return the orderQnt
	 */
	public int getOrderQnt() {
		return orderQnt;
	}

	/**
	 * @return the safetyStock
	 */
	public int getSafetyStock() {
		return safetyStock;
	}

	/**
	 * @return the reorderLevel
	 */
	public int getReorderLevel() {
		return reorderLevel;
	}

	/**
	 * @return the orderUpToLevel
	 */
	public int getOrderUpToLevel() {
		return orderUpToLevel;
	}
}