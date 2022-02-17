import java.util.*;

/**
 * A class to store all relevant information about the supply chain.
 *
 */
public class SC {
	private List<Inventory> inventories;
	private List<Order> orders;
	private List<Material> materials;
	private List<Customer> customers;
	private List<String> arcs;
	private TreeSet<Date> dates;
	private Date splitDate;
	private boolean isDemand;
	private List<String> varNamesS; // 3 variables: demand, inventory position, order quantity --> for stage type DIST
	private List<String> varNamesM; // 4 variables: demand, inventory position, lead time, order quantity --> for stage type WH
	private List<String> varNamesL; // 5 variables: demand, inventory position, lead time, market price, order quantity --> for stage type PLANT
	private List<String> labels3; // List with 3 labels: low, medium, high
	private List<String> labels5; // List with 5 labels: very low, low, medium, high, very high
	private Map<String,WMGeneration> wmgs;
	private Map<String,Double> maxScCost;
	private DataGeneration dataGenerator;
	
	// Maps that contain the fuzzy input data for each inventory-material
	private Map<String,TreeMap<Date,Integer>> demandsMap;
	private Map<String,TreeMap<Date,Integer>> demandChangesMap;
	private Map<String,TreeMap<Date,Integer>> inventoryLevelsMap;
	private Map<String,TreeMap<Date,Integer>> inventoryPositionsMap;
	private Map<String,TreeMap<Date,Integer>> leadTimesMap;
	private Map<String,TreeMap<Date,Double>> marketPricesMap;
	private Map<String,TreeMap<Date,Integer>> orderQuantities;
	private Map<String,List<TreeMap<Date,Integer>>> scSubDemands;
	
	private HashMap<String,Integer> initInvLevels;
	private boolean isSimulation;
	// Maps that contain info that is passed from training to test pahse
	private Map<String,Integer> finalInventoryPositions;
	private Map<String,Integer> finalInventoryLevels;
	private Map<String,TreeMap<Date,Integer>> ordersToArriveInTestPhase;
	
	private List<String> tuples; // List of all inventory-material identifiers
	private List<Inventory> sortedInventories; // List of all active inventories sorted in the order DIST - WH - P
	
	/**
	 * Default constructor.
	 * @param inventories List of all inventories in supply chain
	 * @param orders List of all orders
	 * @param materials List of all materials
	 * @param customers List of all customers
	 * @param arcs List of all arcs of the network
	 * @param dates List of all dates
	 * @param splitDate Last date of training phase
	 * @param isDemand Indicates whether absolute demand or demand change used
	 */
	public SC(List<Inventory> inventories, List<Order> orders, List<Material> materials, List<Customer> customers, 
			List<String> arcs, TreeSet<Date> dates, Date splitDate, boolean isDemand) {
		this.inventories = inventories;
		this.orders = orders;
		this.materials = materials;
		this.customers = customers;
		this.arcs = arcs;
		this.dates = dates;
		this.splitDate = splitDate;
		this.isDemand = isDemand;
		// Initialize variable names
		if(isDemand) {
			this.varNamesS = Arrays.asList("demand", "inventoryPosition", /*"supplierReliability",*/ "orderQuantity");
			this.varNamesM = Arrays.asList("demand", "inventoryPosition", "leadTime", /*"supplierReliability",*/ "orderQuantity");
			this.varNamesL = Arrays.asList("demand", "inventoryPosition", "leadTime", /*"supplierReliability",*/ "marketPrice", "orderQuantity");
		} else {
			this.varNamesS = Arrays.asList("demandChange", "inventoryPosition", /*"supplierReliability",*/ "orderQuantity");
			this.varNamesM = Arrays.asList("demandChange", "inventoryPosition", "leadTime", /*"supplierReliability",*/ "orderQuantity");
			this.varNamesL = Arrays.asList("demandChange", "inventoryPosition", "leadTime", /*"supplierReliability",*/ "marketPrice", "orderQuantity");
		}
		// Initialize label names
		this.labels3 = Arrays.asList("low", "medium", "high");
		this.labels5 = Arrays.asList("veryLow", "low", "medium", "high", "veryHigh");
		this.wmgs = new HashMap<String,WMGeneration>();
		this.maxScCost = new HashMap<>();
		this.demandsMap = new HashMap<>();
		this.demandChangesMap = new HashMap<>();
		this.inventoryLevelsMap = new HashMap<>();
		this.inventoryPositionsMap = new HashMap<>();
		this.leadTimesMap = new HashMap<>();
		this.marketPricesMap = new HashMap<>();
		this.orderQuantities = new HashMap<>();
		this.scSubDemands = new HashMap<>();
		this.initInvLevels = new HashMap<>();
		this.isSimulation = false;
		this.finalInventoryPositions = new HashMap<>();
		this.finalInventoryLevels = new HashMap<>();
		this.ordersToArriveInTestPhase = new HashMap<>();
		this.tuples = new ArrayList<String>();
		this.sortedInventories = new ArrayList<Inventory>();
		// Create all tuples and inventories in the proper order
		// Start with distribution centers
		for(Inventory inventory : inventories) {
			if(inventory.getStageType().equals("DIST")) {
				for(Material material : inventory.getMaterials()) {
					String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
					tuples.add(tuple);
				}
				sortedInventories.add(inventory);
			}
		}
		// 2) warehouses
		for(Inventory inventory : inventories) {
			if(inventory.getStageType().equals("WH")) {
				for(Material material : inventory.getMaterials()) {
					String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
					tuples.add(tuple);
				}
				sortedInventories.add(inventory);
			}
		}
		// 3) plants
		for(Inventory inventory : inventories) {
			if(inventory.getStageType().equals("PLANT")) {
				for(Material material : inventory.getMaterials()) {
					String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
					tuples.add(tuple);
				}
				sortedInventories.add(inventory);
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
	
	/**
	 * Rounds a given value to 2 decimal points.
	 * @param value
	 * @return rounded value
	 */
	public double roundDecimal(double value) {
		return (double)Math.round(value * 100d) / 100d;
	}

	/**
	 * @return the inventories (unsorted)
	 */
	public List<Inventory> getInventories() {
		return inventories;
	}
	
	/**
	 * @param newInventories the inventories to set
	 */
	public void setInventories(List<Inventory> newInventories) {
		this.inventories = newInventories;
	}

	/**
	 * @return the orders
	 */
	public List<Order> getOrders() {
		return orders;
	}

	/**
	 * @param orders the orders to set
	 */
	public void setOrders(List<Order> orders) {
		this.orders = orders;
	}

	/**
	 * @return the materials
	 */
	public List<Material> getMaterials() {
		return materials;
	}

	/**
	 * @param materials the materials to set
	 */
	public void setMaterials(List<Material> materials) {
		this.materials = materials;
	}
	
	/**
	 * @return the customers
	 */
	public List<Customer> getCustomers() {
		return customers;
	}
	
	/**
	 * @return the arcs
	 */
	public List<String> getArcs() {
		return arcs;
	}

	/**
	 * @return the dates
	 */
	public TreeSet<Date> getDates() {
		return dates;
	}

	/**
	 * @return the splitDate
	 */
	public Date getSplitDate() {
		return splitDate;
	}
	
	/**
	 * @return the varNamesS
	 */
	public List<String> getVarNamesS() {
		return varNamesS;
	}
	
	/**
	 * @return the varNamesM
	 */
	public List<String> getVarNamesM() {
		return varNamesM;
	}

	/**
	 * @return the varNamesL
	 */
	public List<String> getVarNamesL() {
		return varNamesL;
	}

	/**
	 * @return the labels3
	 */
	public List<String> getLabels3() {
		return labels3;
	}

	/**
	 * @return the labels5
	 */
	public List<String> getLabels5() {
		return labels5;
	}
	
	/**
	 * @return the isDemand
	 */
	public boolean isDemand() {
		return isDemand;
	}

	/**
	 * @return the wmgs
	 */
	public Map<String, WMGeneration> getWmgs() {
		return wmgs;
	}
	
	/**
	 * @return the maxScCost
	 */
	public Map<String, Double> getMaxScCost() {
		return maxScCost;
	}

	/**
	 * @return set of training dates
	 */
	public TreeSet<Date> getTrainingDates() {
		TreeSet<Date> trainingDates = new TreeSet<Date>();
		for(Date date : dates) {
			if(!date.after(splitDate)) {
				trainingDates.add(date);
			}
		}
		return trainingDates;
	}
	
	/**
	 * @return set of test dates
	 */
	public TreeSet<Date> getTestDates() {
		TreeSet<Date> testDates = new TreeSet<Date>();
		for(Date date : dates) {
			if(date.after(splitDate)) {
				testDates.add(date);
			}
		}
		return testDates;
	}

	/**
	 * @return the demandsMap
	 */
	public Map<String, TreeMap<Date, Integer>> getDemandsMap() {
		return demandsMap;
	}

	/**
	 * @param demandsMap the demandsMap to set
	 */
	public void setDemandsMap(Map<String, TreeMap<Date, Integer>> demandsMap) {
		this.demandsMap = demandsMap;
	}

	/**
	 * @return the demandChangesMap
	 */
	public Map<String, TreeMap<Date, Integer>> getDemandChangesMap() {
		return demandChangesMap;
	}

	/**
	 * @param demandChangesMap the demandChangesMap to set
	 */
	public void setDemandChangesMap(Map<String, TreeMap<Date, Integer>> demandChangesMap) {
		this.demandChangesMap = demandChangesMap;
	}

	/**
	 * @return the inventoryLevelsMap
	 */
	public Map<String, TreeMap<Date, Integer>> getInventoryLevelsMap() {
		return inventoryLevelsMap;
	}

	/**
	 * @param inventoryLevelsMap the inventoryLevelsMap to set
	 */
	public void setInventoryLevelsMap(Map<String, TreeMap<Date, Integer>> inventoryLevelsMap) {
		this.inventoryLevelsMap = inventoryLevelsMap;
	}

	/**
	 * @return the inventoryPositionsMap
	 */
	public Map<String, TreeMap<Date, Integer>> getInventoryPositionsMap() {
		return inventoryPositionsMap;
	}

	/**
	 * @param inventoryPositionsMap the inventoryPositionsMap to set
	 */
	public void setInventoryPositionsMap(Map<String, TreeMap<Date, Integer>> inventoryPositionsMap) {
		this.inventoryPositionsMap = inventoryPositionsMap;
	}

	/**
	 * @return the leadTimesMap
	 */
	public Map<String, TreeMap<Date, Integer>> getLeadTimesMap() {
		return leadTimesMap;
	}

	/**
	 * @param leadTimesMap the leadTimesMap to set
	 */
	public void setLeadTimesMap(Map<String, TreeMap<Date, Integer>> leadTimesMap) {
		this.leadTimesMap = leadTimesMap;
	}

	/**
	 * @return the marketPricesMap
	 */
	public Map<String, TreeMap<Date, Double>> getMarketPricesMap() {
		return marketPricesMap;
	}

	/**
	 * @param marketPricesMap the marketPricesMap to set
	 */
	public void setMarketPricesMap(Map<String, TreeMap<Date, Double>> marketPricesMap) {
		this.marketPricesMap = marketPricesMap;
	}

	/**
	 * @return the orderQuantities
	 */
	public Map<String, TreeMap<Date, Integer>> getOrderQuantities() {
		return orderQuantities;
	}

	/**
	 * @param orderQuantities the orderQuantities to set
	 */
	public void setOrderQuantities(Map<String, TreeMap<Date, Integer>> orderQuantities) {
		this.orderQuantities = orderQuantities;
	}

	/**
	 * @return the scSubDemands
	 */
	public Map<String, List<TreeMap<Date, Integer>>> getScSubDemands() {
		return scSubDemands;
	}

	/**
	 * @param scSubDemands the scSubDemands to set
	 */
	public void setScSubDemands(Map<String, List<TreeMap<Date, Integer>>> scSubDemands) {
		this.scSubDemands = scSubDemands;
	}

	/**
	 * @return the initInvLevels
	 */
	public HashMap<String, Integer> getInitInvLevels() {
		return initInvLevels;
	}

	/**
	 * @param initInvLevels the initInvLevels to set
	 */
	public void setInitInvLevels(HashMap<String, Integer> initInvLevels) {
		this.initInvLevels = initInvLevels;
	}

	/**
	 * @return the isSimulation
	 */
	public boolean isSimulation() {
		return isSimulation;
	}

	/**
	 * @param isSimulation the isSimulation to set
	 */
	public void setIsSimulation(boolean isSimulation) {
		this.isSimulation = isSimulation;
	}

	/**
	 * @return the finalInventoryPositions
	 */
	public Map<String, Integer> getFinalInventoryPositions() {
		return finalInventoryPositions;
	}

	/**
	 * @param finalInventoryPositions the finalInventoryPositions to set
	 */
	public void setFinalInventoryPositions(Map<String, Integer> finalInventoryPositions) {
		this.finalInventoryPositions = finalInventoryPositions;
	}

	/**
	 * @return the finalInventoryLevels
	 */
	public Map<String, Integer> getFinalInventoryLevels() {
		return finalInventoryLevels;
	}

	/**
	 * @param finalInventoryLevels the finalInventoryLevels to set
	 */
	public void setFinalInventoryLevels(Map<String, Integer> finalInventoryLevels) {
		this.finalInventoryLevels = finalInventoryLevels;
	}

	/**
	 * @return the ordersToArriveInTestPhase
	 */
	public Map<String, TreeMap<Date, Integer>> getOrdersToArriveInTestPhase() {
		return ordersToArriveInTestPhase;
	}

	/**
	 * @param ordersToArriveInTestPhase the ordersToArriveInTestPhase to set
	 */
	public void setOrdersToArriveInTestPhase(Map<String, TreeMap<Date, Integer>> ordersToArriveInTestPhase) {
		this.ordersToArriveInTestPhase = ordersToArriveInTestPhase;
	}

	/**
	 * @return the tuples
	 */
	public List<String> getTuples() {
		return tuples;
	}

	/**
	 * @return the sortedInventories
	 */
	public List<Inventory> getSortedInventories() {
		return sortedInventories;
	}

	/**
	 * @return the dataGenerator
	 */
	public DataGeneration getDataGenerator() {
		return dataGenerator;
	}

	/**
	 * @param dataGenerator the dataGenerator to set
	 */
	public void setDataGenerator(DataGeneration dataGenerator) {
		this.dataGenerator = dataGenerator;
	}
}