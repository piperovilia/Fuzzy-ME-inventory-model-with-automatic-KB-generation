import java.util.*;

import net.sourceforge.jFuzzyLogic.FIS;
import net.sourceforge.jFuzzyLogic.FunctionBlock;
import net.sourceforge.jFuzzyLogic.plot.JFuzzyChart;
import net.sourceforge.jFuzzyLogic.rule.Variable;

/**
 * CORE OF PROCEDURE: Every time the performance of knowledge base (KB) or a set of knowledge bases (KB set) is evaluated, a complete time-step simulation of an inventory or the entire SC must be performed.
 * This class offers methods to simulate an inventory or a complete supply chain by applying a particular inventory policy (EOQ, Fuzzy rule-based system) over a time period of interest (training OR test).
 */
public class Simulation {
	private SC sc;
	private List<Inventory> inventories;
	private List<Inventory> plants;
	private List<Inventory> warehouses;
	private List<Inventory> distributionCenters;
	private Map<String,List<TreeMap<Date,Integer>>> subDemandsMaps;
	private TreeSet<Date> dates; // All dates (training + test)
	private Date splitDate; // Last date of training phase
	// The following maps have a key "inventory,material" which identifies a particular inventory and material
	private Map<String,TreeMap<Date,Integer>> orderQuantitiesFIS; // String key = "inventory,material"; This map contains order quantities place by an inventory (inventory-material combination)
	private Map<String,Integer> numTotalOrders; // Contains total number of placed orders
	private Map<String,Integer> numUnmetOrders; // Contains total number of unsatisfied orders
	private Map<String,Integer> totalDemands;
	private Map<String,Double> totalCosts;
	private Map<String,Double> holdingCostMap;
	private Map<String,Double> penaltyCostMap;
	private Map<String,Double> orderCostMap;
	private Map<String,Double> setupCostMap;
	private Map<String,Double> productionCostMap;
	private Map<String,Double> unitCostMap;
	private Map<String,Double> transportationCostMap;
	private Map<String,TreeMap<Date,Integer>> inventoryLevelsMap;
	private Map<String,TreeMap<Date,Integer>> inventoryPositionsMap;
	private int iteration;
    private long startTime;
    private long runTime;
	
	/** Default constructor.
	 * @param sc The supply chain
	 * @param iteration
	 */
	public Simulation(SC sc, int iteration) {
		this.sc = sc;
		this.inventories = sc.getInventories();
		this.plants = new ArrayList<Inventory>();
		this.warehouses = new ArrayList<Inventory>();
		this.distributionCenters = new ArrayList<Inventory>();
		// Create lists of plants, warehouses and distribution centers
		for(Inventory i : inventories) {
			if(i.getStageType().equals("PLANT")) {
				plants.add(i);
			}
			if(i.getStageType().equals("WH")) {
				warehouses.add(i);
			}
			if(i.getStageType().equals("DIST")) {
				distributionCenters.add(i);
			}
		}
		this.subDemandsMaps = new HashMap<>();
		// Initialize sub-demands
		for(Inventory i : sc.getSortedInventories()) {
			for(Material m : i.getMaterials()) {
				String tuple = "" + i.getInventoryId() + "," + m.getMaterialCode();
				if(i.getStageType().equals("DIST")) {
					// If a distribution center --> take sub-demands that were created by the DataGeneration class
					subDemandsMaps.put(tuple, sc.getScSubDemands().get(tuple));
				} else {
					// If not DC --> initialize with empty maps
					List<TreeMap<Date,Integer>> emptyList = new ArrayList<TreeMap<Date,Integer>>();
					for(int j=0; j<i.getSuccessors().size(); j++) {
						TreeMap<Date,Integer> emptyMap = new TreeMap<Date,Integer>();
						for(Date d : sc.getDates()) {
							emptyMap.put(d, 0);
						}
						emptyList.add(emptyMap);
					}
					subDemandsMaps.put(tuple, emptyList);
				}
			}
		}	
		this.dates = sc.getDates();
		this.splitDate = sc.getSplitDate();
		this.orderQuantitiesFIS = new HashMap<String,TreeMap<Date,Integer>>();
		this.numTotalOrders = new HashMap<String,Integer>();
		this.numUnmetOrders = new HashMap<String,Integer>();
		this.totalDemands = new HashMap<String,Integer>();
		this.totalCosts = new HashMap<String,Double>();
		this.holdingCostMap = new HashMap<String,Double>();
		this.penaltyCostMap = new HashMap<String,Double>();
		this.orderCostMap = new HashMap<String,Double>();
		this.setupCostMap = new HashMap<String,Double>();
		this.productionCostMap = new HashMap<String,Double>();
		this.unitCostMap = new HashMap<String,Double>();
		this.transportationCostMap = new HashMap<String,Double>();
		this.inventoryLevelsMap = new HashMap<>();
		this.inventoryPositionsMap = new HashMap<>();
		// Initialize result maps with 0's
		for(Inventory inventory : sc.getSortedInventories()) {
			for(Material material : sc.getMaterials()) {
				String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
				numTotalOrders.put(tuple, 0);
				numUnmetOrders.put(tuple, 0);
				totalDemands.put(tuple, 0);
				totalCosts.put(tuple, 0.0);
				holdingCostMap.put(tuple, 0.0);
				penaltyCostMap.put(tuple, 0.0);
				orderCostMap.put(tuple, 0.0);
				setupCostMap.put(tuple, 0.0);
				productionCostMap.put(tuple, 0.0);
				unitCostMap.put(tuple, 0.0);
				transportationCostMap.put(tuple, 0.0);
				// Initialize with empty maps
				TreeMap<Date,Integer> emptyMap = new TreeMap<Date,Integer>();
				orderQuantitiesFIS.put(tuple, emptyMap);
				TreeMap<Date,Integer> emptyMapIl = new TreeMap<Date,Integer>();
				this.inventoryLevelsMap.put(tuple, emptyMapIl);
				TreeMap<Date,Integer> emptyMapIp = new TreeMap<Date,Integer>();
				this.inventoryPositionsMap.put(tuple, emptyMapIp);
			}
		}
		this.iteration = iteration;
	}
	
	/**
	 * Returns if there is enough raw material on stock at 'predeessorStr' to fulfill a given 'quantity' of an end-product 'endProduct' that is ordered by an inventory-material on date 'd'.
	 * The method directly adds the order to the demands of the raw material predecessors if enough stock is available to fulfill the whole order.
	 * @param inventoryLevels Contains inventory levels of all inventory and materials on date d
	 * @param inventory
	 * @param endProduct
	 * @param quantity Order quantity
	 * @return true if enough raw material available and false if not
	 */
	public boolean isSufficientRawMaterial(Map<String,Integer> inventoryLevels, Inventory inventory, Material endProduct, int quantity, String predecessorStr, Date d, Map<String,TreeMap<Date,Integer>> demands) {
		boolean sufficient = true;
		// Add new order to sub demand of predecessor
		for(Material matROH : endProduct.getIngredients().keySet()) {
			String tuple = "" + inventory.getInventoryId() + "," + matROH.getMaterialCode();
			// Retrieve the order quantity of the raw material
			double amount;
			if(matROH.getMaterialCode().equals("CC-P01") || matROH.getMaterialCode().equals("CC-P02") || matROH.getMaterialCode().equals("CC-P03") || matROH.getMaterialCode().equals("CC-P04")) {
				amount = quantity * endProduct.getIngredients().get(matROH);
			} else {
				amount = quantity * endProduct.getIngredients().get(matROH) * endProduct.getUnitSize();
			}	
			int qntRawMaterial = (int) Math.ceil(amount);
			String rawPredecessorTuple = "" + predecessorStr + "," + matROH.getMaterialCode();
			// Add new order to sub demand of predecessor
			demands.get(rawPredecessorTuple).put(d, qntRawMaterial);
			if(inventoryLevels.get(tuple) < qntRawMaterial) {
				sufficient = false;
			}
		}
		return sufficient;
	}
	
	/**
	 * SIMULATION OF ENTIRE SC BY APPLYING FRBS
	 * This method simulates all inventories of the supply chain over a time horizon (training OR test phase) by applying Fuzzy Sets and Fuzzy Logic (Fuzzy rule-based systems).
	 * Computes KPIs, inventory levels, inventory positions, etc.
	 * @param training Indicates whether the simulation is conducted for the training (true) or test (false) period
	 * @throws Exception
	 */
	public void simulateSupplyChain(boolean training) throws Exception {
		//long startTime = System.nanoTime();
		Map<String,Integer> inventoryLevels = new HashMap<>();
		Map<String,Integer> inventoryPositions = new HashMap<>();
		Map<String,TreeMap<Date,Integer>> ordersToArriveMap = new HashMap<>();
		for(String tuple : sc.getTuples()) {
			if(training) {
				inventoryLevels.put(tuple, sc.getInitInvLevels().get(tuple));
				inventoryPositions.put(tuple, sc.getInitInvLevels().get(tuple));
				TreeMap<Date,Integer> emptyMap = new TreeMap<Date,Integer>();
				for(Date d : sc.getDates()) {
					emptyMap.put(d, 0);
				}
				ordersToArriveMap.put(tuple, emptyMap);
			} else {
				if(sc.isSimulation()) {
					// Start SC from scratch
					inventoryLevels.put(tuple, sc.getInitInvLevels().get(tuple));
					inventoryPositions.put(tuple, sc.getInitInvLevels().get(tuple));
					TreeMap<Date,Integer> emptyMap = new TreeMap<Date,Integer>();
					for(Date d : sc.getDates()) {
						emptyMap.put(d, 0);
					}
					ordersToArriveMap.put(tuple, emptyMap);
				} else {
					// Take over data from training stage
					inventoryLevels.put(tuple,  sc.getFinalInventoryLevels().get(tuple));
					inventoryPositions.put(tuple, sc.getFinalInventoryPositions().get(tuple));
					TreeMap<Date,Integer> ordersToArriveCopy = new TreeMap<Date,Integer>(sc.getOrdersToArriveInTestPhase().get(tuple));
					ordersToArriveMap.put(tuple, ordersToArriveCopy);
				}
			}
		}
		// Keep track of sub-demands at each inventory
		Map<String,List<TreeMap<Date,Integer>>> subDemandMap = new HashMap<>();
		Map<String,TreeMap<Date,Integer>> demands = new HashMap<>();
		for(Inventory i : sc.getSortedInventories()) {
			for(Material m : i.getMaterials()) {
				String tuple = "" + i.getInventoryId() + "," + m.getMaterialCode();
				if(i.getStageType().equals("DIST")) {
					List<TreeMap<Date,Integer>> subDemandCopy = new ArrayList<TreeMap<Date,Integer>>(sc.getScSubDemands().get(tuple));
					subDemandMap.put(tuple, subDemandCopy);
					TreeMap<Date,Integer> demandsCopy = new TreeMap<Date,Integer>(sc.getDemandsMap().get(tuple));
					demands.put(tuple, demandsCopy);
				} else {
					List<TreeMap<Date,Integer>> emptyList = new ArrayList<TreeMap<Date,Integer>>();
					for(int j=0; j<i.getSuccessors().size(); j++) {
						TreeMap<Date,Integer> emptyMap = new TreeMap<Date,Integer>();
						for(Date d : sc.getDates()) {
							emptyMap.put(d, 0);
						}
						emptyList.add(emptyMap);
					}
					subDemandMap.put(tuple, emptyList);
					TreeMap<Date,Integer> emptyMap = new TreeMap<Date,Integer>();
					for(Date d : sc.getDates()) {
						emptyMap.put(d, 0);
					}
					demands.put(tuple, emptyMap);
				}
			}
		}
		
		// Obtain time period of interest
		TreeSet<Date> datesSim;
		if(training) {
			datesSim = sc.getTrainingDates();
		} else {
			datesSim = sc.getTestDates();
		}
		for(Date d : datesSim) {
			// 1) Create demand along the whole chain
			for(Inventory inventory : sc.getSortedInventories()) { 
				for(Material material : inventory.getMaterials()) {
					String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();		
					// Update inventory level and inventory positions
					inventoryLevels.put(tuple, inventoryLevels.get(tuple) + ordersToArriveMap.get(tuple).get(d));
					// Add transport costs for the suborders that can be fulfilled
					for(int i=0; i<subDemandMap.get(tuple).size(); i++) {
						TreeMap<Date,Integer> subDemands = subDemandMap.get(tuple).get(i);
						String successor = inventory.getSuccessors().get(i);
						int subDemand = subDemands.get(d);
						if(subDemand > 0) {
							if(inventory.getStageType().equals("DIST")) {
								numTotalOrders.put(tuple, numTotalOrders.get(tuple) + 1); 
							}	
							if(inventoryLevels.get(tuple) >= subDemand) {
								// Add transportation costs
								if(!inventory.getStageType().equals("PLANT")) {
									transportationCostMap.put(tuple, transportationCostMap.get(tuple) + inventory.getTransportCosts().get(successor).get(material) * subDemand);
								}
								// Sub-order CAN be satisfied with current inventory
								inventoryLevels.put(tuple, inventoryLevels.get(tuple) - subDemand);
								inventoryPositions.put(tuple, inventoryPositions.get(tuple) - subDemand);
							} else {
								// Add penalty costs and increase counter of not-satisfied orders
								if(inventory.getStageType().equals("DIST")) {
									penaltyCostMap.put(tuple, penaltyCostMap.get(tuple) + inventory.getShortageCosts().get(material));
									// Only customer orders matter for the number of unmet orders
									numUnmetOrders.put(tuple, numUnmetOrders.get(tuple) + 1);
								}
							}
						}
					}
					int demand = demands.get(tuple).get(d);
					// Input 3: Lead time (in hours) --> only fuzzy if NOT a distribution center, else 0
					int leadTime = this.getLeadTimeForeCastAvg(inventory, material, d);
					// Input 4: Market price
					double marketPrice = 0.0;
					if(inventory.getStageType().equals("PLANT")) {
						marketPrice = sc.getMarketPricesMap().get(tuple).get(d);
					}
					
					// Retrieve optimal order quantity by applying FIS model
					String predecessorStr = inventory.getPredecessor();
					Inventory predecessor = this.getInventoryByCode(predecessorStr);
					int indexSuccessor = predecessor.getSuccessors().indexOf(inventory.getInventoryId());
					String predecessorTuple = predecessorStr + "," + material.getMaterialCode();
					
					if(inventory.getStageType().equals("DIST")) {				
						double[] inputs = {demand, inventoryPositions.get(tuple)};
						int orderQnt = this.getOrderQuantity(inventory, material, inputs, sc.getVarNamesS()); // 3 labels applied
						demands.get(predecessorTuple).put(d, demands.get(predecessorTuple).get(d) + orderQnt);	
						if(orderQnt > 0) {	
							numTotalOrders.put(predecessorTuple, numTotalOrders.get(predecessorTuple) + 1);
							if(inventoryLevels.get(predecessorTuple) + ordersToArriveMap.get(predecessorTuple).get(d) >= orderQnt) {
								Date newDate;
								if(sc.isSimulation()) {
									newDate = this.addDays(d, (int) sc.getLeadTimesMap().get(tuple).get(d) / 24);
								} else {
									newDate = this.addWorkingDays(d, (int) sc.getLeadTimesMap().get(tuple).get(d) / 24);
								}
								if(newDate.equals(d)) {
									inventoryLevels.put(tuple, inventoryLevels.get(tuple) + orderQnt);
									inventoryPositions.put(tuple, inventoryPositions.get(tuple) + orderQnt);
								} else {
									if(!newDate.after(sc.getDates().last())) {
										ordersToArriveMap.get(tuple).put(newDate, ordersToArriveMap.get(tuple).get(newDate) + orderQnt);
									}		
									inventoryPositions.put(tuple, inventoryPositions.get(tuple) + orderQnt);
								}
								// Add new order to sub demand of predecessor
								int currentQnt = subDemandMap.get(predecessorTuple).get(indexSuccessor).get(d);
								subDemandMap.get(predecessorTuple).get(indexSuccessor).put(d, currentQnt + orderQnt);
								// Add order costs
								orderCostMap.put(tuple, orderCostMap.get(tuple) + inventory.getOrderCosts().get(material));
							} else {
								// Add penalty costs and increase counter of not-satisfied orders
								penaltyCostMap.put(predecessorTuple, penaltyCostMap.get(predecessorTuple) + predecessor.getShortageCosts().get(material));
								// Only customer orders matter for the number of unmet orders
								numUnmetOrders.put(predecessorTuple, numUnmetOrders.get(predecessorTuple) + 1);
							}
						}
					}
					
					if(inventory.getStageType().equals("WH")) {
						double[] inputs = {demand, inventoryPositions.get(tuple), leadTime};
						int orderQnt = this.getOrderQuantity(inventory, material, inputs, sc.getVarNamesM()); // 3 labels applied
						if(orderQnt > 0) {					
							numTotalOrders.put(predecessorTuple, numTotalOrders.get(predecessorTuple) + 1);
							if(this.isSufficientRawMaterial(inventoryLevels, predecessor, material, orderQnt, predecessorStr, d, demands)) {
								Date newDate;
								if(sc.isSimulation()) {
									newDate = this.addDays(d, (int) sc.getLeadTimesMap().get(tuple).get(d) / 24);
								} else {
									newDate = this.addWorkingDays(d, (int) sc.getLeadTimesMap().get(tuple).get(d) / 24);
								}
								if(newDate.equals(d)) {
									inventoryLevels.put(tuple, inventoryLevels.get(tuple) + orderQnt);
									inventoryPositions.put(tuple, inventoryPositions.get(tuple) + orderQnt);
								} else {
									if(!newDate.after(sc.getDates().last())) {
										ordersToArriveMap.get(tuple).put(newDate, ordersToArriveMap.get(tuple).get(newDate) + orderQnt);
									}				
									inventoryPositions.put(tuple, inventoryPositions.get(tuple) + orderQnt);
								}
								for(var entry : material.getIngredients().entrySet()) {
									Material matROH = entry.getKey();
									double amount;
									if(matROH.getMaterialCode().equals("CC-P01") || matROH.getMaterialCode().equals("CC-P02") || matROH.getMaterialCode().equals("CC-P03") || matROH.getMaterialCode().equals("CC-P04")) {
										amount = orderQnt * material.getIngredients().get(matROH);
									} else {
										amount = orderQnt * material.getIngredients().get(matROH) * material.getUnitSize();
									}
									
									int qntRawMaterial = (int) Math.ceil(amount);
									String rawPredecessorTuple = "" + predecessorStr + "," + matROH.getMaterialCode();
									// Add new order to sub demand of predecessor
									int currentQnt = subDemandMap.get(rawPredecessorTuple).get(indexSuccessor).get(d);
									subDemandMap.get(rawPredecessorTuple).get(indexSuccessor).put(d, currentQnt + qntRawMaterial);
								}
								setupCostMap.put(tuple, setupCostMap.get(tuple) + inventory.getSetupCosts().get(material));
								productionCostMap.put(tuple, productionCostMap.get(tuple) + orderQnt * inventory.getProdUnitCosts().get(material));
							} else {
								// Add penalty costs and increase counter of not-satisfied orders
								penaltyCostMap.put(predecessorTuple, penaltyCostMap.get(predecessorTuple) + predecessor.getShortageCosts().get(material));
								// Only customer orders matter for the number of unmet orders
								numUnmetOrders.put(predecessorTuple, numUnmetOrders.get(predecessorTuple) + 1);
							}
						}
					} 		
					if(inventory.getStageType().equals("PLANT")) {				
						double[] inputs = {demand, inventoryPositions.get(tuple), leadTime, marketPrice};
						int orderQnt = this.getOrderQuantity(inventory, material, inputs, sc.getVarNamesL()); // 3 labels applied
						if(orderQnt > 0) {
							// Add order costs
							orderCostMap.put(tuple, orderCostMap.get(tuple) + inventory.getOrderCosts().get(material));
							// Add purchase costs
							unitCostMap.put(tuple, unitCostMap.get(tuple) + orderQnt * marketPrice);
							Date newDate;
							if(sc.isSimulation()) {
								newDate = this.addDays(d, (int) sc.getLeadTimesMap().get(tuple).get(d) / 24);
							} else {
								newDate = this.addWorkingDays(d, (int) sc.getLeadTimesMap().get(tuple).get(d) / 24);
							}
							if(newDate.equals(d)) {
								inventoryLevels.put(tuple, inventoryLevels.get(tuple) + orderQnt);
								inventoryPositions.put(tuple, inventoryPositions.get(tuple) + orderQnt);
							} else {
								if(!newDate.after(sc.getDates().last())) {
									ordersToArriveMap.get(tuple).put(newDate, ordersToArriveMap.get(tuple).get(newDate) + orderQnt);
								}	
								inventoryPositions.put(tuple, inventoryPositions.get(tuple) + orderQnt);
							}
						}
					}
					holdingCostMap.put(tuple, holdingCostMap.get(tuple) + inventoryLevels.get(tuple) * inventory.getHoldingCosts().get(material));
					inventoryLevelsMap.get(tuple).put(d, inventoryLevels.get(tuple));
					inventoryPositionsMap.get(tuple).put(d, inventoryPositions.get(tuple));
				}
			}
			// Update total costs for this date
			for(Inventory inventory : sc.getSortedInventories()) {
				for(Material material : inventory.getMaterials()) {
					String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
					double holdingCosts = holdingCostMap.get(tuple);
					double penaltyCosts = penaltyCostMap.get(tuple);
					double orderCosts = orderCostMap.get(tuple);
					double setupCosts = setupCostMap.get(tuple);
					double productionCosts = productionCostMap.get(tuple);
					double unitCosts = unitCostMap.get(tuple);
					double transportCosts = transportationCostMap.get(tuple);
					double costsToAdd = holdingCosts + penaltyCosts + orderCosts + setupCosts + productionCosts + unitCosts + transportCosts;
					totalCosts.put(tuple, costsToAdd);
				}
			}
		}
		//long runTime = System.nanoTime() - startTime;
		//System.out.println("Runtime algorithm: " + runTime / (Math.pow(10, 6)) + " milliseconds");
	}
	
	/**
	 * SIMULATION OF ENTIRE SC BY APPLYING EOQ (ONLY TEST PHASE)
	 * This method simulates all inventories of the supply chain over the whole test horizon by applying the Economic Order Quantity (EOQ) model.
	 * @param trainingDemands demands Demands of training horizon
	 * @param trainingLeadTimes Lead times of training horizon
	 * @throws Exception
	 */
	public void simulateSupplyChainEOQ(Map<String,TreeMap<Date,Integer>> trainingDemands, Map<String,TreeMap<Date,Integer>> trainingLeadTimes) throws Exception {
		Map<String,Integer> inventoryLevels = new HashMap<>();
		Map<String,Integer> inventoryPositions = new HashMap<>();
		Map<String,TreeMap<Date,Integer>> ordersToArriveMap = new HashMap<>();
		for(String tuple : sc.getTuples()) {
			inventoryLevels.put(tuple, sc.getInitInvLevels().get(tuple));
			inventoryPositions.put(tuple, sc.getInitInvLevels().get(tuple));
			TreeMap<Date,Integer> emptyMap = new TreeMap<Date,Integer>();
			for(Date d : sc.getTestDates()) {
				emptyMap.put(d, 0);
			}
			ordersToArriveMap.put(tuple, emptyMap);
		}
		
		// Keep track of sub-demands at each inventory
		Map<String,List<TreeMap<Date,Integer>>> subDemandMap = new HashMap<>();
		Map<String,TreeMap<Date,Integer>> demands = new HashMap<>();
		for(Inventory i : sc.getSortedInventories()) {
			for(Material m : i.getMaterials()) {
				String tuple = "" + i.getInventoryId() + "," + m.getMaterialCode();
				if(i.getStageType().equals("DIST")) {
					List<TreeMap<Date,Integer>> subDemandCopy = new ArrayList<TreeMap<Date,Integer>>(sc.getScSubDemands().get(tuple));
					subDemandMap.put(tuple, subDemandCopy);
					TreeMap<Date,Integer> demandsCopy = new TreeMap<Date,Integer>(sc.getDemandsMap().get(tuple));
					demands.put(tuple, demandsCopy);
				} else {
					List<TreeMap<Date,Integer>> emptyList = new ArrayList<TreeMap<Date,Integer>>();
					for(int j=0; j<i.getSuccessors().size(); j++) {
						TreeMap<Date,Integer> emptyMap = new TreeMap<Date,Integer>();
						for(Date d : sc.getDates()) {
							emptyMap.put(d, 0);
						}
						emptyList.add(emptyMap);
					}
					subDemandMap.put(tuple, emptyList);
					TreeMap<Date,Integer> emptyMap = new TreeMap<Date,Integer>();
					for(Date d : sc.getDates()) {
						emptyMap.put(d, 0);
					}
					demands.put(tuple, emptyMap);
				}
			}
		}
		
		for(Date d : sc.getTestDates()) {
			// 1) Create demand along the whole chain
			for(Inventory inventory : sc.getSortedInventories()) { 
				for(Material material : inventory.getMaterials()) {
					String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();					
					// Update inventory level and inventory positions
					inventoryLevels.put(tuple, inventoryLevels.get(tuple) + ordersToArriveMap.get(tuple).get(d));
					// Add transport costs for the sub-orders that can be fulfilled
					for(int i=0; i<subDemandMap.get(tuple).size(); i++) {
						TreeMap<Date,Integer> subDemands = subDemandMap.get(tuple).get(i);
						String successor = inventory.getSuccessors().get(i);
						int subDemand = subDemands.get(d);
						if(subDemand > 0) {
							if(inventory.getStageType().equals("DIST")) {
								numTotalOrders.put(tuple, numTotalOrders.get(tuple) + 1);
							}	
							if(inventoryLevels.get(tuple) >= subDemand) {
								// Sub-order CAN be satisfied with current inventory
								inventoryLevels.put(tuple, inventoryLevels.get(tuple) - subDemand);
								inventoryPositions.put(tuple, inventoryPositions.get(tuple) - subDemand);
								// Add transportation costs
								if(!inventory.getStageType().equals("PLANT")) {
									transportationCostMap.put(tuple, transportationCostMap.get(tuple) + inventory.getTransportCosts().get(successor).get(material) * subDemand);
								}
							} else {
								if(inventory.getStageType().equals("DIST")) {
									// Add penalty costs and increase counter of not-satisfied orders
									penaltyCostMap.put(tuple, penaltyCostMap.get(tuple) + inventory.getShortageCosts().get(material));
									// Only customer orders matter for the number of unmet orders
									numUnmetOrders.put(tuple, numUnmetOrders.get(tuple) + 1);
								}
							}
						}
					}
					// Input 3: Lead time (in hours) --> only fuzzy if NOT a distribution center, else 0
					int leadTime = sc.getLeadTimesMap().get(tuple).get(d);
					// Input 4: Market price
					double marketPrice = 0.0;
					if(inventory.getStageType().equals("PLANT")) {
						marketPrice = sc.getMarketPricesMap().get(tuple).get(d);
					}					
					
					// Create EOQ model				
					EOQ eoq = new EOQ(sc, inventory, material, trainingDemands.get(tuple), trainingLeadTimes.get(tuple), 0.95, false);
					int s = eoq.getReorderLevel();
					int S = eoq.getOrderUpToLevel();
					
					// Retrieve optimal order quantity by applying FIS model
					String predecessorStr = inventory.getPredecessor();
					Inventory predecessor = this.getInventoryByCode(predecessorStr);
					int indexSuccessor = predecessor.getSuccessors().indexOf(inventory.getInventoryId());
					String predecessorTuple = predecessorStr + "," + material.getMaterialCode();
					
					if(inventory.getStageType().equals("DIST")) {	
						int orderQnt = 0;
						if(inventoryPositions.get(tuple) <= s) {
							orderQnt = S - inventoryPositions.get(tuple);	
							numTotalOrders.put(predecessorTuple, numTotalOrders.get(predecessorTuple) + 1);
							demands.get(predecessorTuple).put(d, demands.get(predecessorTuple).get(d) + orderQnt);
							if(inventoryLevels.get(predecessorTuple) + ordersToArriveMap.get(predecessorTuple).get(d) >= orderQnt) {
								// Keep track of the date when this order will arrive (important for correct computation of the inventory level)
								Date newDate = this.addDays(d, (int) leadTime/24);
								if(newDate.equals(d)) {
									inventoryLevels.put(tuple, inventoryLevels.get(tuple) + orderQnt);
									inventoryPositions.put(tuple, inventoryPositions.get(tuple) + orderQnt);
								} else {
									if(!newDate.after(sc.getDates().last())) {
										ordersToArriveMap.get(tuple).put(newDate, ordersToArriveMap.get(tuple).get(newDate) + orderQnt);
									}
									inventoryPositions.put(tuple, inventoryPositions.get(tuple) + orderQnt);
								}
								int currentQnt = subDemandMap.get(predecessorTuple).get(indexSuccessor).get(d);
								subDemandMap.get(predecessorTuple).get(indexSuccessor).put(d, currentQnt + orderQnt);
								orderCostMap.put(tuple, orderCostMap.get(tuple) + inventory.getOrderCosts().get(material));	
							} else {
								// Add penalty costs and increase counter of not-satisfied orders
								penaltyCostMap.put(predecessorTuple, penaltyCostMap.get(predecessorTuple) + predecessor.getShortageCosts().get(material));
								// Only customer orders matter for the number of unmet orders
								numUnmetOrders.put(predecessorTuple, numUnmetOrders.get(predecessorTuple) + 1);
							}
						}	
						this.orderQuantitiesFIS.get(tuple).put(d, orderQnt);
					}
					
					if(inventory.getStageType().equals("WH")) {
						int orderQnt = 0;
						if(inventoryPositions.get(tuple) <= s) {
							orderQnt = S - inventoryPositions.get(tuple);
							numTotalOrders.put(predecessorTuple, numTotalOrders.get(predecessorTuple) + 1);
							for(var entry : material.getIngredients().entrySet()) {
								Material matROH = entry.getKey();
								double amount;
								if(matROH.getMaterialCode().equals("CC-P01") || matROH.getMaterialCode().equals("CC-P02") || matROH.getMaterialCode().equals("CC-P03") || matROH.getMaterialCode().equals("CC-P04")) {
									amount = orderQnt * material.getIngredients().get(matROH);
								} else {
									amount = orderQnt * material.getIngredients().get(matROH) * material.getUnitSize();
								}
								int qntRawMaterial = (int) Math.ceil(amount);
								String rawPredecessorTuple = "" + predecessorStr + "," + matROH.getMaterialCode();
								// Add new order to sub demand of predecessor
								demands.get(rawPredecessorTuple).put(d, demands.get(rawPredecessorTuple).get(d) + qntRawMaterial);
							}
							if(this.isSufficientRawMaterial(inventoryLevels, predecessor, material, orderQnt, predecessorStr, d, demands)) {
								Date newDate = this.addDays(d, (int) leadTime / 24);
								if(newDate.equals(d)) {
									inventoryLevels.put(tuple, inventoryLevels.get(tuple) + orderQnt);
									inventoryPositions.put(tuple, inventoryPositions.get(tuple) + orderQnt);
								} else {
									if(!newDate.after(sc.getDates().last())) {
										ordersToArriveMap.get(tuple).put(newDate, ordersToArriveMap.get(tuple).get(newDate) + orderQnt);
									}	
									inventoryPositions.put(tuple, inventoryPositions.get(tuple) + orderQnt);
								}
								for(var entry : material.getIngredients().entrySet()) {
									Material matROH = entry.getKey();
									double amount;
									if(matROH.getMaterialCode().equals("CC-P01") || matROH.getMaterialCode().equals("CC-P02") || matROH.getMaterialCode().equals("CC-P03") || matROH.getMaterialCode().equals("CC-P04")) {
										amount = orderQnt * material.getIngredients().get(matROH);
									} else {
										amount = orderQnt * material.getIngredients().get(matROH) * material.getUnitSize();
									}
									int qntRawMaterial = (int) Math.ceil(amount);
									String rawPredecessorTuple = "" + predecessorStr + "," + matROH.getMaterialCode();
									// Add new order to sub demand of predecessor
									int currentQnt = subDemandMap.get(rawPredecessorTuple).get(indexSuccessor).get(d);
									subDemandMap.get(rawPredecessorTuple).get(indexSuccessor).put(d, currentQnt + qntRawMaterial);
								}
								setupCostMap.put(tuple, setupCostMap.get(tuple) + inventory.getSetupCosts().get(material));
								productionCostMap.put(tuple, productionCostMap.get(tuple) + orderQnt * inventory.getProdUnitCosts().get(material));			
							} else {
								// Add penalty costs and increase counter of not-satisfied orders
								penaltyCostMap.put(predecessorTuple, penaltyCostMap.get(predecessorTuple) + predecessor.getShortageCosts().get(material));
								// Only customer orders matter for the number of unmet orders
								numUnmetOrders.put(predecessorTuple, numUnmetOrders.get(predecessorTuple) + 1);
							}
						}
						this.orderQuantitiesFIS.get(tuple).put(d, orderQnt);
					} 		
					if(inventory.getStageType().equals("PLANT")) {
						int orderQnt = 0;
						if(inventoryPositions.get(tuple) <= s) {
							orderQnt = S - inventoryPositions.get(tuple);
							orderCostMap.put(tuple, orderCostMap.get(tuple) + inventory.getOrderCosts().get(material));
							unitCostMap.put(tuple, unitCostMap.get(tuple) + orderQnt * marketPrice);	
							Date newDate = this.addDays(d, (int) leadTime / 24);
							if(newDate.equals(d)) {
								inventoryLevels.put(tuple, inventoryLevels.get(tuple) + orderQnt);
								inventoryPositions.put(tuple, inventoryPositions.get(tuple) + orderQnt);
							} else {
								if(!newDate.after(sc.getDates().last())) {
									ordersToArriveMap.get(tuple).put(newDate, ordersToArriveMap.get(tuple).get(newDate) + orderQnt);
								}
								inventoryPositions.put(tuple, inventoryPositions.get(tuple) + orderQnt);
							}
						}
						this.orderQuantitiesFIS.get(tuple).put(d, orderQnt);
					}
					// Add holding costs
					holdingCostMap.put(tuple, holdingCostMap.get(tuple) + inventoryLevels.get(tuple) * inventory.getHoldingCosts().get(material));	
				}
			}
		}
		
		// Update total costs for this date
		for(Inventory inventory : sc.getSortedInventories()) {
			for(Material material : inventory.getMaterials()) {
				String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
				double holdingCosts = holdingCostMap.get(tuple);
				double penaltyCosts = penaltyCostMap.get(tuple);
				double orderCosts = orderCostMap.get(tuple);
				double setupCosts = setupCostMap.get(tuple);
				double productionCosts = productionCostMap.get(tuple);
				double unitCosts = unitCostMap.get(tuple);
				double transportCosts = transportationCostMap.get(tuple);
				double costsToAdd = holdingCosts + penaltyCosts + orderCosts + setupCosts + productionCosts + unitCosts + transportCosts;
				totalCosts.put(tuple, costsToAdd);
				sc.getMaxScCost().put(tuple, 5 * costsToAdd);
			}
		}
	}
	
	/** 
	 * SIMULATION OF A SINGLE INVENTORY BY APPLYING FRBS	
	 * This method simulates the inventory of a particular 'inventory'-'material', whereby orders are placed according to a fuzzy rule-based system (FRBS)
	 * @param inventory
	 * @param material
	 * @param dataBase Data base (DB) of FRBS
	 * @param ruleBase Rule base (RB) of FRBS
	 * @param training Indicates whether simulation is applied to training or test phase.
	 * @param multiEchelon Indicates whether method is applied within the multi-echelon FRBS heuristic (true) or solely to this single inventory and material (false).
	 * @throws Exception 
	 */
	public void simulateInventory(Inventory inventory, Material material, double[] dataBase, HashSet<String[]> ruleBase, boolean training, boolean multiEchelon) throws Exception {
		this.startTime = System.nanoTime();
		TreeMap<Date,Integer> orderQuantities = new TreeMap<Date,Integer>();
		String stageType = inventory.getStageType();
		String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
		numTotalOrders.put(tuple, 0);
		numUnmetOrders.put(tuple, 0);
		
		// Cost variables
		double totalHoldingCosts = 0.0;
		double totalPenaltyCosts = 0.0;
		double totalOrderCosts = 0.0;
		double totalSetupCosts = 0.0;
		double totalTransportCosts = 0.0;
		double totalProductionCosts = 0.0;
		double totalUnitCosts = 0.0;
		
		// Retrieve cost parameters
		double holdingCosts = inventory.getHoldingCosts().get(material);
		double penaltyCosts = 0.0;
		if(material.getTypeCode().equals("FERT")) {
			penaltyCosts = inventory.getShortageCosts().get(material);
		}
		// Order costs for plants and distribution centers; Setup costs for warehouses;
		double orderCosts = 0.0;
		double setupCosts = 0.0;
		if(stageType.equals("PLANT") || stageType.equals("DIST")) {
			orderCosts = inventory.getOrderCosts().get(material);
		} else {
			setupCosts = inventory.getSetupCosts().get(material);
		} 
		
		double productionCosts = 0.0;
		if(inventory.getStageType().equals("WH")) {
			productionCosts = inventory.getProdUnitCosts().get(material);
		}

		// Transport costs
		Map<String,Double> transportCosts = new HashMap<String,Double>();
		List<String> successors = inventory.getSuccessors();
		for(String s : successors) {
			HashMap<Material,Double> map = inventory.getTransportCosts().get(s);
			double cost;
			if(material.getTypeCode().equals("FERT")) {
				cost = map.get(material);
			} else {
				cost = 0.0;
			}
			transportCosts.put(s, cost);
		}
		
		int inventoryLevel;
		int inventoryPosition;
		TreeMap<Date,Integer> ordersToArrive;
		if(training) {
			inventoryLevel = sc.getInitInvLevels().get(tuple); // Keeps track of current inventory level
			inventoryPosition = inventoryLevel;
			// A map to keep track when quantities of previously set orders will arrive
			ordersToArrive = new TreeMap<Date,Integer>();
			// Initialize map ordersToArrive with 0's
			for(Date d : sc.getTrainingDates()) {
				ordersToArrive.put(d, 0);
			}
		} else {
			inventoryLevel = sc.getFinalInventoryLevels().get(tuple);
			inventoryPosition = sc.getFinalInventoryPositions().get(tuple);
			TreeMap<Date,Integer> ordersToArriveCopy = new TreeMap<Date,Integer>(sc.getOrdersToArriveInTestPhase().get(tuple));
			ordersToArrive = ordersToArriveCopy;
		}
		
		// Obtain time period of interest
		TreeSet<Date> datesSim;
		if(training) {
			datesSim = sc.getTrainingDates();
		} else {
			datesSim = sc.getTestDates();
		}
		
		List<TreeMap<Date,Integer>> subDemands;
		if(multiEchelon) {
			subDemands = subDemandsMaps.get(tuple);
		} else {
			List<TreeMap<Date,Integer>> subDemandsCopy = new ArrayList<TreeMap<Date,Integer>>(sc.getScSubDemands().get(tuple));
			subDemands = subDemandsCopy;
		}
		// Retrieve training demands, lead times and market prices
		TreeMap<Date,Integer> demands = sc.getDemandsMap().get(tuple);
		// Note: Distribution centers (DIST) do not have fuzzy lead times
		TreeMap<Date,Integer> leadTimes = sc.getLeadTimesMap().get(tuple);
		// Note: Only plants (PLANT) have fuzzy market prices
		TreeMap<Date,Double> marketPrices = new TreeMap<Date,Double>();
		if(inventory.getStageType().equals("PLANT")) {
			marketPrices = sc.getMarketPricesMap().get(tuple);
		}
		
		// Retrieve total demand over the whole time horizon for this inventory and material
		int sumDemand = 0;
		for(Integer dm : demands.values()) {
			sumDemand += dm; 
		}
		totalDemands.put(tuple, sumDemand);
		// Simulate for every date d within the training period
		for(Date d : datesSim) {
			// Retrieve market prices for raw materials for this given day
			double unitCosts = 0.0;
			if(inventory.getStageType().equals("PLANT") && material.getTypeCode().equals("ROH")) {
				if(sc.isSimulation()) {
					unitCosts = marketPrices.get(d);
				} else {
					unitCosts = material.getPrices().get(d);
				}			
			}
			// Increase numTotalOrders by 1 if there is demand to a successor
			for(int i=0; i<subDemands.size(); i++) {
				int subDemand = subDemands.get(i).get(d);
				if(subDemand > 0) {
					numTotalOrders.put(tuple, numTotalOrders.get(tuple) + 1);
				}
			}

			// Input 1: Demand (total demand per day relevant for fuzzy model)
			int demand = demands.get(d);
			
			inventoryLevel += ordersToArrive.get(d);
			// Add transport costs for the suborders that can be fulfilled
			for(int i=0; i<subDemands.size(); i++) {
				String successor = inventory.getSuccessors().get(i);
				int subDemand = subDemands.get(i).get(d);
				double transportCost = transportCosts.get(successor);
				if(inventoryLevel >= subDemand) {	
					// Sub-order CAN be satisfied with current inventory
					// Add transportation costs
					totalTransportCosts += transportCost * subDemand;
					// Update inventory level for demand if order can be satisfied.
					inventoryLevel -= subDemand;
					inventoryPosition -= subDemand;
				} else {
					// Order CANNOT be satisfied immediately
					// Add penalty costs and increase counter of not-satisfied orders
					totalPenaltyCosts += penaltyCosts;
					// Only customer orders matter for the number of unmet orders
					numUnmetOrders.put(tuple, numUnmetOrders.get(tuple) + 1);
				}
			}
			
			// Input 3: Lead time (in hours) --> only fuzzy if NOT a distribution center, else 0
			int leadTime = this.getLeadTimeForeCastAvg(inventory, material, d);
			
			// Input 4: Market price --> only for plants (PLANT)
			double marketPrice = 0.0;
			if(inventory.getStageType().equals("PLANT")) {
				marketPrice = marketPrices.get(d);
			}
			
			// Retrieve optimal order quantity by applying FIS model
			String predecessorStr = inventory.getPredecessor();
			Inventory predecessor = this.getInventoryByCode(predecessorStr);
			int indexSuccessor = predecessor.getSuccessors().indexOf(inventory.getInventoryId());
			String predecessorTuple = predecessorStr + "," + material.getMaterialCode();
			// Retrieve optimal order quantity by applying FIS model
			int orderQnt;
			if(inventory.getStageType().equals("PLANT")) {
				// Input array for fuzzy model
				double[] inputs = {demand, inventoryPosition, leadTime, marketPrice};
				orderQnt = this.getOrderQuantity(inventory, material, inputs, sc.getVarNamesL()); // 3 labels applied
				orderQuantities.put(d, orderQnt);
			} else if(inventory.getStageType().equals("WH")) {
				double[] inputs = {demand, inventoryPosition, leadTime};
				orderQnt = this.getOrderQuantity(inventory, material, inputs, sc.getVarNamesM()); // 3 labels applied
				for(var entry : material.getIngredients().entrySet()) {
					Material matROH = entry.getKey();
					double amount;
					if(matROH.getMaterialCode().equals("CC-P01") || matROH.getMaterialCode().equals("CC-P02") || matROH.getMaterialCode().equals("CC-P03") || matROH.getMaterialCode().equals("CC-P04")) {
						amount = orderQnt * material.getIngredients().get(matROH);
					} else {
						amount = orderQnt * material.getIngredients().get(matROH) * material.getUnitSize();
					}
					int qntRawMaterial = (int) Math.ceil(amount);
					String rawPredecessorTuple = "" + predecessorStr + "," + matROH.getMaterialCode();
					// Add new order to sub demand of predecessor
					int currentQnt = subDemandsMaps.get(rawPredecessorTuple).get(indexSuccessor).get(d);
					subDemandsMaps.get(rawPredecessorTuple).get(indexSuccessor).put(d, currentQnt + qntRawMaterial);
				}
				orderQuantities.put(d, orderQnt);
			} else {				
				double[] inputs = {demand, inventoryPosition};
				orderQnt = this.getOrderQuantity(inventory, material, inputs, sc.getVarNamesS()); // 3 labels applied
				int currentQnt = subDemandsMaps.get(predecessorTuple).get(indexSuccessor).get(d);
				subDemandsMaps.get(predecessorTuple).get(indexSuccessor).put(d, currentQnt + orderQnt);
				orderQuantities.put(d, orderQnt);
			}
			
			if(orderQnt > 0) {
				// Keep track of the date when the new order will arrive
				Date dateToArrive;
				if(sc.isSimulation()) {
					dateToArrive = this.addDays(d, (int) leadTimes.get(d)/24);
				} else {
					dateToArrive = this.addWorkingDays(d, (int) leadTimes.get(d)/24);
				}
				if(dateToArrive.equals(d)) {
					inventoryLevel += orderQnt;
					inventoryPosition += orderQnt;
				} else {
					// Adapt inventory position
					inventoryPosition += orderQnt;
					// Add order to list ordersToArrive
					if(ordersToArrive.containsKey(dateToArrive)) {
						ordersToArrive.put(dateToArrive, ordersToArrive.get(dateToArrive) + orderQnt);
					} else {
						ordersToArrive.put(dateToArrive, orderQnt);
					}
				}				
				if(stageType.equals("PLANT") || stageType.equals("DIST")) {
					// Add order costs if an order is set
					totalOrderCosts += orderCosts;
				} else {
					// Add setup costs if an order is set
					totalSetupCosts += setupCosts;
				}
				
				// Add production unit costs if an order is set (from a WH)
				totalProductionCosts += orderQnt * productionCosts;
				
				// Add purchase costs of raw material
				totalUnitCosts += orderQnt * unitCosts;
			}		
			// Add holding costs 
			totalHoldingCosts += inventoryLevel * holdingCosts;
		}
		// Compute total SC costs for this inventory and material
		double totalCostsInventory = totalHoldingCosts + totalPenaltyCosts + totalOrderCosts + totalSetupCosts + totalTransportCosts + totalProductionCosts + totalUnitCosts;
		orderQuantitiesFIS.put(tuple, orderQuantities);
		// Add total costs to map
		totalCosts.put(tuple, totalCostsInventory);
		holdingCostMap.put(tuple, totalHoldingCosts);
		penaltyCostMap.put(tuple, totalPenaltyCosts);
		orderCostMap.put(tuple, totalOrderCosts);
		setupCostMap.put(tuple, totalSetupCosts);
		transportationCostMap.put(tuple, totalTransportCosts);
		productionCostMap.put(tuple, totalProductionCosts);
		unitCostMap.put(tuple, totalUnitCosts);
		runTime = System.nanoTime() - startTime;
	}
	
	/** Applies a FIS (Fuzzy Inference System) to compute an order quantity considering a set of input values ('inputs').
	 * @param inventory Inventory for which order is placed
	 * @param material Material for which order is placed
	 * @param inputs Array containing fuzzy inputs and fuzzy output
	 * @return a crisp order quantity
	 * @throws Exception
	 */
	public int getOrderQuantity(Inventory inventory, Material material, double[] inputs, List<String> varNames) throws Exception {
		String fileName = "kbGeneration" + inventory.getInventoryId() + "," + material.getMaterialCode() + ".fcl";
		FIS fis = FIS.load(fileName, true);

		if (fis == null) {
			System.err.println("Can't load file: '" + fileName + "'");
			System.exit(1);
		}

		// Get default function block
		FunctionBlock fb = fis.getFunctionBlock(null);

		// Set inputs
		for(int i=0; i<varNames.size()-1; i++) {
			fb.setVariable(varNames.get(i), inputs[i]);
		}

		// Evaluate
		fb.evaluate();

		// Show output variable's chart
		fb.getVariable(varNames.get(varNames.size()-1)).defuzzify();
		
		/*
		JFuzzyChart.get().chart(fb);
		Variable orderQuantity = fb.getVariable(varNames.get(varNames.size()-1));
		JFuzzyChart.get().chart(orderQuantity, orderQuantity.getDefuzzifier(), true);
		System.out.println("");
		System.out.print(inputs[0] + " | " + inputs[1] + " | " + inputs[2] + " | " + inputs[3]);
		System.out.println("");
		 */

		return (int) fb.getVariable(varNames.get(varNames.size()-1)).getValue();
	}
	
	/** Forecasts lead time by taking a 'steps'-moving average.
	 * @return an expected lead time
	 */	
	public int getLeadTimeForeCastMa(Inventory inventory, Material material, Date d, int steps) {
		// Retrieve orders over whole time horizon
		String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
		TreeMap<Date,Integer> ltTimes = sc.getLeadTimesMap().get(tuple);
		// Retrieve lead times that are relevant for the computation of the average
		Map<Date,Integer> leadTimes = new TreeMap<Date,Integer>();
		for (var entry : ltTimes.entrySet()) {
			if(entry.getKey().before(d)) {
				leadTimes.put(entry.getKey(),entry.getValue());
			}
		}
		Set<Date> ltDatesSet = leadTimes.keySet();
		// Sort dates in an increasing order 
		TreeSet<Date> ltDatesTreeSet = new TreeSet<Date>(ltDatesSet);
		ArrayList<Date> ltDates = new ArrayList<Date>(ltDatesTreeSet);
		double sum = 0.0;
		int denominator = 1; // Default value: is changed in the following
		// If less dates than 'step' available --> use only these
		if(ltDates.size() < steps) {
			for (Date date : ltDates) {
				sum += leadTimes.get(date);
			}
			denominator = ltDates.size();
		} else {
			// Retrieve latest 'steps' dates before d
			for(int i=ltDates.size()-1; i>ltDates.size()-1-steps; i--) {
				sum += leadTimes.get(ltDates.get(i));
			}
			denominator = steps;
		}
		int avg = (int) (Math.round((sum/denominator)/24) * 24);
		return avg;
	}
	
	/** Forecasts lead time by taking the sample average.
	 * @return an expected lead time
	 */	
	public int getLeadTimeForeCastAvg(Inventory inventory, Material material, Date d) {
		// Retrieve orders over whole time horizon
		String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
		TreeMap<Date,Integer> ltTimes = sc.getLeadTimesMap().get(tuple);
		// If not information about the lead time is available --> expect it to be 1 day
		int leadTimeForecast = 24;
		if(d.after(dates.first())) {
			int counter = 0;
			int sum = 0;
			for (var entry : ltTimes.entrySet()) {
			    if(entry.getKey().before(d)) {
			    	sum += entry.getValue();
			    	counter += 1;
			    }
			}
			leadTimeForecast = (int) (Math.round((sum/counter)/24) * 24);
		}
		return leadTimeForecast;
	}
	
    /**
     * Adds days to a date and returns new date. If the resulting date is on the weekend, add 2 days.
     * @param date Date to be incremented
     * @param days Number of days to be added
     * @return
     */
    public Date addWorkingDays(Date date, int days) {
        Calendar cal = Calendar.getInstance();
        cal.setTime(date);
        cal.add(Calendar.DATE, days);
        if (cal.get(Calendar.DAY_OF_WEEK) == Calendar.SATURDAY || cal.get(Calendar.DAY_OF_WEEK) == Calendar.SUNDAY) {
        	cal.add(Calendar.DATE, 2);
        }     
        return cal.getTime();
    }
    
	/**
     * Adds days to a date and returns new date.
     * @param date Date to be incremented
     * @param days Number of days to be added
     * @return
     */
    public Date addDays(Date date, int days) {
        Calendar cal = Calendar.getInstance();
        cal.setTime(date);
        cal.add(Calendar.DATE, days);
        return cal.getTime();
    }
    
    /** Returns fill rate of a particular inventory and material (end-product9.
     * @param inventory
     * @param material
     * @return
     */
    public double getFillRate(String tuple) {
    	int numUnmet = numUnmetOrders.get(tuple);
    	int totalNum = numTotalOrders.get(tuple);
    	double fillRate = 1 - ((double) numUnmet / totalNum);
    	return fillRate;
    }
    
    /** Returns fill rate of the whole SC based on the unsatisfied orders for end-customers.
     * @return SC fill rate
     */
    public double getSCFillRate() {
    	int numUnmet = 0;
    	int totalNum = 0;
    	for(Inventory i : inventories) {
    		for(Material m : i.getMaterials()) {
        		if(i.getStageType().equals("DIST")) {
        			String tuple = "" + i.getInventoryId() + "," + m.getMaterialCode();
        			numUnmet += numUnmetOrders.get(tuple);
        			totalNum += numTotalOrders.get(tuple);
        		}
    		}
    	}
    	double fillRate = 1 - ((double) numUnmet / totalNum);
    	return fillRate;
    }

	/**
	 * @return the sc
	 */
	public SC getSc() {
		return sc;
	}

	/**
	 * @return the inventories
	 */
	public List<Inventory> getInventories() {
		return inventories;
	}

	/**
	 * @return the dates
	 */
	public Set<Date> getDates() {
		return dates;
	}

	/**
	 * @return the splitDate
	 */
	public Date getSplitDate() {
		return splitDate;
	}
	
	/**
	 * @return the orderQuantitiesFIS
	 */
	public Map<String, TreeMap<Date, Integer>> getOrderQuantitiesFIS() {
		return orderQuantitiesFIS;
	}

	/**
	 * @return the numTotalOrders
	 */
	public Map<String, Integer> getNumTotalOrders() {
		return numTotalOrders;
	}

	/**
	 * @return the numUnmetOrders
	 */
	public Map<String, Integer> getNumUnmetOrders() {
		return numUnmetOrders;
	}
	
	/**
	 * @return the total costs of the entire supply chain
	 */
	public double getTotalScCosts() {
		double totalScCosts = 0.0;
		for (var entry : totalCosts.entrySet()) {
		    totalScCosts += entry.getValue();
		}
		return totalScCosts;
	}
	
	/**
	 * @return the total holding costs of the entire supply chain
	 */
	public double getTotalHoldingCosts() {
		double totalHoldingCosts = 0.0;
		for (var entry : holdingCostMap.entrySet()) {
			totalHoldingCosts += entry.getValue();
		}
		return totalHoldingCosts;
	}
	
	/**
	 * @return the total penalty costs of the entire supply chain
	 */
	public double getTotalPenaltyCosts() {
		double totalPenaltyCosts = 0.0;
		for (var entry : penaltyCostMap.entrySet()) {
			totalPenaltyCosts += entry.getValue();
		}
		return totalPenaltyCosts;
	}
	
	/**
	 * @return the total order costs of the entire supply chain
	 */
	public double getTotalOrderCosts() {
		double totalOrderCosts = 0.0;
		for (var entry : orderCostMap.entrySet()) {
			totalOrderCosts += entry.getValue();
		}
		return totalOrderCosts;
	}
	
	/**
	 * @return the total setup costs of the entire supply chain
	 */
	public double getTotalSetupCosts() {
		double totalSetupCosts = 0.0;
		for (var entry : setupCostMap.entrySet()) {
			totalSetupCosts += entry.getValue();
		}
		return totalSetupCosts;
	}
	
	/**
	 * @return the total production costs of the entire supply chain
	 */
	public double getTotalProductionCosts() {
		double totalProductionCosts = 0.0;
		for (var entry : productionCostMap.entrySet()) {
			totalProductionCosts += entry.getValue();
		}
		return totalProductionCosts;
	}
	
	/**
	 * @return the total unit costs of the entire supply chain
	 */
	public double getTotalUnitCosts() {
		double totalUnitCosts = 0.0;
		for (var entry : unitCostMap.entrySet()) {
			totalUnitCosts += entry.getValue();
		}
		return totalUnitCosts;
	}
	
	/**
	 * @return the total transportation costs of the entire supply chain
	 */
	public double getTotalTransportationCosts() {
		double totalTransportationCosts = 0.0;
		for (var entry : transportationCostMap.entrySet()) {
			totalTransportationCosts += entry.getValue();
		}
		return totalTransportationCosts;
	}
	
	/**
	 * @return the totalCosts
	 */
	public Map<String, Double> getTotalCosts() {
		return totalCosts;
	}
	
	/**
	 * @return the holdingCostMap
	 */
	public Map<String, Double> getHoldingCostMap() {
		return holdingCostMap;
	}

	/**
	 * @return the penaltyCostMap
	 */
	public Map<String, Double> getPenaltyCostMap() {
		return penaltyCostMap;
	}

	/**
	 * @return the orderCostMap
	 */
	public Map<String, Double> getOrderCostMap() {
		return orderCostMap;
	}

	/**
	 * @return the setupCostMap
	 */
	public Map<String, Double> getSetupCostMap() {
		return setupCostMap;
	}

	/**
	 * @return the productionCostMap
	 */
	public Map<String, Double> getProductionCostMap() {
		return productionCostMap;
	}

	/**
	 * @return the unitCostMap
	 */
	public Map<String, Double> getUnitCostMap() {
		return unitCostMap;
	}

	/**
	 * @return the transportationCostMap
	 */
	public Map<String, Double> getTransportationCostMap() {
		return transportationCostMap;
	}
	
	/**
	 * Returns supply chain costs per unit of a given 'endProduct'
	 * @param inventory
	 * @param material
	 * @return the total SC costs per product "material" demanded at inventory "inventory"
	 */
	public double getScCostPerUnitDemanded(Material endProduct) {
		int totalDemand = 0;
		double totalScCost = 0.0;
		double totalVolume = 0.0;
		
		// Compute the total volume of all customer demand (all end-products) during the whole time horizon
		for(Inventory i : inventories) {
			for(Material m : i.getMaterials()) {
				if(i.getStageType().equals("DIST")) {
					String tuple = "" + i.getInventoryId() + "," + m.getMaterialCode();
					double unitSize = m.getUnitSize();
					totalVolume += totalDemands.get(tuple) * unitSize;
				}
			}
		}
		
		// Compute the total volume of this material
		double totalVolumeProduct = 0.0;
		for(Inventory i : inventories) {
			String tuple = "" + i.getInventoryId() + "," + endProduct.getMaterialCode();
			if(i.getStageType().equals("DIST")) {
				double unitSize = endProduct.getUnitSize();
				totalDemand += totalDemands.get(tuple);
				totalVolumeProduct += totalDemands.get(tuple) * unitSize;
			}
			// SC costs corresponding to actual end-product
			if(!i.getStageType().equals("PLANT") && i.isActiveInv()) {
				totalScCost += totalCosts.get(tuple);
			}
		}
		
		// Compute total SC costs related to raw materials
		double totalPlantCost = 0.0;
		for(Inventory i : inventories) {
			if(i.getStageType().equals("PLANT")) {
				for(Material m : i.getMaterials()) {
					String tuple = "" + i.getInventoryId() + "," + m.getMaterialCode();
					totalPlantCost = totalCosts.get(tuple);
				}
			}
		}
		
		// Add plant costs according to the fraction of this end-product to all end-products that have been sold to a customer 
		totalScCost += (totalVolumeProduct/totalVolume) * totalPlantCost;
		// Return SC cost per product demanded
		return (double) totalScCost / totalDemand;	
	}
	
	/**
	 * Prints objective values and fill rates per inventory and material. 
	 */
	public void printResults() {
		for(Inventory i : inventories) {
			if(i.getStageType().equals("PLANT") || i.getStageType().equals("WH") || i.getStageType().equals("DIST")) {
				for(Material m : i.getMaterials()) {
					String tuple = "" + i.getInventoryId() + "," + m.getMaterialCode();
					System.out.println("Inventory/Material: " + tuple + " | " + "Objective value: " + Math.round(totalCosts.get(tuple) * 100.0) / 100.0 + " | " + "Fill rate: " + Math.round(this.getFillRate(tuple) * 100.0) / 100.0);
				}
			}
		}
		System.out.println("");
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
	 * @return the inventoryLevelsMap
	 */
	public Map<String, TreeMap<Date, Integer>> getInventoryLevelsMap() {
		return inventoryLevelsMap;
	}

	/**
	 * @return the inventoryPositionsMap
	 */
	public Map<String, TreeMap<Date, Integer>> getInventoryPositionsMap() {
		return inventoryPositionsMap;
	}
}