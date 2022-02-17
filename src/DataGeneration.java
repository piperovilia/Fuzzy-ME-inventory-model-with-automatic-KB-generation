import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.*;

import org.apache.commons.math3.distribution.GammaDistribution;

/**
 * This class generates the supply chain data of a German cereal producer that applies a EOQ-based local order policy.
 * 2 Options: 1) Data simulated (randomly sampled; Distributions used: Discrete, Gamma, Normal) or 2) Real data of company --> Which option retrieved by sc.isSimulation()
 */
public class DataGeneration {
	private SC sc;
	private DataReader reader;
	private List<Inventory> inventories;
	
	// Random generators
	private Random randDemandGenerator;
	private GammaDistribution randLtGeneratorWh;
	private GammaDistribution randLtGeneratorP;
	private Map<Material,Random> randomPrGenerators;
	private Random randPrGenerator;
	// Demand distribution data
	private double meanDemand;
	private double stdvDemand;
	private List<Integer> demandsList;
	private List<Double> cdfMapDemand;
	// Parameters Gamma distribution
	private double shapeWh;
	private double scaleWh;
	private double shapeP;
	private double scaleP;
	// Parameters Normal distribution
	private Map<Material,Double> meanPrices;
	private Map<Material,Double> priceStdvs;
	private TreeSet<Date> dates;
	private boolean continuousDemand;
	
	// Fuzzy input data maps
	private Map<String,TreeMap<Date,Integer>> demandsMap;
	private Map<String,TreeMap<Date,Integer>> demandChangesMap;
	private Map<String,TreeMap<Date,Integer>> inventoryLevelsMap;
	private Map<String,TreeMap<Date,Integer>> inventoryPositionsMap;
	private Map<String,TreeMap<Date,Integer>> leadTimesMap;
	private Map<String,TreeMap<Date,Double>> marketPricesMap;
	private Map<String,TreeMap<Date,Integer>> orderQuantities;
	private Map<String,List<TreeMap<Date,Integer>>> subDemandMap;
	private Map<String,Integer> totalDemandsTraining;
	private Map<String,Integer> totalDemandsTest;
	
	// Final results maps
	// TRAINING
	private Map<String,Integer> numTotalOrdersTraining;
	private Map<String,Integer> numUnmetOrdersTraining;
	private Map<String,Double> totalCostsTraining; 
	private Map<String,Double> holdingCostMapTraining;
	private Map<String,Double> penaltyCostMapTraining;
	private Map<String,Double> orderCostMapTraining;
	private Map<String,Double> setupCostMapTraining;
	private Map<String,Double> productionCostMapTraining;
	private Map<String,Double> unitCostMapTraining;
	private Map<String,Double> transportationCostMapTraining;
	//TEST
	private Map<String,Integer> numTotalOrdersTest;
	private Map<String,Integer> numUnmetOrdersTest;
	private Map<String,Double> totalCostsTest; 
	private Map<String,Double> holdingCostMapTest;
	private Map<String,Double> penaltyCostMapTest;
	private Map<String,Double> orderCostMapTest;
	private Map<String,Double> setupCostMapTest;
	private Map<String,Double> productionCostMapTest;
	private Map<String,Double> unitCostMapTest;
	private Map<String,Double> transportationCostMapTest;
	
	private Map<String,Integer> finalInventoryPositions; // Stores IPs on last training date
	private Map<String,Integer> finalInventoryLevels; // Stores ILs on last training date
	private Map<String,TreeMap<Date,Integer>> ordersToArriveInTestPhase; // Stores orders that are ordered in training phase, but arriving in test phase
	
	/**
	 * Constructor used for simulation.
	 * @param sc SC object
	 * @param reader
	 * @param dates All dates (training + test)
	 * @param demandsList List of possible discrete order quantities
	 * @param cdfMapDemand CDF table for demands
	 * @param shapeWh Shape parameter for generating lead times of orders set by the warehouse 
	 * @param scaleWh Scale parameter for generating lead times of orders set by the warehouse 
	 * @param shapeP Shape parameter for generating lead times of orders set by the plant 
	 * @param scaleP Scale parameter for generating lead times of orders set by the plant
	 * @param meanPrices Mean prices per raw material
	 * @param priceStdvs Standard deviations per raw material
	 */
	public DataGeneration(SC sc, DataReader reader, TreeSet<Date> dates, boolean continuousDemand, double meanDemand, double stdvDemand, List<Integer> demandsList, List<Double> cdfMapDemand, double shapeWh, double scaleWh, double shapeP, double scaleP, Map<Material,Double> meanPrices, Map<Material,Double> priceStdvs) {
		this.sc = sc;
		this.reader = reader;
		this.inventories = sc.getInventories();
		this.demandsList = demandsList;
		this.cdfMapDemand = cdfMapDemand;
		this.shapeWh = shapeWh;
		this.scaleWh = scaleWh;
		this.shapeP = shapeP;
		this.scaleP = scaleP;
		this.randDemandGenerator = new Random();
		this.randLtGeneratorWh = new GammaDistribution(shapeWh, scaleWh);
		this.randomPrGenerators = new HashMap<>();
		for(Material material : sc.getMaterials()) {
			if(material.getTypeCode().equals("ROH")) {
				Random random = new Random();
				randomPrGenerators.put(material, random);
			}
		}
		this.randLtGeneratorP = new GammaDistribution(shapeP, scaleP);
		this.randPrGenerator = new Random();
		this.dates = dates;
		this.continuousDemand = continuousDemand;
		this.meanDemand = meanDemand;
		this.stdvDemand = stdvDemand;
		this.meanPrices = meanPrices;
		this.priceStdvs = priceStdvs;
		this.demandsMap = new HashMap<>();
		this.demandChangesMap = new HashMap<>();
		this.inventoryLevelsMap = new HashMap<>();
		this.inventoryPositionsMap = new HashMap<>();
		this.leadTimesMap = new HashMap<>();
		this.marketPricesMap = new HashMap<>();
		this.orderQuantities = new HashMap<>();
		// Initialize with empty maps
		for(String tuple : sc.getTuples()) {
			TreeMap<Date,Integer> emptyMapIL = new TreeMap<Date,Integer>();
			TreeMap<Date,Integer> emptyMapIP = new TreeMap<Date,Integer>();
			TreeMap<Date,Integer> emptyMapLT = new TreeMap<Date,Integer>();
			TreeMap<Date,Double> emptyMapMP = new TreeMap<Date,Double>();
			TreeMap<Date,Integer> emptyMapOQ = new TreeMap<Date,Integer>();
			inventoryLevelsMap.put(tuple, emptyMapIL);
			inventoryPositionsMap.put(tuple, emptyMapIP);
			leadTimesMap.put(tuple, emptyMapLT);
			marketPricesMap.put(tuple, emptyMapMP);
			orderQuantities.put(tuple, emptyMapOQ);
		}
		this.subDemandMap = new HashMap<>();
		this.totalDemandsTraining = new HashMap<>();
		this.totalDemandsTest = new HashMap<>();
		
		this.numTotalOrdersTraining = new HashMap<>();
		this.numUnmetOrdersTraining = new HashMap<>();
		this.totalCostsTraining = new HashMap<>();
		this.holdingCostMapTraining = new HashMap<>();
		this.penaltyCostMapTraining = new HashMap<>();
		this.orderCostMapTraining = new HashMap<>();
		this.setupCostMapTraining = new HashMap<>();
		this.productionCostMapTraining = new HashMap<>();
		this.unitCostMapTraining = new HashMap<>();
		this.transportationCostMapTraining = new HashMap<>();
		
		this.numTotalOrdersTest = new HashMap<>();
		this.numUnmetOrdersTest = new HashMap<>();
		this.totalCostsTest = new HashMap<>();
		this.holdingCostMapTest = new HashMap<>();
		this.penaltyCostMapTest = new HashMap<>();
		this.orderCostMapTest = new HashMap<>();
		this.setupCostMapTest = new HashMap<>();
		this.productionCostMapTest = new HashMap<>();
		this.unitCostMapTest = new HashMap<>();
		this.transportationCostMapTest = new HashMap<>();
		
		// Add 0 to each result map for each tuple (inventory,material) 
		this.addZerosToMapInteger(numTotalOrdersTraining);
		this.addZerosToMapInteger(numUnmetOrdersTraining);
		this.addZerosToMapDouble(totalCostsTraining);
		this.addZerosToMapDouble(holdingCostMapTraining);
		this.addZerosToMapDouble(penaltyCostMapTraining);
		this.addZerosToMapDouble(orderCostMapTraining);
		this.addZerosToMapDouble(setupCostMapTraining);
		this.addZerosToMapDouble(productionCostMapTraining);
		this.addZerosToMapDouble(unitCostMapTraining);
		this.addZerosToMapDouble(transportationCostMapTraining);
		
		this.addZerosToMapInteger(numTotalOrdersTest);
		this.addZerosToMapInteger(numUnmetOrdersTest);
		this.addZerosToMapDouble(totalCostsTest);
		this.addZerosToMapDouble(holdingCostMapTest);
		this.addZerosToMapDouble(penaltyCostMapTest);
		this.addZerosToMapDouble(orderCostMapTest);
		this.addZerosToMapDouble(setupCostMapTest);
		this.addZerosToMapDouble(productionCostMapTest);
		this.addZerosToMapDouble(unitCostMapTest);
		this.addZerosToMapDouble(transportationCostMapTest);
		
		this.finalInventoryPositions = new HashMap<>();
		this.finalInventoryLevels = new HashMap<>();
		this.ordersToArriveInTestPhase = new HashMap<>();
		// Initialize ordersToArriveInTestPhase map
		for(Inventory i : inventories) {
			if(i.isActiveInv()) {
				for(Material m : i.getMaterials()) {
					String tuple = "" + i.getInventoryId() + "," + m.getMaterialCode();
					TreeMap<Date,Integer> emptyMap = new TreeMap<Date,Integer>();
					for(Date testDate : sc.getTestDates()) {
						emptyMap.put(testDate, 0);
					}
					ordersToArriveInTestPhase.put(tuple, emptyMap);
				}
			}
		}
	}
	
	/**
	 * Constructor used for real data set.
	 * @param sc
	 * @param reader
	 * @param dates
	 */
	public DataGeneration(SC sc, DataReader reader, TreeSet<Date> dates) {
		this.sc = sc;
		this.reader = reader;
		this.inventories = sc.getInventories();
		this.dates = dates;
		this.demandsMap = new HashMap<>();
		this.demandChangesMap = new HashMap<>();
		this.inventoryLevelsMap = new HashMap<>();
		this.inventoryPositionsMap = new HashMap<>();
		this.leadTimesMap = new HashMap<>();
		this.marketPricesMap = new HashMap<>();
		this.orderQuantities = new HashMap<>();
		// Initialize with empty maps
		for(String tuple : sc.getTuples()) {
			TreeMap<Date,Integer> emptyMapIL = new TreeMap<Date,Integer>();
			TreeMap<Date,Integer> emptyMapIP = new TreeMap<Date,Integer>();
			TreeMap<Date,Integer> emptyMapLT = new TreeMap<Date,Integer>();
			TreeMap<Date,Double> emptyMapMP = new TreeMap<Date,Double>();
			TreeMap<Date,Integer> emptyMapOQ = new TreeMap<Date,Integer>();
			inventoryLevelsMap.put(tuple, emptyMapIL);
			inventoryPositionsMap.put(tuple, emptyMapIP);
			leadTimesMap.put(tuple, emptyMapLT);
			marketPricesMap.put(tuple, emptyMapMP);
			orderQuantities.put(tuple, emptyMapOQ);
		}
		this.subDemandMap = new HashMap<>();
		this.totalDemandsTraining = new HashMap<>();
		this.totalDemandsTest = new HashMap<>();
		
		this.numTotalOrdersTraining = new HashMap<>();
		this.numUnmetOrdersTraining = new HashMap<>();
		this.totalCostsTraining = new HashMap<>();
		this.holdingCostMapTraining = new HashMap<>();
		this.penaltyCostMapTraining = new HashMap<>();
		this.orderCostMapTraining = new HashMap<>();
		this.setupCostMapTraining = new HashMap<>();
		this.productionCostMapTraining = new HashMap<>();
		this.unitCostMapTraining = new HashMap<>();
		this.transportationCostMapTraining = new HashMap<>();
		
		this.numTotalOrdersTest = new HashMap<>();
		this.numUnmetOrdersTest = new HashMap<>();
		this.totalCostsTest = new HashMap<>();
		this.holdingCostMapTest = new HashMap<>();
		this.penaltyCostMapTest = new HashMap<>();
		this.orderCostMapTest = new HashMap<>();
		this.setupCostMapTest = new HashMap<>();
		this.productionCostMapTest = new HashMap<>();
		this.unitCostMapTest = new HashMap<>();
		this.transportationCostMapTest = new HashMap<>();
		
		// Add 0's to result maps
		this.addZerosToMapInteger(numTotalOrdersTraining);
		this.addZerosToMapInteger(numUnmetOrdersTraining);
		this.addZerosToMapDouble(totalCostsTraining);
		this.addZerosToMapDouble(holdingCostMapTraining);
		this.addZerosToMapDouble(penaltyCostMapTraining);
		this.addZerosToMapDouble(orderCostMapTraining);
		this.addZerosToMapDouble(setupCostMapTraining);
		this.addZerosToMapDouble(productionCostMapTraining);
		this.addZerosToMapDouble(unitCostMapTraining);
		this.addZerosToMapDouble(transportationCostMapTraining);
		
		this.addZerosToMapInteger(numTotalOrdersTest);
		this.addZerosToMapInteger(numUnmetOrdersTest);
		this.addZerosToMapDouble(totalCostsTest);
		this.addZerosToMapDouble(holdingCostMapTest);
		this.addZerosToMapDouble(penaltyCostMapTest);
		this.addZerosToMapDouble(orderCostMapTest);
		this.addZerosToMapDouble(setupCostMapTest);
		this.addZerosToMapDouble(productionCostMapTest);
		this.addZerosToMapDouble(unitCostMapTest);
		this.addZerosToMapDouble(transportationCostMapTest);
		
		this.finalInventoryPositions = new HashMap<>();
		this.finalInventoryLevels = new HashMap<>();
		// Initialize ordersToArriveInTestPhase map
		this.ordersToArriveInTestPhase = new HashMap<>();
		for(Inventory i : inventories) {
			if(i.isActiveInv()) {
				for(Material m : i.getMaterials()) {
					String tuple = "" + i.getInventoryId() + "," + m.getMaterialCode();
					TreeMap<Date,Integer> emptyMap = new TreeMap<Date,Integer>();
					for(Date testDate : sc.getTestDates()) {
						emptyMap.put(testDate, 0);
					}
					ordersToArriveInTestPhase.put(tuple, emptyMap);
				}
			}
		}
	}
	
	/**
	 * Adds a 0.0 for each tuple (inventory, material) to a given map.
	 * @param mapToFill
	 */
	public void addZerosToMapDouble(Map<String,Double> mapToFill) {
		for(Inventory inventory : inventories) {
			for(Material material : sc.getMaterials()) {
				String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
				mapToFill.put(tuple, 0.0);
			}
		}
	}
	
	/**
	 * Adds a 0 for each tuple (inventory, material) to a given map.
	 * @param mapToFill
	 */
	public void addZerosToMapInteger(Map<String,Integer> mapToFill) {
		for(Inventory inventory : inventories) {
			for(Material material : sc.getMaterials()) {
				String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
				mapToFill.put(tuple, 0);
			}
		}
	}
	
	/**
	 * Returns if the current 'inventoryLevels' at a plant are sufficient to fulfill 'quantity' units of an end product 'endProduct'. 
	 * @param inventoryLevels Current inventory levels per tuple (inventory, material combination)
	 * @param inventory
	 * @param endProduct
	 * @param quantity
	 * @return
	 */
	public boolean isSufficientRawMaterial(Map<String,Integer> inventoryLevels, Inventory inventory, Material endProduct, int quantity) {
		boolean sufficient = true;
		// Go through all ingredients of this endProduct
		for(Material matROH : endProduct.getIngredients().keySet()) {
			String tuple = "" + inventory.getInventoryId() + "," + matROH.getMaterialCode();
			// Retrieve required raw material quantity for fulfilling the complete order
			double amount;
			if(matROH.getMaterialCode().equals("CC-P01") || matROH.getMaterialCode().equals("CC-P02") || matROH.getMaterialCode().equals("CC-P03") || matROH.getMaterialCode().equals("CC-P04")) {
				amount = quantity * endProduct.getIngredients().get(matROH);
			} else {
				amount = quantity * endProduct.getIngredients().get(matROH) * endProduct.getUnitSize();
			}	
			int qntRawMaterial = (int) Math.ceil(amount);
			// Check if enough of matROH available
			if(inventoryLevels.get(tuple) < qntRawMaterial) {
				sufficient = false;
			}
		}
		return sufficient;
	}
	

	/**
	 * Generates SC data for complete time horizon (training + test).
	 * @throws Exception
	 */
	public void generateSCData() throws Exception {
		Map<String,Integer> inventoryLevels = new HashMap<>();
		Map<String,Integer> inventoryPositions = new HashMap<>();
		Map<String,TreeMap<Date,Integer>> ordersToArriveMap = new HashMap<>();
		// Initialize IL and IP maps
		for(String tuple : sc.getTuples()) {
			inventoryLevels.put(tuple, sc.getInitInvLevels().get(tuple));
			inventoryPositions.put(tuple, sc.getInitInvLevels().get(tuple));
			TreeMap<Date,Integer> emptyMap = new TreeMap<Date,Integer>();
			for(Date d : sc.getDates()) {
				emptyMap.put(d, 0);
			}
			ordersToArriveMap.put(tuple, emptyMap);
		}
		// Initialize sub demands maps
		for(String tuple : sc.getTuples()) {
			List<TreeMap<Date,Integer>> emptyList = new ArrayList<TreeMap<Date,Integer>>();
			subDemandMap.put(tuple, emptyList);
		}
		// Initialize total demands maps
		Map<String,TreeMap<Date,Integer>> demands = new HashMap<>();
		for(String tuple : sc.getTuples()) {
			TreeMap<Date,Integer> emptyMap = new TreeMap<Date,Integer>();
			for(Date d : sc.getDates()) {
				emptyMap.put(d, 0);
			}
			demands.put(tuple, emptyMap);
		}
		
		// Create initial demands respectively sub-demands
		for(Inventory i : sc.getSortedInventories()) {
			for(Material m : i.getMaterials()) {
				String tuple = "" + i.getInventoryId() + "," + m.getMaterialCode();
				// If DC --> 1) sample demand from distribution (Simulation) OR 2) take known customer orders (Real application)
				if(i.getStageType().equals("DIST")) {
					for(int j=0; j<i.getSuccessors().size(); j++) {
						TreeMap<Date,Integer> customerDemand;
						if(sc.isSimulation()) {
							if(continuousDemand) {
								customerDemand = this.generateNormalDemand();
							} else {
								customerDemand = this.generateRandomDemand();
							}							
						} else {
							customerDemand = sc.getScSubDemands().get(tuple).get(j);
						}	
						subDemandMap.get(tuple).add(customerDemand);
						for(Date date : dates) {
							int oldQnt = demands.get(tuple).get(date);
							demands.get(tuple).put(date, oldQnt + customerDemand.get(date));
						}
					}
				} else {
					// If P or WH --> 0 demand
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
		// Store demands and lead times until this date --> for EOQ model
		Map<String,TreeMap<Date,Integer>> currentDemands = new HashMap<>();
		Map<String,TreeMap<Date,Integer>> currentLeadTimes = new HashMap<>();
		for(String tuple : sc.getTuples()) {
			TreeMap<Date,Integer> emptyMap1 = new TreeMap<Date,Integer>();
			TreeMap<Date,Integer> emptyMap2 = new TreeMap<Date,Integer>();
			currentDemands.put(tuple, emptyMap1);
			currentLeadTimes.put(tuple, emptyMap2);
		}
		
		// For each date (training + test) --> Create demand along the chain (from right to left: DC --> WH --> P)
		for(Date d : dates) {
			for(Inventory inventory : sc.getSortedInventories()) { 
				for(Material material : inventory.getMaterials()) {
					String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();					
					// Add orders that arrive today to IL
					inventoryLevels.put(tuple, inventoryLevels.get(tuple) + ordersToArriveMap.get(tuple).get(d));
					// Check if customer orders can be satisfied
					for(int i=0; i<subDemandMap.get(tuple).size(); i++) {
						TreeMap<Date,Integer> subDemands = subDemandMap.get(tuple).get(i);
						String successor = inventory.getSuccessors().get(i);
						int subDemand = subDemands.get(d);
						if(subDemand > 0) {
							if(inventory.getStageType().equals("DIST")) {
								// Increment total number of orders by 1
								if(!d.after(sc.getSplitDate())) {
									numTotalOrdersTraining.put(tuple, numTotalOrdersTraining.get(tuple) + 1);
								} else {
									numTotalOrdersTest.put(tuple, numTotalOrdersTest.get(tuple) + 1);
								} 
							}	
							if(inventoryLevels.get(tuple) >= subDemand) {
								// Sub-order CAN be satisfied with current inventory
								// Deduct order quantity from IL and IP
								inventoryLevels.put(tuple, inventoryLevels.get(tuple) - subDemand);
								inventoryPositions.put(tuple, inventoryPositions.get(tuple) - subDemand);
								// Add transportation costs
								if(!inventory.getStageType().equals("PLANT")) {
									if(!d.after(sc.getSplitDate())) {
										transportationCostMapTraining.put(tuple, transportationCostMapTraining.get(tuple) + inventory.getTransportCosts().get(successor).get(material) * subDemand);
									} else {
										transportationCostMapTest.put(tuple, transportationCostMapTest.get(tuple) + inventory.getTransportCosts().get(successor).get(material) * subDemand);
									}
								}
							} else {
								// Sub-order CANNOT be satisfied with current inventory
								if(inventory.getStageType().equals("DIST")) {
									// Add penalty costs and increase counter of unsatisfied orders
									if(!d.after(sc.getSplitDate())) {									
										penaltyCostMapTraining.put(tuple, penaltyCostMapTraining.get(tuple) + inventory.getShortageCosts().get(material));
										numUnmetOrdersTraining.put(tuple, numUnmetOrdersTraining.get(tuple) + 1);
									} else {
										penaltyCostMapTest.put(tuple, penaltyCostMapTest.get(tuple) + inventory.getShortageCosts().get(material));
										numUnmetOrdersTest.put(tuple, numUnmetOrdersTest.get(tuple) + 1);
									}
								}
							}
						}
					}
					// Retrieve DEMAND on this date
					int demand = demands.get(tuple).get(d);
					currentDemands.get(tuple).put(d, demand);
					// Generate LEAD TIME on this date
					double randomLeadTime = 0.0;
					// If DC --> LT is 24 hours (1 day)
					if(inventory.getStageType().equals("DIST")) {
						if(sc.isSimulation()) {
							randomLeadTime = 1;
						} else {
							randomLeadTime = 24;
						}
						
					}
					// If WH --> sample random LT or retrieve from given data
					if(inventory.getStageType().equals("WH")) {
						if(sc.isSimulation()) {
							randomLeadTime = randLtGeneratorWh.sample();
						} else {
							randomLeadTime = sc.getLeadTimesMap().get(tuple).get(d);
						}					
					}
					// If P --> sample random LT or retrieve from given data
					if(inventory.getStageType().equals("PLANT")) {
						if(sc.isSimulation()) {
							randomLeadTime = randLtGeneratorP.sample();
						} else {
							randomLeadTime = sc.getLeadTimesMap().get(tuple).get(d);
						}						
					}
					int leadTime;
					if(sc.isSimulation()) {
						int leadTimeDays = (int) Math.round(randomLeadTime);
						leadTime = leadTimeDays * 24;
					} else {
						leadTime = (int) randomLeadTime;
					}
					// Store newly created lead time
					this.leadTimesMap.get(tuple).put(d, leadTime);
					currentLeadTimes.get(tuple).put(d, leadTime);
					// Generate market price
					double marketPrice = 0.0;
					if(inventory.getStageType().equals("PLANT")) {
						// Sample or retrieve from given data
						if(sc.isSimulation()) {
							double mean = meanPrices.get(material);
							double stdv = priceStdvs.get(material);
							marketPrice = randomPrGenerators.get(material).nextGaussian() * stdv + mean;
						} else {
							marketPrice = reader.getMarketPrices().get(tuple).get(d);
						}
						// Store new market price
						this.marketPricesMap.get(tuple).put(d, marketPrice);
					}
					
					this.inventoryLevelsMap.get(tuple).put(d, inventoryLevels.get(tuple));
					this.inventoryPositionsMap.get(tuple).put(d, inventoryPositions.get(tuple));
					
					// Create EOQ model
					EOQ eoq = new EOQ(sc, inventory, material, currentDemands.get(tuple), currentLeadTimes.get(tuple), 0.95, false);
					int s = eoq.getReorderLevel();
					int S = eoq.getOrderUpToLevel();
					
					// Retrieve predecessor
					String predecessorStr = inventory.getPredecessor();
					Inventory predecessor = this.getInventoryByCode(predecessorStr);
					int indexSuccessor = predecessor.getSuccessors().indexOf(inventory.getInventoryId());
					String predecessorTuple = predecessorStr + "," + material.getMaterialCode();
					
					// IF DISTRIBUTION CENTER
					if(inventory.getStageType().equals("DIST")) {	
						int orderQnt = 0;
						// Place order if IP <= reorder level
						if(inventoryPositions.get(tuple) <= s) {
							// Order-up quantity 
							orderQnt = S - inventoryPositions.get(tuple);	
							// Increment total number of orders by 1
							if(!d.after(sc.getSplitDate())) {
								numTotalOrdersTraining.put(predecessorTuple, numTotalOrdersTraining.get(predecessorTuple) + 1);
							} else {
								numTotalOrdersTest.put(predecessorTuple, numTotalOrdersTest.get(predecessorTuple) + 1);
							}
							// Add order to total demand of predecessor
							demands.get(predecessorTuple).put(d, demands.get(predecessorTuple).get(d) + orderQnt);
							// Check if predecessor has enough stock
							if(inventoryLevels.get(predecessorTuple) + ordersToArriveMap.get(predecessorTuple).get(d) >= orderQnt) {
								// Predecessor has enough stock
								// Order is set --> retrieve arrival date
								Date newDate;
								if(sc.isSimulation()) {
									newDate = this.addDays(d, (int) leadTime/24);
								} else {
									newDate = this.addWorkingDays(d, (int) leadTime/24);
								}
								if(newDate.equals(d)) {
									// Order arrives on same date --> add order quantity directly to IL and IP
									inventoryLevels.put(tuple, inventoryLevels.get(tuple) + orderQnt);
									inventoryPositions.put(tuple, inventoryPositions.get(tuple) + orderQnt);
								} else {
									// Order arrives another day but before last test date --> increase IP
									if(!newDate.after(sc.getDates().last())) {
										ordersToArriveMap.get(tuple).put(newDate, ordersToArriveMap.get(tuple).get(newDate) + orderQnt);
									}
									// Increase IP
									inventoryPositions.put(tuple, inventoryPositions.get(tuple) + orderQnt);
								}
								// Keep track of products that are ordered in the training phase, but delivered in the test phase
								if(!d.after(sc.getSplitDate()) && newDate.after(sc.getSplitDate())) {
									this.ordersToArriveInTestPhase.get(tuple).put(newDate, orderQnt);
								}
								int currentQnt = subDemandMap.get(predecessorTuple).get(indexSuccessor).get(d);
								// // Add new order to sub demand of predecessor
								subDemandMap.get(predecessorTuple).get(indexSuccessor).put(d, currentQnt + orderQnt);
								// Add order costs
								if(!d.after(sc.getSplitDate())) {
									orderCostMapTraining.put(tuple, orderCostMapTraining.get(tuple) + inventory.getOrderCosts().get(material));
								} else {
									orderCostMapTest.put(tuple, orderCostMapTest.get(tuple) + inventory.getOrderCosts().get(material));
								}	
							} else {
								// Add penalty costs and increase counter of unsatisfied orders
								if(!d.after(sc.getSplitDate())) {								
									penaltyCostMapTraining.put(predecessorTuple, penaltyCostMapTraining.get(predecessorTuple) + predecessor.getShortageCosts().get(material));
									numUnmetOrdersTraining.put(predecessorTuple, numUnmetOrdersTraining.get(predecessorTuple) + 1);
								} else {
									penaltyCostMapTest.put(predecessorTuple, penaltyCostMapTest.get(predecessorTuple) + inventory.getShortageCosts().get(material));
									numUnmetOrdersTest.put(predecessorTuple, numUnmetOrdersTest.get(predecessorTuple) + 1);
								}
							}
						}
						// Store order quantity
						this.orderQuantities.get(tuple).put(d, orderQnt);
					}
					
					// IF WAREHOUSE
					if(inventory.getStageType().equals("WH")) {
						int orderQnt = 0;
						// Place order if IP <= reorder level
						if(inventoryPositions.get(tuple) <= s) {
							// Order-up quantity 
							orderQnt = S - inventoryPositions.get(tuple);
							// Increment total number of orders by 1
							if(!d.after(sc.getSplitDate())) {
								numTotalOrdersTraining.put(predecessorTuple, numTotalOrdersTraining.get(predecessorTuple) + 1);
							} else {
								numTotalOrdersTest.put(predecessorTuple, numTotalOrdersTest.get(predecessorTuple) + 1);
							}
							// Count raw material quantities as demand --> even if not sufficient stock --> important for input data set 
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
								// Add quantity of this raw material to predecessor
								demands.get(rawPredecessorTuple).put(d, demands.get(rawPredecessorTuple).get(d) + qntRawMaterial);
							}
							// Check if there is enough raw material for this order
							if(this.isSufficientRawMaterial(inventoryLevels, predecessor, material, orderQnt)) {
								// Enough raw material available
								// Order is set --> retrieve arrival date
								Date newDate;
								if(sc.isSimulation()) {
									newDate = this.addDays(d, (int) leadTime / 24);
								} else {
									newDate = this.addWorkingDays(d, (int) leadTime / 24);
								}
								if(newDate.equals(d)) {
									// Order arrives on same date --> add order quantity directly to IL and IP
									inventoryLevels.put(tuple, inventoryLevels.get(tuple) + orderQnt);
									inventoryPositions.put(tuple, inventoryPositions.get(tuple) + orderQnt);
								} else {
									// Order arrives another day but before last test date
									// Store order to arrive
									if(!newDate.after(sc.getDates().last())) {
										ordersToArriveMap.get(tuple).put(newDate, ordersToArriveMap.get(tuple).get(newDate) + orderQnt);
									}	
									// Increase IP
									inventoryPositions.put(tuple, inventoryPositions.get(tuple) + orderQnt);
								}
								// Keep track of products that are ordered in the training phase, but delivered in the test phase
								if(!d.after(sc.getSplitDate()) && newDate.after(sc.getSplitDate())) {
									this.ordersToArriveInTestPhase.get(tuple).put(newDate, orderQnt);
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
								// Add setup and production costs
								if(!d.after(sc.getSplitDate())) {
									setupCostMapTraining.put(tuple, setupCostMapTraining.get(tuple) + inventory.getSetupCosts().get(material));
									productionCostMapTraining.put(tuple, productionCostMapTraining.get(tuple) + orderQnt * inventory.getProdUnitCosts().get(material));
								} else {
									setupCostMapTest.put(tuple, setupCostMapTest.get(tuple) + inventory.getSetupCosts().get(material));
									productionCostMapTest.put(tuple, productionCostMapTest.get(tuple) + orderQnt * inventory.getProdUnitCosts().get(material));
								}			
							} else {
								// Add penalty costs and increase counter of unsatisfied orders
								if(!d.after(sc.getSplitDate())) {									
									penaltyCostMapTraining.put(predecessorTuple, penaltyCostMapTraining.get(predecessorTuple) + predecessor.getShortageCosts().get(material));
									numUnmetOrdersTraining.put(predecessorTuple, numUnmetOrdersTraining.get(predecessorTuple) + 1);
								} else {
									penaltyCostMapTest.put(predecessorTuple, penaltyCostMapTest.get(predecessorTuple) + inventory.getShortageCosts().get(material));
									numUnmetOrdersTest.put(predecessorTuple, numUnmetOrdersTest.get(predecessorTuple) + 1);
								}
							}
						}
						// Store order quantity
						this.orderQuantities.get(tuple).put(d, orderQnt);
					} 	
					
					// IF PLANT
					if(inventory.getStageType().equals("PLANT")) {
						int orderQnt = 0;
						// Place order if IP <= reorder level
						if(inventoryPositions.get(tuple) <= s) {	
							// Order-up quantity 
							orderQnt = S - inventoryPositions.get(tuple);
							// Add order and purchase costs
							if(!d.after(sc.getSplitDate())) {
								orderCostMapTraining.put(tuple, orderCostMapTraining.get(tuple) + inventory.getOrderCosts().get(material));
								unitCostMapTraining.put(tuple, unitCostMapTraining.get(tuple) + orderQnt * marketPrice);
							} else {
								orderCostMapTest.put(tuple, orderCostMapTest.get(tuple) + inventory.getOrderCosts().get(material));
								unitCostMapTest.put(tuple, unitCostMapTest.get(tuple) + orderQnt * marketPrice);
							}	
							// Order is set --> retrieve arrival date
							Date newDate;
							if(sc.isSimulation()) {
								newDate = this.addDays(d, (int) leadTime / 24);
							} else {
								newDate = this.addWorkingDays(d, (int) leadTime / 24);
							}
							
							if(newDate.equals(d)) {
								// Order arrives on same date --> add order quantity directly to IL and IP
								inventoryLevels.put(tuple, inventoryLevels.get(tuple) + orderQnt);
								inventoryPositions.put(tuple, inventoryPositions.get(tuple) + orderQnt);
							} else {
								// Order arrives another day but before last test date
								// Store order to arrive
								if(!newDate.after(sc.getDates().last())) {
									ordersToArriveMap.get(tuple).put(newDate, ordersToArriveMap.get(tuple).get(newDate) + orderQnt);
								}
								// Increase IP
								inventoryPositions.put(tuple, inventoryPositions.get(tuple) + orderQnt);
							}
							// Keep track of products that are ordered in the training phase, but delivered in the test phase
							if(!d.after(sc.getSplitDate()) && newDate.after(sc.getSplitDate())) {
								this.ordersToArriveInTestPhase.get(tuple).put(newDate, orderQnt);
							}
						}
						// Store order quantity 
						this.orderQuantities.get(tuple).put(d, orderQnt);
					}
					// Add holding costs
					if(!d.after(sc.getSplitDate())) {	
						holdingCostMapTraining.put(tuple, holdingCostMapTraining.get(tuple) + inventoryLevels.get(tuple) * inventory.getHoldingCosts().get(material));
					} else {
						holdingCostMapTest.put(tuple, holdingCostMapTest.get(tuple) + inventoryLevels.get(tuple) * inventory.getHoldingCosts().get(material));
					}	
					// Store final training IL and IP if last training date
					if(d.equals(sc.getSplitDate())) {
						finalInventoryPositions.put(tuple, inventoryPositions.get(tuple));
						finalInventoryLevels.put(tuple, inventoryLevels.get(tuple));
					}
				}
			}
		}
		// Store demand input for fuzzy model
		this.demandsMap = demands;
		
		// Update total costs for this date
		for(Inventory inventory : sc.getSortedInventories()) {
			for(Material material : inventory.getMaterials()) {
				this.generateDemandChanges(inventory, material);
				String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
				double holdingCostsTraining = holdingCostMapTraining.get(tuple);
				double penaltyCostsTraining = penaltyCostMapTraining.get(tuple);
				double orderCostsTraining = orderCostMapTraining.get(tuple);
				double setupCostsTraining = setupCostMapTraining.get(tuple);
				double productionCostsTraining = productionCostMapTraining.get(tuple);
				double unitCostsTraining = unitCostMapTraining.get(tuple);
				double transportCostsTraining = transportationCostMapTraining.get(tuple);
				double costsToAddTraining = holdingCostsTraining + penaltyCostsTraining + orderCostsTraining + setupCostsTraining + productionCostsTraining + unitCostsTraining + transportCostsTraining;
				totalCostsTraining.put(tuple, costsToAddTraining);
				// Max SC costs are 5 times total costs of EOQ
				sc.getMaxScCost().put(tuple, 5 * costsToAddTraining);
				
				double holdingCostsTest = holdingCostMapTest.get(tuple);
				double penaltyCostsTest = penaltyCostMapTest.get(tuple);
				double orderCostsTest = orderCostMapTest.get(tuple);
				double setupCostsTest = setupCostMapTest.get(tuple);
				double productionCostsTest = productionCostMapTest.get(tuple);
				double unitCostsTest = unitCostMapTest.get(tuple);
				double transportCostsTest = transportationCostMapTest.get(tuple);
				double costsToAddTest = holdingCostsTest + penaltyCostsTest + orderCostsTest + setupCostsTest + productionCostsTest + unitCostsTest + transportCostsTest;
				totalCostsTest.put(tuple, costsToAddTest);
			}
		}
		// Write .txt files with fuzzy input data
		this.writeData();
		// Add fuzzy input data to SC object
		this.addDataToSc();
		// Store info that is passed to test stage
		sc.setFinalInventoryPositions(finalInventoryPositions);
		sc.setFinalInventoryLevels(finalInventoryLevels);
		sc.setOrdersToArriveInTestPhase(ordersToArriveInTestPhase);
		// Add this data generator object to the SC
		sc.setDataGenerator(this);
	}
	
	/** 
	 * Generates demands of training and test phase. 
	 * Interarrival times of demands follow exponential distribution. 
	 * Demand quantities are sampled from given discrete probability distribution.
	 */
	public TreeMap<Date,Integer> generateRandomDemand() {
		Random randomExp = new Random();
		TreeMap<Date,Integer> demands = new TreeMap<Date,Integer>();
		for(Date d : dates) {
			demands.put(d, 0);
		}
		Date startDate = sc.getDates().first();
		Date currentDate = startDate;
		while(!currentDate.after(sc.getDates().last())) {
			double randomTime = Math.log(1-randomExp.nextDouble())/(-1);
			int randomTimeDays = (int) Math.round(randomTime);
			int randomDemand = 0;
			double randDouble = this.randDemandGenerator.nextDouble();
			for(int i=0; i<cdfMapDemand.size(); i++) {
			    if(i == 0) {
			    	if(randDouble <= cdfMapDemand.get(i)) {
			    		randomDemand =  demandsList.get(i);
			    	}
			    } else {
			    	if(randDouble <= cdfMapDemand.get(i) && randDouble > cdfMapDemand.get(i-1)) {
			    		randomDemand = demandsList.get(i);
			    	}
			    }
			}
			Date date;
			if(sc.isSimulation()) {
				date = this.addDays(currentDate, randomTimeDays);
			} else {
				date = this.addWorkingDays(currentDate, randomTimeDays);
			}
			if(!date.after(sc.getDates().last()) ) {
				demands.put(date, randomDemand);
			}
			currentDate = date;
		}
		return demands;
	}
	
	/** 
	 * Generates demands of training and test phase according to a normal distribution with a given mean and standard deviation.
	 * Interarrival times of demands follow exponential distribution. 
	 */
	public TreeMap<Date,Integer> generateNormalDemand() {
		Random randomExp = new Random();
		Random randomDemandGnr = new Random();
		TreeMap<Date,Integer> demands = new TreeMap<Date,Integer>();
		for(Date d : dates) {
			demands.put(d, 0);
		}
		Date startDate = sc.getDates().first();
		Date currentDate = startDate;
		while(!currentDate.after(sc.getDates().last())) {
			double randomTime = Math.log(1-randomExp.nextDouble())/(-1);
			int randomTimeDays = (int) Math.round(randomTime);
			double randomDemand = randomDemandGnr.nextGaussian() * stdvDemand + meanDemand;
			int demand = (int) Math.round(randomDemand);
			Date date;
			if(sc.isSimulation()) {
				date = this.addDays(currentDate, randomTimeDays);
			} else {
				date = this.addWorkingDays(currentDate, randomTimeDays);
			}
			if(!date.after(sc.getDates().last()) ) {
				demands.put(date, demand);
			}
			currentDate = date;
		}
		return demands;
	}
	
	/**
	 * Generates demand changes (current demand - last demand) corresponding to the generated demands.
	 */
	public void generateDemandChanges(Inventory inventory, Material material) {
		String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
		TreeMap<Date,Integer> demandChanges = new TreeMap<Date,Integer>();
		Map.Entry<Date,Integer> firstEntry = demandsMap.get(tuple).firstEntry();
		int lastDemand = firstEntry.getValue();
		for(Date d : dates) {
			if(d.equals(dates.first())) {
				demandChanges.put(d, 0);
			} else {
				int newDemand = demandsMap.get(tuple).get(d);
				int demandChange = newDemand - lastDemand;
				demandChanges.put(d, demandChange);
				lastDemand = newDemand;
			}
		}
		demandChangesMap.put(tuple, demandChanges);
	}
	
	/**
	 * Generates lead times of training and test phase according to a GAMMA(shape,scale) distribution.
	 * Mean = shape * scale
	 * Variance = shape * scale^2
	 * @param warehouse is 1 if warehouse, 0 if plant.
	 */
	public TreeMap<Date,Integer> generateLeadTimes(boolean warehouse) {
		TreeMap<Date,Integer> leadTimes = new TreeMap<Date,Integer>();
		for(Date d : dates) {
			double randomLeadTime;
			if(warehouse) {
				randomLeadTime = randLtGeneratorWh.sample();
			} else {
				randomLeadTime = randLtGeneratorP.sample();
			}
			int leadTime = (int) Math.round(randomLeadTime);
			int ltInHours = leadTime * 24;
			leadTimes.put(d, ltInHours);
		}
		return leadTimes;
	}
	
	/** 
	 * Generates market prices of training and test phase according to a normal distribution with a given mean and standard deviation.
	 * @param rawMaterial Material of interest
	 */
	public TreeMap<Date,Double> generateMarketPrices(Material rawMaterial) {
		TreeMap<Date,Double> marketPrices = new TreeMap<Date,Double>();
		double mean = meanPrices.get(rawMaterial);
		double stdv = priceStdvs.get(rawMaterial);
		for(Date d : dates) {
			double randomPrice = randPrGenerator.nextGaussian() * stdv + mean;
			if(randomPrice >= 0.0) {
				marketPrices.put(d, randomPrice);
			} else {
				marketPrices.put(d, 0.0);
			}		
		}
		return marketPrices;
	}
	
	/**
	 * @return the dates
	 */
	public TreeSet<Date> getDates() {
		return dates;
	}

	/**
	 * @return the demandsMap
	 */
	public Map<String, TreeMap<Date, Integer>> getDemandsMap() {
		return demandsMap;
	}

	/**
	 * @return the demandChangesMap
	 */
	public Map<String, TreeMap<Date, Integer>> getDemandChangesMap() {
		return demandChangesMap;
	}

	/**
	 * @return the inventoryLevelsMap
	 */
	public Map<String, TreeMap<Date, Integer>> getInventoryLevelsMap() {
		return inventoryLevelsMap;
	}

	/**
	 * @return the leadTimesMap
	 */
	public Map<String, TreeMap<Date, Integer>> getLeadTimesMap() {
		return leadTimesMap;
	}

	/**
	 * @return the marketPricesMap
	 */
	public Map<String, TreeMap<Date, Double>> getMarketPricesMap() {
		return marketPricesMap;
	}

	/**
	 * @return the orderQuantitiesEOQ
	 */
	public Map<String, TreeMap<Date, Integer>> getOrderQuantities() {
		return orderQuantities;
	}
	
    /**
	 * @return the inventoryPositionsMap
	 */
	public Map<String, TreeMap<Date, Integer>> getInventoryPositionsMap() {
		return inventoryPositionsMap;
	}

	/**
	 * @return the scSubDemands
	 */
	public Map<String, List<TreeMap<Date, Integer>>> getSubDemandMap() {
		return subDemandMap;
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
	
    /**
     * Writes a .txt file with all training input data for a given inventory and material.
     */
    public void writeInputData(Inventory inventory, Material material, TreeMap<Date,Integer> demands, TreeMap<Date,Integer> demandChanges,
    		TreeMap<Date,Integer> invPositions, TreeMap<Date,Integer> leadTimes, TreeMap<Date,Double> marketPrices, TreeMap<Date,Integer> orderQuantities) {
        try {
            FileWriter writer = new FileWriter("inputData" + inventory.getInventoryId() + "," + material.getMaterialCode() + ".txt", true);
            // Write field names
            if(inventory.getStageType().equals("PLANT")) {
            	writer.write("date demand demandChange inventoryPosition leadTime marketPrice orderQuantity");
            } else if(inventory.getStageType().equals("WH")) {
            	writer.write("date demand demandChange inventoryPosition leadTime orderQuantity");
            } else { 
            	writer.write("date demand demandChange inventoryPosition orderQuantity");
            }
            
            writer.write("\n");
            // Write input data
    		for(Date d : dates) {
    			String newDateStr = new SimpleDateFormat("yyyy-MM-dd").format(d);
    			writer.write(newDateStr + " ");
    			writer.write(demands.get(d) + " ");
    			writer.write(demandChanges.get(d) + " ");
    			writer.write(invPositions.get(d) + " ");
    			if(inventory.getStageType().equals("PLANT") || inventory.getStageType().equals("WH")) {
    				writer.write(leadTimes.get(d) + " ");
    			}
    			if(inventory.getStageType().equals("PLANT")) {
    				writer.write(marketPrices.get(d) + " ");
    			}
    			writer.write(orderQuantities.get(d) + " ");
    			writer.write("\n");
    		}
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
	
	/**
	 * Writes a .txt file with all training input data for all inventories and materials.
	 */
	public void writeData() {
		for(Inventory inventory : inventories) {
			if(inventory.isActiveInv()) {
				for(Material material : inventory.getMaterials()) {
					String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
					TreeMap<Date,Integer> demands = demandsMap.get(tuple);
					TreeMap<Date,Integer> demandChanges = demandChangesMap.get(tuple);
					TreeMap<Date,Integer> invPositions = inventoryPositionsMap.get(tuple);
					TreeMap<Date,Integer> leadTimes = leadTimesMap.get(tuple);
					TreeMap<Date,Double> marketPrices = marketPricesMap.get(tuple);
					TreeMap<Date,Integer> orderQnts = orderQuantities.get(tuple);
					this.writeInputData(inventory, material, demands, demandChanges, invPositions, leadTimes, marketPrices, orderQnts);
				}
			}
		}
	}
	
	/**
	 * Adds data to SC object.
	 */
	public void addDataToSc() {
		// Add copies of the maps --> Input data maps of this class remain unchanged
		Map<String,TreeMap<Date,Integer>> initDemandsMap = new HashMap<String,TreeMap<Date,Integer>>(demandsMap);
		Map<String,TreeMap<Date,Integer>> initDemandChangesMap = new HashMap<String,TreeMap<Date,Integer>>(demandChangesMap);
		Map<String,TreeMap<Date,Integer>> initInventoryLevelsMap = new HashMap<String,TreeMap<Date,Integer>>(inventoryLevelsMap);
		Map<String,TreeMap<Date,Integer>> initInventoryPositionsMap = new HashMap<String,TreeMap<Date,Integer>>(inventoryPositionsMap);
		Map<String,TreeMap<Date,Integer>> initLeadTimesMap = new HashMap<String,TreeMap<Date,Integer>>(leadTimesMap);
		Map<String,TreeMap<Date,Double>> initMarketPricesMap = new HashMap<String,TreeMap<Date,Double>>(marketPricesMap);
		Map<String,TreeMap<Date,Integer>> initOrderQuantities = new HashMap<String,TreeMap<Date,Integer>>(orderQuantities);
		Map<String,List<TreeMap<Date,Integer>>> initSubDemandMap = new HashMap<String,List<TreeMap<Date,Integer>>>(subDemandMap);
		sc.setDemandsMap(initDemandsMap);
		sc.setDemandChangesMap(initDemandChangesMap);
		sc.setInventoryLevelsMap(initInventoryLevelsMap);
		sc.setInventoryPositionsMap(initInventoryPositionsMap);
		sc.setLeadTimesMap(initLeadTimesMap);
		sc.setMarketPricesMap(initMarketPricesMap);
		sc.setOrderQuantities(initOrderQuantities);
		sc.setScSubDemands(initSubDemandMap);
	}
	
    /** Returns training fill rate for a particular inventory-material combination (tuple).
     * @param tuple String that identifies the stock of a material at a particular inventory 
     * @return fillRate
     */
    public double getTrainingFillRate(String tuple) {
    	int numUnmet = numUnmetOrdersTraining.get(tuple);
    	int totalNum = numTotalOrdersTraining.get(tuple);
    	double fillRate = 1 - ((double) numUnmet / totalNum);
    	return fillRate;
    }
    
    /** Returns training fill rate of the whole SC based on the unmet orders for end-customers.
     * @return SC fill rate
     */
    public double getTrainingSCFillRate() {
    	int numUnmet = 0;
    	int totalNum = 0;
    	for(Inventory i : inventories) {
    		for(Material m : i.getMaterials()) {
        		if(i.getStageType().equals("DIST")) {
        			String tuple = "" + i.getInventoryId() + "," + m.getMaterialCode();
        			numUnmet += numUnmetOrdersTraining.get(tuple);
        			totalNum += numTotalOrdersTraining.get(tuple);
        		}
    		}
    	}
    	double fillRate = 1 - ((double) numUnmet / totalNum);
    	return fillRate;
    }
    
    /** Returns test fill rate after applying a cost simulation.
     * @param tuple String that identifies the stock of a material at a particular inventory
     * @return
     */
    public double getTestFillRate(String tuple) {
    	int numUnmet = numUnmetOrdersTest.get(tuple);
    	int totalNum = numTotalOrdersTest.get(tuple);
    	double fillRate = 1 - ((double) numUnmet / totalNum);
    	return fillRate;
    }
    
    /** Returns test fill rate of the whole SC based on the unmet orders for end-customers.
     * @return SC fill rate
     */
    public double getTestSCFillRate() {
    	int numUnmet = 0;
    	int totalNum = 0;
    	for(Inventory i : inventories) {
    		for(Material m : i.getMaterials()) {
        		if(i.getStageType().equals("DIST")) {
        			String tuple = "" + i.getInventoryId() + "," + m.getMaterialCode();
        			numUnmet += numUnmetOrdersTest.get(tuple);
        			totalNum += numTotalOrdersTest.get(tuple);
        		}
    		}
    	}
    	double fillRate = 1 - ((double) numUnmet / totalNum);
    	return fillRate;
    }
    
	/**
	 * @return the total training costs of the entire supply chain
	 */
	public double getTrainingTotalScCosts() {
		double totalScCosts = 0.0;
		for (var entry : totalCostsTraining.entrySet()) {
		    totalScCosts += entry.getValue();
		}
		return totalScCosts;
	}
	
	/**
	 * @return the total test costs of the entire supply chain
	 */
	public double getTestTotalScCosts() {
		double totalScCosts = 0.0;
		for (var entry : totalCostsTest.entrySet()) {
		    totalScCosts += entry.getValue();
		}
		return totalScCosts;
	}
	
	/**
	 * @param endProduct
	 * @return the total training SC costs per product "material" demanded at inventory "inventory"
	 */
	public double getTrainingScCostPerUnitDemanded(Material endProduct) {
		int totalDemand = 0;
		double totalScCost = 0.0;
		double totalVolume = 0.0;
		
		// Compute the total volume of all customer demand (all end-products) during the whole time horizon
		for(Inventory i : inventories) {
			for(Material m : i.getMaterials()) {
				if(i.getStageType().equals("DIST")) {
					String tuple = "" + i.getInventoryId() + "," + m.getMaterialCode();
					double unitSize = m.getUnitSize();
					totalVolume += totalDemandsTraining.get(tuple) * unitSize;
				}
			}
		}
		
		// Compute the total volume of this material
		double totalVolumeProduct = 0.0;
		for(Inventory i : inventories) {
			String tuple = "" + i.getInventoryId() + "," + endProduct.getMaterialCode();
			if(i.getStageType().equals("DIST")) {
				double unitSize = endProduct.getUnitSize();
				totalDemand += totalDemandsTraining.get(tuple);
				totalVolumeProduct += totalDemandsTraining.get(tuple) * unitSize;
			}
			// SC costs corresponding to actual end-product
			if(!i.getStageType().equals("PLANT") && i.isActiveInv()) {
				totalScCost += totalCostsTraining.get(tuple);
			}
		}
		
		// Compute total SC costs related to raw materials
		double totalPlantCost = 0.0;
		for(Inventory i : inventories) {
			if(i.getStageType().equals("PLANT")) {
				for(Material m : i.getMaterials()) {
					String tuple = "" + i.getInventoryId() + "," + m.getMaterialCode();
					totalPlantCost = totalCostsTraining.get(tuple);
				}
			}
		}
		
		// Add plant costs according to the fraction of this end-product to all end-products that have been sold to a customer 
		totalScCost += (totalVolumeProduct/totalVolume) * totalPlantCost;
		// Return SC cost per product demanded
		return (double) totalScCost / totalDemand;	
	}
	
	/**
	 * @param inventory
	 * @param material
	 * @return the total test SC costs per product "material" demanded at inventory "inventory"
	 */
	public double getTestScCostPerUnitDemanded(Material endProduct) {
		int totalDemand = 0;
		double totalScCost = 0.0;
		double totalVolume = 0.0;
		
		// Compute the total volume of all customer demand (all end-products) during the whole time horizon
		for(Inventory i : inventories) {
			for(Material m : i.getMaterials()) {
				if(i.getStageType().equals("DIST")) {
					String tuple = "" + i.getInventoryId() + "," + m.getMaterialCode();
					double unitSize = m.getUnitSize();
					totalVolume += totalDemandsTest.get(tuple) * unitSize;
				}
			}
		}
		
		// Compute the total volume of this material
		double totalVolumeProduct = 0.0;
		for(Inventory i : inventories) {
			String tuple = "" + i.getInventoryId() + "," + endProduct.getMaterialCode();
			if(i.getStageType().equals("DIST")) {
				double unitSize = endProduct.getUnitSize();
				totalDemand += totalDemandsTest.get(tuple);
				totalVolumeProduct += totalDemandsTest.get(tuple) * unitSize;
			}
			// SC costs corresponding to actual end-product
			if(!i.getStageType().equals("PLANT") && i.isActiveInv()) {
				totalScCost += totalCostsTest.get(tuple);
			}
		}
		
		// Compute total SC costs related to raw materials
		double totalPlantCost = 0.0;
		for(Inventory i : inventories) {
			if(i.getStageType().equals("PLANT")) {
				for(Material m : i.getMaterials()) {
					String tuple = "" + i.getInventoryId() + "," + m.getMaterialCode();
					totalPlantCost = totalCostsTest.get(tuple);
				}
			}
		}
		
		// Add plant costs according to the fraction of this end-product to all end-products that have been sold to a customer 
		totalScCost += (totalVolumeProduct/totalVolume) * totalPlantCost;
		// Return SC cost per product demanded
		return (double) totalScCost / totalDemand;	
	}
	
	/**
	 * @return the total training holding costs of the entire supply chain
	 */
	public double getTrainingTotalHoldingCosts() {
		double totalHoldingCosts = 0.0;
		for (var entry : holdingCostMapTraining.entrySet()) {
			totalHoldingCosts += entry.getValue();
		}
		return totalHoldingCosts;
	}
	
	/**
	 * @return the total test holding costs of the entire supply chain
	 */
	public double getTestTotalHoldingCosts() {
		double totalHoldingCosts = 0.0;
		for (var entry : holdingCostMapTest.entrySet()) {
			totalHoldingCosts += entry.getValue();
		}
		return totalHoldingCosts;
	}
	
	/**
	 * @return the total training penalty costs of the entire supply chain
	 */
	public double getTrainingTotalPenaltyCosts() {
		double totalPenaltyCosts = 0.0;
		for (var entry : penaltyCostMapTraining.entrySet()) {
			totalPenaltyCosts += entry.getValue();
		}
		return totalPenaltyCosts;
	}
	
	/**
	 * @return the total test penalty costs of the entire supply chain
	 */
	public double getTestTotalPenaltyCosts() {
		double totalPenaltyCosts = 0.0;
		for (var entry : penaltyCostMapTest.entrySet()) {
			totalPenaltyCosts += entry.getValue();
		}
		return totalPenaltyCosts;
	}
	
	/**
	 * @return the total training order costs of the entire supply chain
	 */
	public double getTrainingTotalOrderCosts() {
		double totalOrderCosts = 0.0;
		for (var entry : orderCostMapTraining.entrySet()) {
			totalOrderCosts += entry.getValue();
		}
		return totalOrderCosts;
	}
	
	/**
	 * @return the total test order costs of the entire supply chain
	 */
	public double getTestTotalOrderCosts() {
		double totalOrderCosts = 0.0;
		for (var entry : orderCostMapTest.entrySet()) {
			totalOrderCosts += entry.getValue();
		}
		return totalOrderCosts;
	}
	
	/**
	 * @return the total training setup costs of the entire supply chain
	 */
	public double getTrainingTotalSetupCosts() {
		double totalSetupCosts = 0.0;
		for (var entry : setupCostMapTraining.entrySet()) {
			totalSetupCosts += entry.getValue();
		}
		return totalSetupCosts;
	}
	
	/**
	 * @return the total test setup costs of the entire supply chain
	 */
	public double getTestTotalSetupCosts() {
		double totalSetupCosts = 0.0;
		for (var entry : setupCostMapTest.entrySet()) {
			totalSetupCosts += entry.getValue();
		}
		return totalSetupCosts;
	}
	
	/**
	 * @return the total training production costs of the entire supply chain
	 */
	public double getTrainingTotalProductionCosts() {
		double totalProductionCosts = 0.0;
		for (var entry : productionCostMapTraining.entrySet()) {
			totalProductionCosts += entry.getValue();
		}
		return totalProductionCosts;
	}
	
	/**
	 * @return the total test production costs of the entire supply chain
	 */
	public double getTestTotalProductionCosts() {
		double totalProductionCosts = 0.0;
		for (var entry : productionCostMapTest.entrySet()) {
			totalProductionCosts += entry.getValue();
		}
		return totalProductionCosts;
	}
	
	/**
	 * @return the total training unit costs of the entire supply chain
	 */
	public double getTrainingTotalUnitCosts() {
		double totalUnitCosts = 0.0;
		for (var entry : unitCostMapTraining.entrySet()) {
			totalUnitCosts += entry.getValue();
		}
		return totalUnitCosts;
	}
	
	/**
	 * @return the total test unit costs of the entire supply chain
	 */
	public double getTestTotalUnitCosts() {
		double totalUnitCosts = 0.0;
		for (var entry : unitCostMapTest.entrySet()) {
			totalUnitCosts += entry.getValue();
		}
		return totalUnitCosts;
	}
	
	/**
	 * @return the total training transportation costs of the entire supply chain
	 */
	public double getTrainingTotalTransportationCosts() {
		double totalTransportationCosts = 0.0;
		for (var entry : transportationCostMapTraining.entrySet()) {
			totalTransportationCosts += entry.getValue();
		}
		return totalTransportationCosts;
	}
	
	/**
	 * @return the total test transportation costs of the entire supply chain
	 */
	public double getTestTotalTransportationCosts() {
		double totalTransportationCosts = 0.0;
		for (var entry : transportationCostMapTest.entrySet()) {
			totalTransportationCosts += entry.getValue();
		}
		return totalTransportationCosts;
	}

	/**
	 * @return the finalInventoryPositions
	 */
	public Map<String, Integer> getFinalInventoryPositions() {
		return finalInventoryPositions;
	}

	/**
	 * @return the finalInventoryLevels
	 */
	public Map<String, Integer> getFinalInventoryLevels() {
		return finalInventoryLevels;
	}

	/**
	 * @return the ordersToArriveInTestPhase
	 */
	public Map<String, TreeMap<Date, Integer>> getOrdersToArriveInTestPhase() {
		return ordersToArriveInTestPhase;
	}
	
	/**
	 * @return the totalCostsTraining
	 */
	public Map<String, Double> getTotalCostsTraining() {
		return totalCostsTraining;
	}

	/**
	 * @return the holdingCostMapTraining
	 */
	public Map<String, Double> getHoldingCostMapTraining() {
		return holdingCostMapTraining;
	}

	/**
	 * @return the penaltyCostMapTraining
	 */
	public Map<String, Double> getPenaltyCostMapTraining() {
		return penaltyCostMapTraining;
	}

	/**
	 * @return the orderCostMapTraining
	 */
	public Map<String, Double> getOrderCostMapTraining() {
		return orderCostMapTraining;
	}

	/**
	 * @return the setupCostMapTraining
	 */
	public Map<String, Double> getSetupCostMapTraining() {
		return setupCostMapTraining;
	}

	/**
	 * @return the productionCostMapTraining
	 */
	public Map<String, Double> getProductionCostMapTraining() {
		return productionCostMapTraining;
	}

	/**
	 * @return the unitCostMapTraining
	 */
	public Map<String, Double> getUnitCostMapTraining() {
		return unitCostMapTraining;
	}

	/**
	 * @return the transportationCostMapTraining
	 */
	public Map<String, Double> getTransportationCostMapTraining() {
		return transportationCostMapTraining;
	}

	/**
	 * @return the totalCostsTest
	 */
	public Map<String, Double> getTotalCostsTest() {
		return totalCostsTest;
	}

	/**
	 * @return the holdingCostMapTest
	 */
	public Map<String, Double> getHoldingCostMapTest() {
		return holdingCostMapTest;
	}

	/**
	 * @return the penaltyCostMapTest
	 */
	public Map<String, Double> getPenaltyCostMapTest() {
		return penaltyCostMapTest;
	}

	/**
	 * @return the orderCostMapTest
	 */
	public Map<String, Double> getOrderCostMapTest() {
		return orderCostMapTest;
	}

	/**
	 * @return the setupCostMapTest
	 */
	public Map<String, Double> getSetupCostMapTest() {
		return setupCostMapTest;
	}

	/**
	 * @return the productionCostMapTest
	 */
	public Map<String, Double> getProductionCostMapTest() {
		return productionCostMapTest;
	}

	/**
	 * @return the unitCostMapTest
	 */
	public Map<String, Double> getUnitCostMapTest() {
		return unitCostMapTest;
	}

	/**
	 * @return the transportationCostMapTest
	 */
	public Map<String, Double> getTransportationCostMapTest() {
		return transportationCostMapTest;
	}

	/**
	 * Prints relevant results.
	 */
	public void printResults() {
		System.out.println("TRAINING: ");
		System.out.println("Total SC costs: " + this.getTrainingTotalScCosts());
		System.out.println("Fill rate: " + this.getTrainingSCFillRate());
		System.out.println("Holding costs: " + this.getTrainingTotalHoldingCosts());
		System.out.println("");
		System.out.println("TEST: ");
		System.out.println("Total SC costs: " + this.getTestTotalScCosts());
		System.out.println("Fill rate: " + this.getTestSCFillRate());
		System.out.println("Holding costs: " + this.getTestTotalHoldingCosts());
	}
	
    /** Returns the inventory object corresponding to an inventoryId.
     * @param inventoryId
     * @return inventory
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
}