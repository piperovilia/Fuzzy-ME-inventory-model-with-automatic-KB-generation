import java.io.File;
import java.io.FileNotFoundException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 * A class to read and store data.
 *
 */
public class DataReader {
	private List<Customer> customers;
	private List<Material> materials;
	private List<Inventory> inventories;
	private List<String> arcs;
	private List<Order> orders;
	private SC supplyChain;
	private TreeSet<Date> dates;
	private List<String> invIds; // List of inventory IDs
	private Map<String,String> stageTypes;
	private Map<String,String> descriptions;
	private Map<String,HashMap<Material,Double>> holdingCosts;
	private Map<String,HashMap<Material,Double>> shortageCosts;
	private Map<String,HashMap<Material,Double>> setupCosts;
	private Map<String,HashMap<Material,Double>> orderCosts;
	private Map<String,HashMap<Material,Double>> prodUnitCosts;
	private Map<String,HashMap<Material,Double>> transportCosts; // String is combination "from,to"
	private Map<String,String> predecessor; // Map value is successor of Map key
	private Map<String,List<String>> successors;
	private Map<String,List<TreeMap<Date,Integer>>> demands;
	private Map<String,TreeMap<Date,Integer>> leadTimes;
	private Map<String,ArrayList<Double>> cumProbsMap;
	private Map<String,ArrayList<Integer>> leadTimeValues;
	private Map<String,TreeMap<Date,Double>> marketPrices;
	
	private double splitRatio;
	private boolean isDemand;
	private Date splitDate;
	private Random random;
	
    /**
     * Default constructor.
     * @param splitRatio Fraction of dates within training stage --> [0,1]
     * @param isDemand is 1 if demands are used and 0 if demand changes are used
     * @throws FileNotFoundException
     * @throws ParseException
     */
    public DataReader(double splitRatio, boolean isDemand) throws FileNotFoundException, ParseException {
        this.customers = new ArrayList<Customer>();
        this.materials = new ArrayList<Material>();
        this.inventories = new ArrayList<Inventory>();
        this.arcs = new ArrayList<String>();
        this.orders = new ArrayList<Order>();
        
        // Auxiliary lists for creation of inventory objects
        this.dates = new TreeSet<Date>();
        this.invIds = new ArrayList<String>();
        this.stageTypes = new HashMap<String,String>();
        this.descriptions = new HashMap<String,String>();
        this.holdingCosts = new HashMap<String,HashMap<Material,Double>>();
        this.shortageCosts = new HashMap<String,HashMap<Material,Double>>();
        this.setupCosts = new HashMap<String,HashMap<Material,Double>>();
        this.orderCosts = new HashMap<String,HashMap<Material,Double>>();
        this.prodUnitCosts = new HashMap<String,HashMap<Material,Double>>();
        this.transportCosts = new HashMap<String,HashMap<Material,Double>>();
        this.predecessor = new HashMap<String,String>();
        this.successors = new HashMap<String,List<String>>();
        this.demands = new HashMap<>();
        this.leadTimes = new HashMap<>();
        this.cumProbsMap = new HashMap<>();
        this.leadTimeValues = new HashMap<>();
        this.marketPrices = new HashMap<>();
        
        this.splitRatio = splitRatio;
        this.isDemand = isDemand;
        this.random = new Random(0);
    }
    
    /**
     * Reads data for Simulation.
     * @param startDate Start date of simulation
     * @param numDates Number of dates in simulation
     * @throws ParseException
     * @throws FileNotFoundException
     */
    public void readDataSimulation(Date startDate, int numDates) throws ParseException, FileNotFoundException {     
        // Create dates
		dates.add(startDate);
		Date currentDate = startDate;
		for(int i=1; i<numDates; i++) {
			Date newDate = this.addDays(currentDate, 1);
			dates.add(newDate);
			currentDate = newDate;
		}
		
        // Retrieve split date
        List<Date> datesList = new ArrayList<Date>(dates);
        int index = (int) (dates.size() * splitRatio);
        this.splitDate = datesList.get(index);
    	
        // Read materials
    	File fileMaterials = new File("Materials_Sim.txt");
        try (Scanner s = new Scanner(fileMaterials)) {
        	s.nextLine();
        	while(s.hasNext()) {
        		String materialCode = s.next();
        		String typeCode = s.next();
        		String unitSizeStr = s.next();
        		double unitSize = Double.parseDouble(unitSizeStr);
        		String description = s.next();
        		Material material = new Material(materialCode, typeCode, unitSize, description, null);
        		materials.add(material);
        	}
        }
        
        // Create recipes for end-products
        File fileIngredients = new File("BOM_Sim.txt");;    
        try (Scanner s = new Scanner(fileIngredients)) {
        	s.nextLine(); // Skip line with field names
        	while(s.hasNext()) {
        		String materialCode = s.next();
        		Material material = this.getMaterialByCode(materialCode);    		
        		
        		// CC-R01
        		String amountR01Str = s.next(); // percentage
        		double amountR01 = Double.parseDouble(amountR01Str); // expressed in kg or pcs.
        		Material r01 = this.getMaterialByCode("CC-R01");
        		material.addIngredient(r01, amountR01);
        		
        		// CC-R02
        		String amountR02Str = s.next(); // percentage
        		double amountR02 = Double.parseDouble(amountR02Str); // expressed in kg or pcs.
        		Material r02 = this.getMaterialByCode("CC-R02");
        		material.addIngredient(r02, amountR02);
        	}   
        }
        
        // Read arcs
        File fileArcs = new File("Arcs.txt");     
        try (Scanner s = new Scanner(fileArcs)) {
        	s.nextLine(); // Skip line with field names
        	while(s.hasNext()) {
            	String from = s.next();
            	String to = s.next();
            	String arc = "" + from + "," + to;
            	arcs.add(arc);
        	}   
        }
        
        // Read inventory objects
        // Read inventory id's and stage types
        File fileNodes = new File("Nodes.txt"); 
        try (Scanner s = new Scanner(fileNodes)) {
        	s.nextLine();
        	s.nextLine();        	
        	while(s.hasNext()) {
        		String id = s.next();
        		invIds.add(id);
        		String stageType = s.next();
        		stageTypes.put(id, stageType);
        		String desc = s.next();
        		descriptions.put(id, desc);
        	}
        }
        
        // Read holding, shortage and order costs
        File fileCosts = new File("Costs_Sim.txt");
        // Add empty maps for all inventories
        for(int i=0; i<5; i++) {
            for(String inventory : invIds) {
            	HashMap<Material,Double> emptyMapH = new HashMap<Material,Double>();
            	HashMap<Material,Double> emptyMapB = new HashMap<Material,Double>();
            	HashMap<Material,Double> emptyMapS = new HashMap<Material,Double>();
            	HashMap<Material,Double> emptyMapO = new HashMap<Material,Double>();
            	HashMap<Material,Double> emptyMapP = new HashMap<Material,Double>();
            	if(i==0) {
            		holdingCosts.put(inventory, emptyMapH);
            	}
            	if(i==1) {
            		shortageCosts.put(inventory, emptyMapB);
            	}
            	if(i==2) {
            		setupCosts.put(inventory, emptyMapS);
            	}
            	if(i==3) {
            		orderCosts.put(inventory, emptyMapO);
            	}
            	if(i==4) {
            		prodUnitCosts.put(inventory, emptyMapP);
            	}
            }
        }

        try (Scanner s = new Scanner(fileCosts)) {
        	s.nextLine(); // Skip line with field names
        	while(s.hasNext()) {
            	String invStr = s.next();
            	String matCodeStr = s.next();
            	Material mat = this.getMaterialByCode(matCodeStr);
            	String costType = s.next();
            	String costStr = s.next();
            	double cost = Double.parseDouble(costStr);
            	s.nextLine();
            	// Add corresponding cost to map
            	if(costType.equals("h")) {
            		holdingCosts.get(invStr).put(mat, cost);
            	}
            	if(costType.equals("b")) {
            		shortageCosts.get(invStr).put(mat, cost);
            	}
            	if(costType.equals("c_s")) {
            		setupCosts.get(invStr).put(mat, cost);
            	}
            	if(costType.equals("c_o")) {
            		orderCosts.get(invStr).put(mat, cost);
            	}
            	if(costType.equals("c_p")) {
            		prodUnitCosts.get(invStr).put(mat, cost);
            	}
        	}   
        }
        
        // Read transportation costs
        File fileTransportCosts = new File("Transport_Costs_Sim.txt");
        for(String arc : arcs) {
        	HashMap<Material,Double> emptyMap = new HashMap<Material,Double>();
        	transportCosts.put(arc, emptyMap);
        }
        try (Scanner s = new Scanner(fileTransportCosts)) {
        	s.nextLine(); // Skip line with field names
        	while(s.hasNext()) {
            	String from = s.next();
            	String to = s.next();
            	String str = "" + from + "," + to;
            	String matCodeStr = s.next();
            	Material mat = this.getMaterialByCode(matCodeStr);
            	String transportCostStr = s.next();
            	double transportCost = Double.parseDouble(transportCostStr);
            	transportCosts.get(str).put(mat, transportCost);
            	s.nextLine();
        	}   
        }
        
        // Read successor and predecessors of inventories
        // Add empty maps with predecessors
        for(String invId : this.invIds) {
        	ArrayList<String> emptyList = new ArrayList<String>();
        	successors.put(invId, emptyList);
        }
        for(String arc : arcs) {
        	String from = this.getFromFromArc(arc);
        	String to = this.getToFromArc(arc);
        	predecessor.put(to, from);
        	successors.get(from).add(to);
        }
        
        // Create actual inventory objects
        for(String invId : invIds) {
        	HashMap<String,HashMap<Material,Double>> invTransportCosts = new HashMap<String,HashMap<Material,Double>>();
        	for(String key : transportCosts.keySet()) {
        		String from = this.getFromFromArc(key);
        		String to = this.getToFromArc(key);
        		if(invId.equals(from)) {
        			invTransportCosts.put(to, transportCosts.get(key));
        		}
        	}
        	// Create list of materials/products used at particular inventory
        	List<Material> relevantMaterials = new ArrayList<Material>();
        	if(stageTypes.get(invId).equals("PLANT")) {
            	for(Material material : materials) {
            		if(material.getTypeCode().equals("ROH")) {
            			relevantMaterials.add(material);
            		}
            	}
        	} else {
        		for(Material material : materials) {
        			if(material.getTypeCode().equals("FERT")) {
        				relevantMaterials.add(material);
        			}
        		}
        	}
        	
        	// Only inventories of plants, warehouses or DCs are controlled (= are active)
        	boolean isActiveInv = false;
        	if(stageTypes.get(invId).equals("PLANT") || stageTypes.get(invId).equals("WH") || stageTypes.get(invId).equals("DIST")) {
        		isActiveInv = true;
        	}
        	
        	// Create inventory object
        	Inventory newInventory = new Inventory(invId, stageTypes.get(invId), descriptions.get(invId), holdingCosts.get(invId),
        			shortageCosts.get(invId), invTransportCosts, setupCosts.get(invId), orderCosts.get(invId), prodUnitCosts.get(invId), predecessor.get(invId),
        			successors.get(invId), relevantMaterials, splitDate, isActiveInv, dates);
        	inventories.add(newInventory);
        }
        
        // Read initial inventory levels
        HashMap<String,Integer> initInvLevels = new HashMap<String,Integer>();
        File fileInitInvLevels = new File("Init_Inv_Levels_Sim.txt");
        try (Scanner s = new Scanner(fileInitInvLevels)) {
        	while(s.hasNext()) {
        		s.nextLine();
				String tuple = s.next();
				int initInvLevel = s.nextInt();
        		initInvLevels.put(tuple, initInvLevel);
        	}
        }
        
        // Create supply chain object
        supplyChain = new SC(inventories, orders, materials, customers, arcs, dates, splitDate, isDemand);
        supplyChain.setInitInvLevels(initInvLevels);
    }
    
    /**
     * Reads data for real-world application.
     * @throws ParseException
     * @throws FileNotFoundException
     */
    public void readData() throws ParseException, FileNotFoundException {     
        // Create dates
    	int numDates;
    	File fileDates = new File("Dates.txt");
        try (Scanner s = new Scanner(fileDates)) {
        	numDates = s.nextInt();
        	s.next();
        	while(s.hasNext()) {
				String dateStr = s.next();
				Date date = new SimpleDateFormat("yyyy-MM-dd").parse(dateStr);
        		dates.add(date);
        	}
        }
		
        // Retrieve split date
        List<Date> datesList = new ArrayList<Date>(dates);
        int index = (int) (dates.size() * splitRatio);
        this.splitDate = datesList.get(index);
    	
        // Read materials
    	File fileMaterials = new File("Materials.txt");
        try (Scanner s = new Scanner(fileMaterials)) {
        	s.nextLine();
        	while(s.hasNext()) {
        		String materialCode = s.next();
        		String typeCode = s.next();
        		String unitSizeStr = s.next();
        		double unitSize = Double.parseDouble(unitSizeStr);
        		String description = s.next();
        		TreeMap<Date,Double> prices = new TreeMap<Date,Double>();
        		// Retrieve prices of raw material on different dates
        		if(typeCode.equals("ROH")) {
        			File filePurchases = new File("Purchases.txt");;
        			try (Scanner s2 = new Scanner(filePurchases)) {
        				// Create Map of prices for the raw material of interest
        				while(s2.hasNext()) {
        					s2.nextLine();
        					String matCode = s2.next();
        					if(matCode.equals(materialCode)) {
        						String orderDate = s2.next();
        						Date priceDate = new SimpleDateFormat("yyyy-MM-dd").parse(orderDate);
        						String matPriceStr = s2.next();
        						double matPrice = Double.parseDouble(matPriceStr);
        						prices.put(priceDate, matPrice);
        					}
        				}
        			}
        			// Fill missing prices --> If no price information available, take last known price
        			Date startDate = prices.firstKey();
        			double lastPrice = prices.get(startDate);
        			for(Date date : dates) {
        				if(prices.containsKey(date)) {
        					lastPrice = prices.get(date);
        				} else {
        					prices.put(date, lastPrice);
        				}
        			}
        		} else {
        			prices = null;
        		}
        		Material material = new Material(materialCode, typeCode, unitSize, description, prices);
        		materials.add(material);
        	}
        }
        
        // Create recipes for end-products
        File fileIngredients = new File("BOM.txt");;    
        try (Scanner s = new Scanner(fileIngredients)) {
        	s.nextLine(); // Skip line with field names
        	while(s.hasNext()) {
        		String materialCode = s.next();
        		Material material = this.getMaterialByCode(materialCode);           	
        		
        		// CC-P01
        		String amountP01Str = s.next(); // percentage
        		double amountP01 = Double.parseDouble(amountP01Str); // expressed in kg or pcs.
        		Material p01 = this.getMaterialByCode("CC-P01");
        		material.addIngredient(p01, amountP01);

        		// CC-P02
        		String amountP02Str = s.next(); // percentage
        		double amountP02 = Double.parseDouble(amountP02Str); // expressed in kg or pcs.
        		Material p02 = this.getMaterialByCode("CC-P02");
        		material.addIngredient(p02, amountP02);

        		// CC-P03
        		String amountP03Str = s.next(); // percentage      		
        		double amountP03 = Double.parseDouble(amountP03Str); // expressed in kg or pcs.
        		Material p03 = this.getMaterialByCode("CC-P03");
        		material.addIngredient(p03, amountP03);
        		
        		// CC-P04
        		String amountP04Str = s.next(); // percentage
        		double amountP04 = Double.parseDouble(amountP04Str); // expressed in kg or pcs.
        		Material p04 = this.getMaterialByCode("CC-P04");
        		material.addIngredient(p04, amountP04);
        		
        		// CC-R01
        		String amountR01Str = s.next(); // percentage
        		double amountR01 = Double.parseDouble(amountR01Str); // expressed in kg or pcs.
        		Material r01 = this.getMaterialByCode("CC-R01");
        		material.addIngredient(r01, amountR01);
        		
        		// CC-R02
        		String amountR02Str = s.next(); // percentage
        		double amountR02 = Double.parseDouble(amountR02Str); // expressed in kg or pcs.
        		Material r02 = this.getMaterialByCode("CC-R02");
        		material.addIngredient(r02, amountR02);
        		
        		// CC-R03
        		String amountR03Str = s.next(); // percentage
        		double amountR03 = Double.parseDouble(amountR03Str); // expressed in kg or pcs.
        		Material r03 = this.getMaterialByCode("CC-R03");
        		material.addIngredient(r03, amountR03);
        		
        		// CC-R04
        		String amountR04Str = s.next(); // percentage
        		double amountR04 = Double.parseDouble(amountR04Str); // expressed in kg or pcs.
        		Material r04 = this.getMaterialByCode("CC-R04");
        		material.addIngredient(r04, amountR04);
        		
        		// CC-R05
        		String amountR05Str = s.next(); // percentage
        		double amountR05 = Double.parseDouble(amountR05Str); // expressed in kg or pcs.
        		Material r05 = this.getMaterialByCode("CC-R05");
        		material.addIngredient(r05, amountR05);
        		
        		// CC-R06
        		String amountR06Str = s.next(); // percentage
        		double amountR06 = Double.parseDouble(amountR06Str); // expressed in kg or pcs.
        		Material r06 = this.getMaterialByCode("CC-R06");
        		material.addIngredient(r06, amountR06);
        	}   
        }
        
        // Read arcs
        File fileArcs = new File("Arcs.txt");     
        try (Scanner s = new Scanner(fileArcs)) {
        	s.nextLine(); // Skip line with field names
        	while(s.hasNext()) {
            	String from = s.next();
            	String to = s.next();
            	String arc = "" + from + "," + to;
            	arcs.add(arc);
        	}   
        }
        
        // Read inventory objects
        // Read inventory id's and stage types
        File fileNodes = new File("Nodes.txt"); 
        try (Scanner s = new Scanner(fileNodes)) {
        	s.nextLine();
        	s.nextLine();        	
        	while(s.hasNext()) {
        		String id = s.next();
        		invIds.add(id);
        		String stageType = s.next();
        		stageTypes.put(id, stageType);
        		String desc = s.next();
        		descriptions.put(id, desc);
        	}
        }
        
        // Read holding, shortage and order costs
        File fileCosts = new File("Costs.txt");
        // Add empty maps for all inventories
        for(int i=0; i<5; i++) {
            for(String inventory : invIds) {
            	HashMap<Material,Double> emptyMapH = new HashMap<Material,Double>();
            	HashMap<Material,Double> emptyMapB = new HashMap<Material,Double>();
            	HashMap<Material,Double> emptyMapS = new HashMap<Material,Double>();
            	HashMap<Material,Double> emptyMapO = new HashMap<Material,Double>();
            	HashMap<Material,Double> emptyMapP = new HashMap<Material,Double>();
            	if(i==0) {
            		holdingCosts.put(inventory, emptyMapH);
            	}
            	if(i==1) {
            		shortageCosts.put(inventory, emptyMapB);
            	}
            	if(i==2) {
            		setupCosts.put(inventory, emptyMapS);
            	}
            	if(i==3) {
            		orderCosts.put(inventory, emptyMapO);
            	}
            	if(i==4) {
            		prodUnitCosts.put(inventory, emptyMapP);
            	}
            }
        }

        try (Scanner s = new Scanner(fileCosts)) {
        	s.nextLine(); // Skip line with field names
        	while(s.hasNext()) {
            	String invStr = s.next();
            	String matCodeStr = s.next();
            	Material mat = this.getMaterialByCode(matCodeStr);
            	String costType = s.next();
            	String costStr = s.next();
            	double cost = Double.parseDouble(costStr);
            	s.nextLine();
            	// Add corresponding cost to map
            	if(costType.equals("h")) {
            		holdingCosts.get(invStr).put(mat, cost);
            	}
            	if(costType.equals("b")) {
            		shortageCosts.get(invStr).put(mat, cost);
            	}
            	if(costType.equals("c_s")) {
            		setupCosts.get(invStr).put(mat, cost);
            	}
            	if(costType.equals("c_o")) {
            		orderCosts.get(invStr).put(mat, cost);
            	}
            	if(costType.equals("c_p")) {
            		prodUnitCosts.get(invStr).put(mat, cost);
            	}
        	}   
        }
        
        // Read transportation costs
        File fileTransportCosts = new File("Transport_Costs.txt");
        for(String arc : arcs) {
        	HashMap<Material,Double> emptyMap = new HashMap<Material,Double>();
        	transportCosts.put(arc, emptyMap);
        }
        try (Scanner s = new Scanner(fileTransportCosts)) {
        	s.nextLine(); // Skip line with field names
        	while(s.hasNext()) {
            	String from = s.next();
            	String to = s.next();
            	String str = "" + from + "," + to;
            	String matCodeStr = s.next();
            	Material mat = this.getMaterialByCode(matCodeStr);
            	String transportCostStr = s.next();
            	double transportCost = Double.parseDouble(transportCostStr);
            	transportCosts.get(str).put(mat, transportCost);
            	s.nextLine();
        	}   
        }
        
        // Read successor and predecessors of inventories
        // Add empty maps with predecessors
        for(String invId : this.invIds) {
        	ArrayList<String> emptyList = new ArrayList<String>();
        	successors.put(invId, emptyList);
        }
        for(String arc : arcs) {
        	String from = this.getFromFromArc(arc);
        	String to = this.getToFromArc(arc);
        	predecessor.put(to, from);
        	successors.get(from).add(to);
        }
        
        // Create actual inventory objects
        for(String invId : invIds) {
        	HashMap<String,HashMap<Material,Double>> invTransportCosts = new HashMap<String,HashMap<Material,Double>>();
        	for(String key : transportCosts.keySet()) {
        		String from = this.getFromFromArc(key);
        		String to = this.getToFromArc(key);
        		if(invId.equals(from)) {
        			invTransportCosts.put(to, transportCosts.get(key));
        		}
        	}
        	// Create list of materials/products used at particular inventory
        	List<Material> relevantMaterials = new ArrayList<Material>();
        	if(stageTypes.get(invId).equals("PLANT")) {
            	for(Material material : materials) {
            		if(material.getTypeCode().equals("ROH")) {
            			relevantMaterials.add(material);
            		}
            	}
        	} else {
        		for(Material material : materials) {
        			if(material.getTypeCode().equals("FERT")) {
        				relevantMaterials.add(material);
        			}
        		}
        	}
        	
        	// Only inventories of plants, warehouses or DCs are controlled (= are active)cfc
        	boolean isActiveInv = false;
        	if(stageTypes.get(invId).equals("PLANT") || stageTypes.get(invId).equals("WH") || stageTypes.get(invId).equals("DIST")) {
        		isActiveInv = true;
        	}
        	
        	// Create inventory object
        	Inventory newInventory = new Inventory(invId, stageTypes.get(invId), descriptions.get(invId), holdingCosts.get(invId),
        			shortageCosts.get(invId), invTransportCosts, setupCosts.get(invId), orderCosts.get(invId), prodUnitCosts.get(invId), predecessor.get(invId),
        			successors.get(invId), relevantMaterials, splitDate, isActiveInv, dates);
        	inventories.add(newInventory);
        }
        
        // Retrieve market prices maps
        for(Inventory i : inventories) {
        	if(i.isActiveInv()) {
        		for(Material m : i.getMaterials()) {
        			String tuple = "" + i.getInventoryId() + "," + m.getMaterialCode();
        			TreeMap<Date,Double> materialPrices = m.getPrices();
        			marketPrices.put(tuple, materialPrices);
        		}
        	}
        }
        
        // Initialize lead time maps with 0's
        for(Inventory i : inventories) {
    		for(Material m : i.getMaterials()) {
    			String tuple = "" + i.getInventoryId() + "," + m.getMaterialCode();
    			TreeMap<Date,Integer> emptyMap = new TreeMap<Date,Integer>();
    			for(Date date : dates) {
    				emptyMap.put(date, 0);
    			}
    			leadTimes.put(tuple, emptyMap);
    		}
        }
        
        // Read demands -> Extracted from orders table
        File fileDemands = new File("Orders.txt");
        // Add empty maps with 0's
        for(Inventory i : inventories) {
        	if(i.isActiveInv()) {
        		for(Material m : i.getMaterials()) {
        			String tuple = "" + i.getInventoryId() + "," + m.getMaterialCode();
        			// Related to demand
        			List<TreeMap<Date,Integer>> emptyList = new ArrayList<TreeMap<Date,Integer>>();
        			for(String successor : i.getSuccessors()) {
        				TreeMap<Date,Integer> emptyMap = new TreeMap<Date,Integer>();
        				for(Date date : dates) {
        					emptyMap.put(date, 0);
        				}
        				emptyList.add(emptyMap);
        			}
        			demands.put(tuple, emptyList);
        			// Related to lead times
        			TreeMap<Date,Integer> emptyMapLt = new TreeMap<Date,Integer>();
        			leadTimes.put(tuple, emptyMapLt);
        		}
        	}
        }
        try (Scanner s = new Scanner(fileDemands)) {
        	s.nextLine(); // Skip line with field names
            while(s.hasNext()) {
            	// Skip order_number, stage_type
            	s.next();
            	s.next();
            	String from = s.next();
            	String to = s.next();
            	Inventory inv = this.getInventoryByCode(to);
            	String materialCode = s.next();
        		int quantity = s.nextInt();
        		String startDateStr = s.next();
        		Date startDate = new SimpleDateFormat("yyyy-MM-dd").parse(startDateStr);
        		String endDateStr = s.next();
        		Date endDate = new SimpleDateFormat("yyyy-MM-dd").parse(endDateStr);
        		int leadTime = s.nextInt();
            	if(!to.equals("V1") && !to.equals("V2") && inv.getStageType().equals("DIST")) {
                	// Demands are orders in reversed order --> Order: from -> to; Demand: to -> from; 
            		String tuple = to + "," + materialCode;
        			if(!endDate.after(dates.last())) {
        				int indexSuccessor = inv.getSuccessors().indexOf(from);
                		demands.get(tuple).get(indexSuccessor).put(endDate, quantity);
        			}
            	}
            	String str = from + "," + materialCode;
            	leadTimes.get(str).put(startDate, leadTime);
        		s.nextLine();
            }
        }
        
        // Fill missing lead times --> according to distribution derived from past data
        this.createLeadTimeCdf();
        for (var entry : leadTimes.entrySet()) {
            for(Date date : dates) {
            	if(!entry.getValue().containsKey(date)) {
            		int randomLeadTime = this.getRandomLeadTime(entry.getKey());
            		entry.getValue().put(date, randomLeadTime);
            	}
            }
        }
        
        // Read initial inventory levels
        HashMap<String,Integer> initInvLevels = new HashMap<String,Integer>();
        File fileInitInvLevels = new File("Init_Inv_Levels.txt");
        try (Scanner s = new Scanner(fileInitInvLevels)) {
        	while(s.hasNext()) {
        		s.nextLine();
				String tuple = s.next();
				int initInvLevel = s.nextInt();
        		initInvLevels.put(tuple, initInvLevel);
        	}
        }
        
        // Create supply chain object
        supplyChain = new SC(inventories, orders, materials, customers, arcs, dates, splitDate, isDemand);
        supplyChain.setScSubDemands(demands);
        supplyChain.setLeadTimesMap(leadTimes);
        supplyChain.setMarketPricesMap(marketPrices);
        supplyChain.setInitInvLevels(initInvLevels);
    }
    
	/**
	 * Creates a CDF for the lead time based on past data.
	 */
	public void createLeadTimeCdf() {
		for (var entry : leadTimes.entrySet()) {
			// Map that captures how many times a particular lead time has occurred in training set (Key, Value) = (lead time in hours, number of occurrences)
			HashMap<Integer,Integer> counters = new HashMap<Integer,Integer>();
			for(Integer leadTime : entry.getValue().values()) {
				if(counters.containsKey(leadTime)) {
					// Increase counter by 1
					counters.put(leadTime, counters.get(leadTime) + 1);
				} else {
					// Set counter to 1 if observed for the first time
					counters.put(leadTime, 1);
				}
			}
			// Get total sum of counted lead times in counters
			int sumCounters = 0;
			for (var c : counters.entrySet()) {
				sumCounters += c.getValue();
			}
			// Define cumulative probability mass function (discrete)
			Set<Integer> leadTimesSet = counters.keySet();
			ArrayList<Integer> leadTimes = new ArrayList<Integer>(leadTimesSet);
			ArrayList<Double> cumProbs = new ArrayList<Double>();
			double currentSum = 0.0;
			for (int i=0; i<leadTimes.size(); i++) {
				currentSum += (double) counters.get(leadTimes.get(i)) / sumCounters;
				cumProbs.add(currentSum); 
			}
			leadTimeValues.put(entry.getKey(), leadTimes);
			cumProbsMap.put(entry.getKey(), cumProbs);
		}
	}
	
	/**
	 * Returns a random lead time sampled from a discrete distribution that is derived from past data.
	 * @param tuple This String identifies inventory-material combination
	 */
	public int getRandomLeadTime(String tuple) {
		ArrayList<Integer> leadTimes = leadTimeValues.get(tuple);
		ArrayList<Double> cumProbs = cumProbsMap.get(tuple);
		int randomLeadTime = -1;
		// Create random number generator
		double randDouble = random.nextDouble();
		for(int i=0; i<cumProbs.size(); i++) {
			if(i == 0) {
				if(randDouble <= cumProbs.get(i)) {
					randomLeadTime = leadTimes.get(i);
				}
			} else {
				if(randDouble <= cumProbs.get(i) && randDouble > cumProbs.get(i-1)) {
					randomLeadTime = leadTimes.get(i);
				}
			}
		}
		return randomLeadTime;
	}
	
	/** Returns the material object corresponding to a materialCode.
     * @param materialCode
     * @return
     */
    public Material getMaterialByCode(String materialCode) {
    	Material material = null;
    	for(Material m : materials) {
    		if(m.getMaterialCode().equals(materialCode)) {
    			material = m;
    			break;
    		}
    	}
    	return material;
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
    
    /** Returns the inventory id of the "from"-node of an arc.
     * @param arc String representation of both inventory id's being part of an arc
     * @return
     */
    public String getFromFromArc(String arc) {
    	String str = "";
    	int indexComma = arc.indexOf(",");
        for (int i=0; i<indexComma; i++) {
            str = str + arc.charAt(i);
        }
        return str;
    }
    
    /** Returns the inventory id of the "to"-node of an arc.
     * @param arc String representation of both inventory id's being part of an arc
     * @return
     */
    public String getToFromArc(String arc) {
    	String str = "";
    	int indexComma = arc.indexOf(",");
        for (int i=indexComma+1; i<arc.length(); i++) {
            str = str + arc.charAt(i);
        }
        return str;
    }
    
	/**
	 * @return the supplyChain
	 */
	public SC getSupplyChain() {
		return supplyChain;
	}

	/**
	 * @return the dates
	 */
	public TreeSet<Date> getDates() {
		return dates;
	}

	/**
	 * @return the demands
	 */
	public Map<String, List<TreeMap<Date, Integer>>> getDemands() {
		return demands;
	}

	/**
	 * @return the leadTimes
	 */
	public Map<String, TreeMap<Date, Integer>> getLeadTimes() {
		return leadTimes;
	}

	/**
	 * @return the marketPrices
	 */
	public Map<String, TreeMap<Date, Double>> getMarketPrices() {
		return marketPrices;
	}

	/**
	 * @return the splitRatio
	 */
	public double getSplitRatio() {
		return splitRatio;
	}

	/**
	 * @return the splitDate
	 */
	public Date getSplitDate() {
		return splitDate;
	}

	/**
	 * @return the customers
	 */
	public List<Customer> getCustomers() {
		return customers;
	}
	
	/**
	 * @return the materials
	 */
	public List<Material> getMaterials() {
		return materials;
	}
    
    /**
	 * @return the inventories
	 */
	public List<Inventory> getInventories() {
		return inventories;
	}

	/**
	 * @return the arcs
	 */
	public List<String> getArcs() {
		return arcs;
	}

	/**
	 * @return the orders
	 */
	public List<Order> getOrders() {
		return orders;
	}

	/**
	 * @return the invIds
	 */
	public List<String> getInvIds() {
		return invIds;
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
}