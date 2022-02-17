import java.io.File;
import java.io.FileNotFoundException;
import java.text.ParseException;
import java.util.*;
import java.util.concurrent.TimeUnit;

/**
 * WANG-MENDEL METHOD:
 * A class to automatically establish an initial fuzzy knowledge base based out of a given data set.
 * The universes of discourse of the fuzzy sets are defined by the MIN respectively MAX values observed in the past data.
 * Unique rules are defined by applying Wang-Mendel method.
 * The method is applied for a particular inventory and material.
 */
public class WMGeneration {
	private SC sc;
	private Inventory inventory;
	private Material material;
	private String tuple; // String identifying inventory and material
	private int numLabels; // Number of labels per fuzzy variable
	private double[] dataBase; // Stores definition points of membership functions
	private HashSet<String[]> ruleBase; // Stores unique Mamdani-Type (IF-THEN) fuzzy rules
	private KnowledgeBase wmKb; // The knowledge base after applying the WM method
	// Fuzzy input data sets
	private TreeMap<Date,Integer> demands;
	private TreeMap<Date,Integer> demandChanges;
	private TreeMap<Date,Integer> invPositions;
	private TreeMap<Date,Integer> leadTimes;
	private TreeMap<Date,Double> marketPrices;
	private TreeMap<Date,Integer> orderQnts;
	private boolean isDemand;
		
	/** 
	 * Default constructor.
	 * @param sc The supply chain
	 * @param inventory The inventory of interest
	 * @param material The material of interest
	 * @param numLabels The number of fuzzy labels representing a fuzzy variable
	 * @throws FileNotFoundException
	 */
	public WMGeneration(SC sc, Inventory inventory, Material material, int numLabels) throws FileNotFoundException {
		this.sc = sc;
		this.numLabels = numLabels;
		this.isDemand = sc.isDemand();
		
		// Initialize objects, variables, lists and maps
		this.inventory = inventory;
		this.material = material;
		this.tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
		this.dataBase = new double[inventory.getNumVars() * numLabels * 3];
		this.ruleBase = new HashSet<String[]>();
		this.wmKb = null;
		this.demands = new TreeMap<Date,Integer>();
		this.invPositions = new TreeMap<Date,Integer>();
		this.demandChanges = new TreeMap<Date,Integer>();
		this.leadTimes = new TreeMap<Date,Integer>();
		this.marketPrices = new TreeMap<Date,Double>();
		this.orderQnts = new TreeMap<Date,Integer>();
		
		// Extract input data from .txt input data files
		this.createInputDataFromFile();
	}
	
	/**
	 * Creates input data by reading data from .txt input data files.
	 * @throws FileNotFoundException
	 */
	public void createInputDataFromFile() throws FileNotFoundException {
        File fileInput = new File("inputData" + inventory.getInventoryId() + "," + material.getMaterialCode() + ".txt");
        try (Scanner s = new Scanner(fileInput)) {
        	s.nextLine(); // Skip line with field names
        	for(Date d : sc.getDates()) {
        		if(!d.after(sc.getSplitDate())) {
                	// Skip date
            		s.next();
            		
            		int demand = s.nextInt();
            		demands.put(d, demand);
            		
            		int demandChange = s.nextInt();
            		demandChanges.put(d, demandChange);
            		
            		int inventoryPosition = s.nextInt();
            		invPositions.put(d, inventoryPosition);
            		
            		// Distribution centers --> constant lead time (no fuzzy variable)
            		if(inventory.getStageType().equals("PLANT") || inventory.getStageType().equals("WH")) {
                		int leadTime = s.nextInt();
                		leadTimes.put(d, leadTime);
            		}
            		
            		// Market prices are only relevant for plants
            		if(inventory.getStageType().equals("PLANT")) {
            			String marketPriceStr = s.next();
            			double marketPrice = Double.parseDouble(marketPriceStr);
            			marketPrices.put(d, marketPrice);
            		}
            		
    				int orderQnt = s.nextInt();
    				orderQnts.put(d, orderQnt);
        		}
        	}
        }
	}
		
	/**
	 * Creates a data base by defining the definitions points of all membership functions.
	 */
	public void createDataBase() {
		// List with the universes of discourses (ranges) per fuzzy variable
		List<double[]> ranges = this.getRanges(isDemand);

		// Create array containing definition points
		// One variable consists of "numLabels" (odd number) labels. Every label (represented by a MF) is explained by 3 definition points.
		// All definition points (of all variables and labels) are stored in one array
		for(int i=0; i<inventory.getNumVars(); i++) {
			// Retrieve relevant universe of discourse (range)
			double[] range = ranges.get(i);
			double distRange = range[1] - range[0];
			double shift = distRange / 2;
			double min = range[0] - shift; // Left border of universe of discourse
			double max = range[1] + shift; // Right border of universe of discourse
			double univDiscourse = max - min; // Whole range
			int numParts = numLabels + 1; // The centers (middle definition points) of MFs are defined by splitting the universe of discourse into (number of labels + 1) equally sized parts
			// Retrieve distance to the next definition point (equal in this case)
			double distToNextPos = univDiscourse / numParts;
			// Start on the far left of the universe of discourse
			double position = min;
			for(int j=0; j<numLabels; j++) {
				for(int k=0; k<3; k++) {
					// Compute right index in array
					dataBase[(i*numLabels*3 + j*3 + k)] = position;
					// Update current position within universe of discourse
					position += distToNextPos;
				}
				// MF of next label starts at center point of previous label --> Go to previous definition point
				position -= 2*distToNextPos;
			}
		}
	}

	/**
	 * Creates an initial rule base by applying the Wang-Mendel algorithm.
	 * @throws Exception 
	 */
	public void createRuleBase() throws Exception {
		// Create Wang-Mendel RB
		ArrayList<String[]> ruleBaseRaw = new ArrayList<String[]>();
		ArrayList<Double> ruleDegrees = new ArrayList<Double>();
		
		List<String> labels = sc.getLabels3();
		if(numLabels == 5) {
			labels = sc.getLabels5();
		}

		for(Date d : invPositions.keySet()) {
			// DEMAND
			// 1) Extract relevant antecedents from data
			double demand;
			if(isDemand) {
				demand = demands.get(d);
			} else {
				demand = demandChanges.get(d);
			}
			// 2) Retrieve fuzzified values
			double[] fuzzDemand;
			if(isDemand) {
				fuzzDemand = this.getFuzzified("demand", demand);
			} else {
				fuzzDemand = this.getFuzzified("demandChange", demand);
			}
			// 3) Retrieve antecedents of the rule: Choose by highest membership degree
			int maxAtDemand = 0;
			for (int i=0; i<fuzzDemand.length; i++) {
			    maxAtDemand = fuzzDemand[i] > fuzzDemand[maxAtDemand] ? i : maxAtDemand;
			}
			String antecedentDemand = labels.get(maxAtDemand);
		
			// INVENTORY POSITION
			// 1) Extract relevant antecedents from data
			double iPositions = invPositions.get(d);
			// 2) Retrieve fuzzified values
			double[] fuzzInvPosition = this.getFuzzified("inventoryPosition", iPositions);
			// 3) Retrieve antecedents of the rule: Choose by highest membership degree
			int maxAtInvPosition = 0;
			for (int i=0; i<fuzzInvPosition.length; i++) {
				maxAtInvPosition = fuzzInvPosition[i] > fuzzInvPosition[maxAtInvPosition] ? i : maxAtInvPosition;
			}
			String antecedentIlLevel = labels.get(maxAtInvPosition);
			
			
			// LEAD TIME (only in case of plants (PLANT) or warehouses (WH))
			double leadTime = 0.0;
			double[] fuzzLeadTime = new double[numLabels];
			int maxAtLeadTime = 0;
			String antecedentLeadTime = "";
			if(inventory.getStageType().equals("PLANT") || inventory.getStageType().equals("WH")) {
				// 1) Extract relevant antecedents from data
				leadTime = leadTimes.get(d);
				// 2) Retrieve fuzzified values
				fuzzLeadTime = this.getFuzzified("leadTime", leadTime);
				// 3) Retrieve antecedents of the rule: Choose by highest membership degree
				for (int i=0; i<fuzzLeadTime.length; i++) {
					maxAtLeadTime = fuzzLeadTime[i] > fuzzLeadTime[maxAtLeadTime] ? i : maxAtLeadTime;
				}
				antecedentLeadTime = labels.get(maxAtLeadTime);
			}
			
			// MARKET PRICE (only in case of plants (PLANT))
			double marketPrice = 0.0;
			double[] fuzzMarketPrice = new double[numLabels];
			int maxAtMarketPrice = 0;
			String antecedentMarketPrice = "";
			if(inventory.getStageType().equals("PLANT")) {
				// 1) Extract relevant antecedents from data
				marketPrice = marketPrices.get(d);
				// 2) Retrieve fuzzified values
				fuzzMarketPrice = this.getFuzzified("marketPrice", marketPrice);
				// 3) Retrieve antecedents of the rule: Choose by highest membership degree
				for (int i=0; i<fuzzMarketPrice.length; i++) {
					maxAtMarketPrice = fuzzMarketPrice[i] > fuzzMarketPrice[maxAtMarketPrice] ? i : maxAtMarketPrice;
				}
				antecedentMarketPrice = labels.get(maxAtMarketPrice);
			}
			
			// ORDER QUANTITY
			// 1) Extract relevant antecedents from data
			double orderQnt = orderQnts.get(d);
			// 2) Retrieve fuzzified values
			double[] fuzzOrderQnt = this.getFuzzified("orderQuantity", orderQnt);		
			// 3) Retrieve antecedents of the rule: Choose by highest membership degree
			int maxAtOrderQnt = 0;
			for (int i=0; i<fuzzOrderQnt.length; i++) {
				maxAtOrderQnt = fuzzOrderQnt[i] > fuzzOrderQnt[maxAtOrderQnt] ? i : maxAtOrderQnt;
			}
			String antecedentOrderQnt = labels.get(maxAtOrderQnt);	
			
			// Create rule
			String[] rule = new String[inventory.getNumVars()];
			if(inventory.getStageType().equals("PLANT")) {
				rule[0] = antecedentDemand;
				rule[1] = antecedentIlLevel;
				rule[2] = antecedentLeadTime;
				rule[3] = antecedentMarketPrice;
				rule[4] = antecedentOrderQnt;
			} else if(inventory.getStageType().equals("WH")) {
				rule[0] = antecedentDemand;
				rule[1] = antecedentIlLevel;
				rule[2] = antecedentLeadTime;
				rule[3] = antecedentOrderQnt;
			} else {
				rule[0] = antecedentDemand;
				rule[1] = antecedentIlLevel;
				rule[2] = antecedentOrderQnt;
			}
			ruleBaseRaw.add(rule);
			// Compute the degree of the rule
			double degree;
			if(inventory.getStageType().equals("PLANT")) {
				degree = fuzzDemand[maxAtDemand] * fuzzInvPosition[maxAtInvPosition] * fuzzLeadTime[maxAtLeadTime] * fuzzMarketPrice[maxAtMarketPrice] * fuzzOrderQnt[maxAtOrderQnt];
			} else if(inventory.getStageType().equals("WH")) {
				degree = fuzzDemand[maxAtDemand] * fuzzInvPosition[maxAtInvPosition] * fuzzLeadTime[maxAtLeadTime] * fuzzOrderQnt[maxAtOrderQnt];
			} else {
				degree = fuzzDemand[maxAtDemand] * fuzzInvPosition[maxAtInvPosition] * fuzzOrderQnt[maxAtOrderQnt];
			}
			// Add degree to list
			ruleDegrees.add(degree);
		}		
		// Choose unique rules with the highest degree and add to the rule base
		ArrayList<List<String>> addedRules = new ArrayList<List<String>>();
		ruleBase.addAll(this.getUniqueRules(ruleBaseRaw, ruleDegrees, addedRules));	
		// Add this WMGeneration object to the SC
		sc.getWmgs().put(tuple, this);
	}

	/**
	 * Returns a rule base containing the rule with the highest degree of importance for each partition.
	 * @param rulesList List of preliminary rules
	 * @param ruleDegreesList List of degrees of preliminary rules
	 * @param addedRulesList Empty list to keep track of already defined rules
	 * @return
	 */
	public List<String[]> getUniqueRules(ArrayList<String[]> rulesList, ArrayList<Double> ruleDegreesList, ArrayList<List<String>> addedRulesList) {
		List<String[]> uniqueRules = new ArrayList<String[]>();
		for(int i=0; i<rulesList.size(); i++) {
			String[] rule1 = rulesList.get(i);
			// Only compare antecedents
			List<String> rule1List = new LinkedList<String>(Arrays.asList(rule1));
			rule1List.remove(rule1List.size()-1);
			if(!addedRulesList.contains(rule1List)) {
				double maxDegree = ruleDegreesList.get(i);
				int indexMax = i;
				for(int j=i+1; j<rulesList.size(); j++) {
					if(i != j) {
						String[] rule2 = rulesList.get(j);
						List<String> rule2List = new LinkedList<String>(Arrays.asList(rule2));
						rule2List.remove(rule2List.size()-1);
						if(rule1List.equals(rule2List)) {
							double degree2 = ruleDegreesList.get(j);
							if(degree2 > maxDegree) {
								maxDegree = degree2;
								indexMax = j;
							}
						}
					}
				}
				addedRulesList.add(rule1List);
				uniqueRules.add(rulesList.get(indexMax));
			}
		}
		return uniqueRules;
	}
	
	/**
	 * Fuzzifies a given value --> Returns an array containing the memberships of a given value to the fuzzy sets of a given variable 'varName'.
	 * @param varName One variable name out of: "demand", "inventory position", "lead time", "supplier reliability", "market price", "order quantity"
	 * @param x Input value
	 * @return array that contains the memberships of a value to each fuzzy set
	 */
	public double[] getFuzzified(String varName, double x) {
		double[] memberships = new double[numLabels];
		// Retrieve data base
		double[] dB = dataBase;
		// Retrieve starting point in data base array
		int numVar = inventory.getNumVars();
		int startPos = 0;
		// Go to relevant variable in data base
		if(numVar == 3) {
			startPos = sc.getVarNamesS().indexOf(varName) * numLabels * 3;
			
		} else if(numVar == 4) {
			startPos = sc.getVarNamesM().indexOf(varName) * numLabels * 3;
			
		} else {
			startPos = sc.getVarNamesL().indexOf(varName) * numLabels * 3;
		}
		// Compute degree of membership for each MF 
		for(int i=0; i<numLabels; i++) {
			// Retrieve the 3 definition points of MF
			double a = dB[startPos + i*3 + 0];
			double m = dB[startPos + i*3 + 1];
			double b = dB[startPos + i*3 + 2];
			// Compute degree of membership according to input
			double membership = 0.0;
			if(a < x && x <= m) {
				membership = (x-a) / (m-a);
			}
			if(m < x && x < b) {
				membership = (b-x) / (b-m);
			}
			memberships[i] = membership;
		}
		return memberships;
	}
	
	/**
	 * Prints rule base. IF x_i IS ... THEN order quantity IS ...
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
	 * Prints data base.
	 */
	public void printDb() {
		for(Double a : dataBase) {
			System.out.print(a + " | ");
		}
		System.out.println("");
	}
	
	/** Creates knowledge base obtained from Wang-Mendel method.
	 * @param multiEchelon
	 * @throws Exception 
	 */
	public void createKnowledgeBase(boolean multiEchelon) throws Exception {
		wmKb = new KnowledgeBase(sc, inventory, material, dataBase, ruleBase, multiEchelon, 1.0, 1.0);
	}
	
	/** 
	 * @return the knowledge base created by applying Wang-Mendel algorithm
	 * @throws Exception
	 */
	public KnowledgeBase getWmKnowledgeBase() throws Exception {
		return wmKb;
	}

	/**
     * Returns difference between two dates in preferred time unit.
     * @param date1 the newer date
     * @param date2 the older date
     * @param timeUnit the unit in which the difference should be expressed
     * @return the difference between the two dates
     */
    public long getDateDiff(Date date1, Date date2, TimeUnit timeUnit) {
        long diffInMillies = date1.getTime() - date2.getTime();
        return timeUnit.convert(diffInMillies,TimeUnit.MILLISECONDS);
    }
    
	/** Returns a list with the universe of discourses of all fuzzy variables.
	 * @param material
	 * @param demand Indicates if demand or demand change is used
	 * @return
	 */
	public List<double[]> getRanges(boolean isDemand) {
		List<double[]> ranges = new ArrayList<double[]>();
		double[] dRange;
		if(isDemand) {
			dRange = this.getDemandRange();
		} else {
			dRange = this.getDemandChangeRange();
		}
		double[] ipRange = this.getInventoryPositionRange();
		double[] oqRange = this.getOrderQuantityRange();
		
		// Add ranges to list
		ranges.add(dRange);
		ranges.add(ipRange);
		if(inventory.getStageType().equals("PLANT") || inventory.getStageType().equals("WH")) {
			ranges.add(this.getLeadTimeRange());
		}
		// Add market price range if this inventory is a plant
		if(inventory.getStageType().equals("PLANT")) {
			ranges.add(this.getMarketPriceRange());
		}
		ranges.add(oqRange);
		
		// Return list of ranges
		return ranges;
	}
	
	/** Returns universe of discourse of order quantity.
	 * @return
	 */
	public double[] getOrderQuantityRange() {
		int max = 0;
		for (var entry : orderQnts.entrySet()) {
			if(entry.getValue() > max && !entry.getKey().after(sc.getSplitDate())) {
				max = entry.getValue();
			}
		}
		double[] range = new double[2];
		range[0] = 0;
		range[1] = max;
		return range;
	}
	
	/** Returns universe of discourse of variable lead time.
	 * @return
	 */
	public double[] getLeadTimeRange() {
		int min = Integer.MAX_VALUE;
		int max = 0;
		for (var entry : leadTimes.entrySet()) {
			if (entry.getValue() < min && !entry.getKey().after(sc.getSplitDate())) {
				min = entry.getValue();
			} else if(entry.getValue() > max && !entry.getKey().after(sc.getSplitDate())) {
				max = entry.getValue();
			}
		}
		double[] range = new double[2];
		range[0] = min;
		range[1] = max;
		return range;
	}
	
	/** Returns universe of discourse of variable demand.
	 * @return
	 */
	public double[] getDemandRange() {
		int min = 0;
		int max = 0;
		for (var entry : demands.entrySet()) {
			if(entry.getValue() > max && !entry.getKey().after(sc.getSplitDate())) {
				max = entry.getValue();
			}
		}
		double[] range = new double[2];
		range[0] = min;
		range[1] = max;
		return range;
	}
	
	/** Returns universe of discourse of variable demand change.
	 * @return
	 */
	public double[] getDemandChangeRange() {
		int min = Integer.MAX_VALUE;
		int max = 0;
		for (var entry : demandChanges.entrySet()) {
			if (entry.getValue() < min && !entry.getKey().after(sc.getSplitDate())) {
				min = entry.getValue();
			} else if(entry.getValue() > max && !entry.getKey().after(sc.getSplitDate())) {
				max = entry.getValue();
			}
		}
		double[] range = new double[2];
		range[0] = min;
		range[1] = max;
		return range;
	}
	
	/** Returns universe of discourse of variable inventory position.
	 * @return
	 */
	public double[] getInventoryPositionRange() {
		int min = Integer.MAX_VALUE;
		int max = 0;
		for (var entry : invPositions.entrySet()) {
			if (entry.getValue() > max && !entry.getKey().after(sc.getSplitDate())) {
				max = entry.getValue();
			} else if(entry.getValue() < min && !entry.getKey().after(sc.getSplitDate())) {
				min = entry.getValue();
			}
		}
		double[] range = new double[2];
		range[0] = 0;
		range[1] = max;
		return range;
	}
	
	/** Returns universe of discourse of variable market price.
	 * @return
	 */
	public double[] getMarketPriceRange() {
		double min = Double.POSITIVE_INFINITY;
		double max = 0.0;
		for (var entry : marketPrices.entrySet()) {
			if (entry.getValue() < min && !entry.getKey().after(sc.getSplitDate())) {
				min = entry.getValue();
			} else if(entry.getValue() > max && !entry.getKey().after(sc.getSplitDate())) {
				max = entry.getValue();
			}
		}
		double[] range = new double[2];
		range[0] = min;
		range[1] = max;
		return range;
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
}