import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * This class offers a collection of methods for writing information into .txt files.
 *
 */
public class Writer {
	private SC sc;
	private File file;

    /**
     * @param sc SC object
     */
    public Writer(SC sc) {
    	this.sc = sc;
    }
    
    /**
     * Writes run time (in seconds) of a method into separate .txt file.
     * @param runTime 
     * @param heuristic is 1 if heuristic and 0 if exact model
     * @param scenario Indicates which scenarios was solved
     */
    public void writeRuntime(long runTime, boolean heuristic, int scenario) {
        try {
        	String modelType = "";
        	if(heuristic) {
        		modelType = "Heuristic";
        	} else {
        		modelType = "Exact";
        	}
            FileWriter writer = new FileWriter("runtime" + modelType + "Scenario" + scenario +".txt", true);
            writer.write("Runtime:");
            writer.write("\n");
            writer.write((runTime / (Math.pow(10, 9))) + " ");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Writes data bases and rule bases of a given knowledge base set into a .txt file.
     */
    public void writeDBsAndRBs(boolean heuristic, KnowledgeBaseSet kbSet, int scenario) {
    	// One .txt file per KB
    	for (var entry : kbSet.getKnowledgeBases().entrySet()) {
    		// Data base 
            try {
            	String modelType = "";
            	if(heuristic) {
            		modelType = "Heuristic";
            	} else {
            		modelType = "Exact";
            	}
                FileWriter writer = new FileWriter("dataBase" + entry.getKey() + modelType + "Scenario" + scenario +".txt", true);
                writer.write("MFParameters");
                writer.write("\n");
                double[] dataBase = entry.getValue().getDataBase();
                for (int i=0; i<dataBase.length; i++) {
                	writer.write(dataBase[i] + " ");
                }
                writer.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
            // Rule base
            try {
            	String modelType = "";
            	if(heuristic) {
            		modelType = "Heuristic";
            	} else {
            		modelType = "Exact";
            	}
                FileWriter writer = new FileWriter("ruleBase" + entry.getKey() + modelType + "Scenario" + scenario +".txt", true);
                writer.write("Rules");
                writer.write("\n");
                HashSet<String[]> ruleBase = entry.getValue().getRuleBase();
                for (String[] rule : ruleBase) {
                	for(int i=0; i<rule.length; i++) {
                		writer.write(rule[i] + " ");
                	}
                	writer.write("\n");
                }
                writer.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
    	}
    }
    
    /**
     * Writes the final training results (overall and per inventory) of a model into a .txt file.
     * @param kbSet Set of knowledge bases (if kbSet is null --> data generator not null and vice versa)
     * @param dataGeneration
     * @param modelType EOQ (0), Exact (1) or Heuristic (2)
     * @param scenario
     */
    public void writeFinalResultsTraining(KnowledgeBaseSet kbSet, DataGeneration dataGeneration, int modelType, int scenario) {
    	// Overall results 
    	try {
    		String modelTypeStr = "EOQ";
    		if(modelType == 1) {
    			modelTypeStr = "Exact";
    		}
    		if(modelType == 2) {
    			modelTypeStr = "Heuristic";
    		}
    		if(modelType == 3) {
    			modelTypeStr = "WM";
    		}
    		FileWriter writer = new FileWriter("finalResults" + modelTypeStr + "Training" + "Scenario" + scenario + ".txt", true);
    		writer.write("ScCost FR holdingCosts penaltyCosts orderCosts setupCosts productionCosts purchaseCosts transportCosts");
    		writer.write("\n");
    		if(kbSet != null && dataGeneration == null) {
    			Simulation simulation = kbSet.getSimulation();
    			writer.write(simulation.getTotalScCosts() + " ");
    			writer.write(simulation.getSCFillRate() + " ");
    			writer.write(simulation.getTotalHoldingCosts() + " ");
    			writer.write(simulation.getTotalPenaltyCosts() + " ");
    			writer.write(simulation.getTotalOrderCosts() + " ");
    			writer.write(simulation.getTotalSetupCosts() + " ");
    			writer.write(simulation.getTotalProductionCosts() + " ");
    			writer.write(simulation.getTotalUnitCosts() + " ");
    			writer.write(simulation.getTotalTransportationCosts() + " ");
    		}
    		if(dataGeneration != null && kbSet == null) {
    			writer.write(dataGeneration.getTrainingTotalScCosts() + " ");
    			writer.write(dataGeneration.getTrainingSCFillRate() + " ");
    			writer.write(dataGeneration.getTrainingTotalHoldingCosts() + " ");
    			writer.write(dataGeneration.getTrainingTotalPenaltyCosts() + " ");
    			writer.write(dataGeneration.getTrainingTotalOrderCosts() + " ");
    			writer.write(dataGeneration.getTrainingTotalSetupCosts() + " ");
    			writer.write(dataGeneration.getTrainingTotalProductionCosts() + " ");
    			writer.write(dataGeneration.getTrainingTotalUnitCosts() + " ");
    			writer.write(dataGeneration.getTrainingTotalTransportationCosts() + " ");
    		}
    		writer.close();
    	} catch (IOException e) {
    		e.printStackTrace();
    	}

    	// Results per inventory and material
    	try {
    		String modelTypeStr = "EOQ";
    		if(modelType == 1) {
    			modelTypeStr = "Exact";
    		}
    		if(modelType == 2) {
    			modelTypeStr = "Heuristic";
    		}
    		if(modelType == 3) {
    			modelTypeStr = "WM";
    		}
    		FileWriter writer = new FileWriter("finalResultsPerInventory" + modelTypeStr + "Training" + "Scenario" + scenario + ".txt", true);
    		writer.write("Tuple ScCost FR holdingCosts penaltyCosts orderCosts setupCosts productionCosts purchaseCosts transportCosts");
    		writer.write("\n");
    		for(Inventory inventory : sc.getInventories()) {
    			if(inventory.isActiveInv()) {
    				for(Material material : sc.getMaterials()) {
    					String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
    					if(kbSet != null && dataGeneration == null) {
    						Simulation simulation = kbSet.getSimulation();
    						writer.write(tuple + " ");
    						writer.write(simulation.getTotalCosts().get(tuple) + " ");
    						writer.write(simulation.getFillRate(tuple) + " ");
    						writer.write(simulation.getHoldingCostMap().get(tuple) + " ");
    						writer.write(simulation.getPenaltyCostMap().get(tuple) + " ");
    						writer.write(simulation.getOrderCostMap().get(tuple) + " ");
    						writer.write(simulation.getSetupCostMap().get(tuple) + " ");
    						writer.write(simulation.getProductionCostMap().get(tuple) + " ");
    						writer.write(simulation.getUnitCostMap().get(tuple) + " ");
    						writer.write(simulation.getTransportationCostMap().get(tuple) + " ");
    						writer.write("\n");
    					}
    					if(dataGeneration != null && kbSet == null) {
    						writer.write(tuple + " ");
    						writer.write(dataGeneration.getTotalCostsTraining().get(tuple) + " ");
    						writer.write(dataGeneration.getTrainingFillRate(tuple) + " ");
    						writer.write(dataGeneration.getHoldingCostMapTraining().get(tuple) + " ");
    						writer.write(dataGeneration.getPenaltyCostMapTraining().get(tuple) + " ");
    						writer.write(dataGeneration.getOrderCostMapTraining().get(tuple) + " ");
    						writer.write(dataGeneration.getSetupCostMapTraining().get(tuple) + " ");
    						writer.write(dataGeneration.getProductionCostMapTraining().get(tuple) + " ");
    						writer.write(dataGeneration.getUnitCostMapTraining().get(tuple) + " ");
    						writer.write(dataGeneration.getTransportationCostMapTraining().get(tuple) + " ");
    						writer.write("\n");
    					}
    				}
    			}
    		}
    		writer.close();
    	} catch (IOException e) {
    		e.printStackTrace();
    	}
    }
    
    /**
     * Writes the final test results (overall and per inventory) of a model into a .txt file.
     * @param kbSet Set of knowledge bases (if kbSet is null --> data generator not null and vice versa)
     * @param dataGeneration
     * @param modelType EOQ (0), Exact (1) or Heuristic (2)
     * @param scenario
     */
    public void writeFinalResultsTest(Simulation simulationFIS, DataGeneration dataGeneration, int modelType, int scenario) {
        try {
        	String modelTypeStr = "EOQ";
        	if(modelType == 1) {
        		modelTypeStr = "Exact";
        	}
        	if(modelType == 2) {
        		modelTypeStr = "Heuristic";
        	}
            FileWriter writer = new FileWriter("finalResults" + modelTypeStr + "Test" + "Scenario" + scenario + ".txt", true);
            writer.write("ScCost FR holdingCosts penaltyCosts orderCosts setupCosts productionCosts purchaseCosts transportCosts");
            writer.write("\n");
            if(simulationFIS != null && dataGeneration == null) {
                writer.write(simulationFIS.getTotalScCosts() + " ");
                writer.write(simulationFIS.getSCFillRate() + " ");
                writer.write(simulationFIS.getTotalHoldingCosts() + " ");
                writer.write(simulationFIS.getTotalPenaltyCosts() + " ");
                writer.write(simulationFIS.getTotalOrderCosts() + " ");
                writer.write(simulationFIS.getTotalSetupCosts() + " ");
                writer.write(simulationFIS.getTotalProductionCosts() + " ");
                writer.write(simulationFIS.getTotalUnitCosts() + " ");
                writer.write(simulationFIS.getTotalTransportationCosts() + " ");
            }
            if(dataGeneration != null && simulationFIS == null) {
                writer.write(dataGeneration.getTestTotalScCosts() + " ");
                writer.write(dataGeneration.getTestSCFillRate() + " ");
                writer.write(dataGeneration.getTestTotalHoldingCosts() + " ");
                writer.write(dataGeneration.getTestTotalPenaltyCosts() + " ");
                writer.write(dataGeneration.getTestTotalOrderCosts() + " ");
                writer.write(dataGeneration.getTestTotalSetupCosts() + " ");
                writer.write(dataGeneration.getTestTotalProductionCosts() + " ");
                writer.write(dataGeneration.getTestTotalUnitCosts() + " ");
                writer.write(dataGeneration.getTestTotalTransportationCosts() + " ");
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        
        // Write results for each inventory and material
    	try {
    		String modelTypeStr = "EOQ";
    		if(modelType == 1) {
    			modelTypeStr = "Exact";
    		}
    		if(modelType == 2) {
    			modelTypeStr = "Heuristic";
    		}
    		FileWriter writer = new FileWriter("finalResultsPerInventory" + modelTypeStr + "Test" + "Scenario" + scenario + ".txt", true);
    		writer.write("Tuple ScCost FR holdingCosts penaltyCosts orderCosts setupCosts productionCosts purchaseCosts transportCosts");
    		writer.write("\n");
    		for(Inventory inventory : sc.getInventories()) {
    			if(inventory.isActiveInv()) {
    				for(Material material : sc.getMaterials()) {
    					String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
    					if(simulationFIS != null && dataGeneration == null) {
    						writer.write(tuple + " ");
    						writer.write(simulationFIS.getTotalCosts().get(tuple) + " ");
    						writer.write(simulationFIS.getFillRate(tuple) + " ");
    						writer.write(simulationFIS.getHoldingCostMap().get(tuple) + " ");
    						writer.write(simulationFIS.getPenaltyCostMap().get(tuple) + " ");
    						writer.write(simulationFIS.getOrderCostMap().get(tuple) + " ");
    						writer.write(simulationFIS.getSetupCostMap().get(tuple) + " ");
    						writer.write(simulationFIS.getProductionCostMap().get(tuple) + " ");
    						writer.write(simulationFIS.getUnitCostMap().get(tuple) + " ");
    						writer.write(simulationFIS.getTransportationCostMap().get(tuple) + " ");
    						writer.write("\n");
    					}
    					if(dataGeneration != null && simulationFIS == null) {
    						writer.write(tuple + " ");
    						writer.write(dataGeneration.getTotalCostsTest().get(tuple) + " ");
    						writer.write(dataGeneration.getTestFillRate(tuple) + " ");
    						writer.write(dataGeneration.getHoldingCostMapTest().get(tuple) + " ");
    						writer.write(dataGeneration.getPenaltyCostMapTest().get(tuple) + " ");
    						writer.write(dataGeneration.getOrderCostMapTest().get(tuple) + " ");
    						writer.write(dataGeneration.getSetupCostMapTest().get(tuple) + " ");
    						writer.write(dataGeneration.getProductionCostMapTest().get(tuple) + " ");
    						writer.write(dataGeneration.getUnitCostMapTest().get(tuple) + " ");
    						writer.write(dataGeneration.getTransportationCostMapTest().get(tuple) + " ");
    						writer.write("\n");
    					}
    				}
    			}
    		}
    		writer.close();
    	} catch (IOException e) {
    		e.printStackTrace();
    	}
    }
    
    /**
     * Writes training inventory levels (per date) for each inventory/material into .txt file. 
     * For each inventory/material there is one row (Date1, Date2, ..., DateN)
     * @param simulationFIS
     * @param dataGeneration
     * @param modelType EOQ (0), Exact (1) or Heuristic (2)
     */
    public void writeInventoryLevelsTraining(Simulation simulationFIS, DataGeneration dataGeneration, int modelType) {
        // Write results for each inventory and material
    	try {
    		String modelTypeStr = "EOQ";
    		if(modelType == 1) {
    			modelTypeStr = "Exact";
    		}
    		if(modelType == 2) {
    			modelTypeStr = "Heuristic";
    		}
    		FileWriter writer = new FileWriter("inventoryLevels" + modelTypeStr + "Training" + ".txt", true);
    		writer.write("Tuple ");
    		for(Date trainingDate: sc.getTrainingDates()) {
    			writer.write(trainingDate + " ");
    		}
    		writer.write("\n");
    		for(Inventory inventory : sc.getInventories()) {
    			if(inventory.isActiveInv()) {
    				for(Material material : inventory.getMaterials()) {
    					String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
    					writer.write(tuple + " ");
    					if(simulationFIS != null && dataGeneration == null) {
    						TreeMap<Date,Integer> inventoryLevels = simulationFIS.getInventoryLevelsMap().get(tuple);
    						for (var entry : inventoryLevels.entrySet()) {
    							if(!entry.getKey().after(sc.getSplitDate())) {
    								writer.write(entry.getValue() + " ");
    							}
    						}
    						writer.write("\n");
    					}
    					if(dataGeneration != null && simulationFIS == null) {
    						TreeMap<Date,Integer> inventoryLevels = dataGeneration.getInventoryLevelsMap().get(tuple);
    						for (var entry : inventoryLevels.entrySet()) {
    							if(!entry.getKey().after(sc.getSplitDate())) {
    								writer.write(entry.getValue() + " ");
    							}
    						}
    						writer.write("\n");
    					}
    				}
    			}
    		}
    		writer.close();
    	} catch (IOException e) {
    		e.printStackTrace();
    	}
    }
    
    /**
     * Writes test inventory levels (per date) for each inventory/material into .txt file. 
     * For each inventory/material there is one row (Date1, Date2, ..., DateN)
     * @param simulationFIS
     * @param dataGeneration
     * @param modelType EOQ (0), Exact (1) or Heuristic (2)
     */
    public void writeInventoryLevelsTest(Simulation simulationFIS, DataGeneration dataGeneration, int modelType) {
        // Write results for each inventory and material
    	try {
    		String modelTypeStr = "EOQ";
    		if(modelType == 1) {
    			modelTypeStr = "Exact";
    		}
    		if(modelType == 2) {
    			modelTypeStr = "Heuristic";
    		}
    		FileWriter writer = new FileWriter("inventoryLevels" + modelTypeStr + "Test" + ".txt", true);
    		writer.write("Tuple ");
    		for(Date testDate: sc.getTestDates()) {
    			writer.write(testDate + " ");
    		}
    		writer.write("\n");
    		for(Inventory inventory : sc.getInventories()) {
    			if(inventory.isActiveInv()) {
    				for(Material material : inventory.getMaterials()) {
    					String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
    					writer.write(tuple + " ");
    					if(simulationFIS != null && dataGeneration == null) {
    						TreeMap<Date,Integer> inventoryLevels = simulationFIS.getInventoryLevelsMap().get(tuple);
    						for (var entry : inventoryLevels.entrySet()) {
    							if(entry.getKey().after(sc.getSplitDate())) {
    								writer.write(entry.getValue() + " ");
    							}
    						}
    						writer.write("\n");
    					}
    					if(dataGeneration != null && simulationFIS == null) {
    						TreeMap<Date,Integer> inventoryLevels = dataGeneration.getInventoryLevelsMap().get(tuple);
    						for (var entry : inventoryLevels.entrySet()) {
    							if(entry.getKey().after(sc.getSplitDate())) {
    								writer.write(entry.getValue() + " ");
    							}
    						}
    						writer.write("\n");
    					}
    				}
    			}
    		}
    		writer.close();
    	} catch (IOException e) {
    		e.printStackTrace();
    	}
    }
    
    /**
     * Writes training inventory positions (per date) for each inventory/material into .txt file. 
     * For each inventory/material there is one row (Date1, Date2, ..., DateN)
     * @param simulationFIS
     * @param dataGeneration
     * @param modelType EOQ (0), Exact (1) or Heuristic (2)
     */
    public void writeInventoryPositionsTraining(Simulation simulationFIS, DataGeneration dataGeneration, int modelType) {
        // Write results for each inventory and material
    	try {
    		String modelTypeStr = "EOQ";
    		if(modelType == 1) {
    			modelTypeStr = "Exact";
    		}
    		if(modelType == 2) {
    			modelTypeStr = "Heuristic";
    		}
    		FileWriter writer = new FileWriter("inventoryPositions" + modelTypeStr + "Training" + ".txt", true);
    		writer.write("Tuple ");
    		for(Date trainingDate : sc.getTrainingDates()) {
    			writer.write(trainingDate + " ");
    		}
    		writer.write("\n");
    		for(Inventory inventory : sc.getInventories()) {
    			if(inventory.isActiveInv()) {
    				for(Material material : inventory.getMaterials()) {
    					String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
    					writer.write(tuple + " ");
    					if(simulationFIS != null && dataGeneration == null) {
    						TreeMap<Date,Integer> inventoryPositions = simulationFIS.getInventoryPositionsMap().get(tuple);
    						for (var entry : inventoryPositions.entrySet()) {
    							if(!entry.getKey().after(sc.getSplitDate())) {
    								writer.write(entry.getValue() + " ");
    							}
    						}
    						writer.write("\n");
    					}
    					if(dataGeneration != null && simulationFIS == null) {
    						TreeMap<Date,Integer> inventoryPositions = dataGeneration.getInventoryPositionsMap().get(tuple);
    						for (var entry : inventoryPositions.entrySet()) {
    							if(!entry.getKey().after(sc.getSplitDate())) {
    								writer.write(entry.getValue() + " ");
    							}
    						}
    						writer.write("\n");
    					}
    				}
    			}
    		}
    		writer.close();
    	} catch (IOException e) {
    		e.printStackTrace();
    	}
    }
    
    /**
     * Writes test inventory positions (per date) for each inventory/material into .txt file. 
     * For each inventory/material there is one row (Date1, Date2, ..., DateN)
     * @param simulationFIS
     * @param dataGeneration
     * @param modelType EOQ (0), Exact (1) or Heuristic (2)
     */
    public void writeInventoryPositionsTest(Simulation simulationFIS, DataGeneration dataGeneration, int modelType) {
        // Write results for each inventory and material
    	try {
    		String modelTypeStr = "EOQ";
    		if(modelType == 1) {
    			modelTypeStr = "Exact";
    		}
    		if(modelType == 2) {
    			modelTypeStr = "Heuristic";
    		}
    		FileWriter writer = new FileWriter("inventoryPositions" + modelTypeStr + "Test" + ".txt", true);
    		writer.write("Tuple ");
    		for(Date testDate : sc.getTestDates()) {
    			writer.write(testDate + " ");
    		}
    		writer.write("\n");
    		for(Inventory inventory : sc.getInventories()) {
    			if(inventory.isActiveInv()) {
    				for(Material material : inventory.getMaterials()) {
    					String tuple = "" + inventory.getInventoryId() + "," + material.getMaterialCode();
    					writer.write(tuple + " ");
    					if(simulationFIS != null && dataGeneration == null) {
    						TreeMap<Date,Integer> inventoryPositions = simulationFIS.getInventoryPositionsMap().get(tuple);
    						for (var entry : inventoryPositions.entrySet()) {
    							if(entry.getKey().after(sc.getSplitDate())) {
    								writer.write(entry.getValue() + " ");
    							}
    						}
    						writer.write("\n");
    					}
    					if(dataGeneration != null && simulationFIS == null) {
    						TreeMap<Date,Integer> inventoryPositions = dataGeneration.getInventoryPositionsMap().get(tuple);
    						for (var entry : inventoryPositions.entrySet()) {
    							if(entry.getKey().after(sc.getSplitDate())) {
    								writer.write(entry.getValue() + " ");
    							}
    						}
    						writer.write("\n");
    					}
    				}
    			}
    		}
    		writer.close();
    	} catch (IOException e) {
    		e.printStackTrace();
    	}
    }
    
    /** Writes the interim results after each step of the exact 3-step-procedure into a .txt file.
     * @param kbSetWM kbSet after Wang Mendel method
     * @param kbSetSA kbSet after Simulated Annealing
     * @param kbSetLA kbSet after LA-Tuning algorithm (genetic algorithm)
     * @param scenario
     */
    public void writeInterimTrainingResultsExact(KnowledgeBaseSet kbSetWM, KnowledgeBaseSet kbSetSA, KnowledgeBaseSet kbSetLA, int scenario) {
        try {
            FileWriter writer = new FileWriter("InterimTrainingResultsExact" + "Scenario" + scenario + ".txt", true);
            writer.write("FitnessValue ScCosts FR");
            writer.write("\n");
            writer.write(kbSetWM.getFitnessValue() + " " + kbSetWM.getObjVal() + " " + kbSetWM.getFillRate());
            writer.write("\n");
            writer.write(kbSetSA.getFitnessValue() + " " + kbSetSA.getObjVal() + " " + kbSetSA.getFillRate());
            writer.write("\n");
            writer.write(kbSetLA.getFitnessValue() + " " + kbSetLA.getObjVal() + " " + kbSetLA.getFillRate());
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    /** Writes the interim results after each step of the local (heuristic) 3-step-procedure into a .txt file.
     * @param kbSetWM kbSet after Wang-Mendel method
     * @param kbSetSA kbSet after Simulated Annealing
     * @param kbSetLA kbSet after LA-Tuning algorithm (genetic algorithm)
     * @param scenario
     */
    public void writeInterimTrainingResultsLocal(String tuple, KnowledgeBase kbWM, KnowledgeBase kbSA, KnowledgeBase kbLA, int scenario) {
        try {
            FileWriter writer = new FileWriter("InterimTrainingResultsLocal" + tuple + "Scenario" + scenario + ".txt", true);
            writer.write("FitnessValue ScCosts FR");
            writer.write("\n");
            writer.write(kbWM.getFitnessValue() + " " + kbWM.getObjVal() + " " + kbWM.getFillRate());
            writer.write("\n");
            writer.write(kbSA.getFitnessValue() + " " + kbSA.getObjVal() + " " + kbSA.getFillRate());
            writer.write("\n");
            writer.write(kbLA.getFitnessValue() + " " + kbLA.getObjVal() + " " + kbLA.getFillRate());
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    /** Writes results of global (exact method) LA tuning algorithm per iteration into a .txt file.
     * @param fitnessValues List of fitness values
     * @param objValues List of objective values
     * @param fillRates List of fill rates
     * @throws IOException 
     */
    public void writeLAIterationValuesGlobal(List<Double> fitnessValues, List<Double> objValues, List<Double> fillRates, int scenario) throws IOException {
        try {
            FileWriter writer = new FileWriter("ValuesPerIterationLA" + "Scenario" + scenario + ".txt", true);
            writer.write("Iteration FitnessValue ScCosts FR");
            writer.write("\n");
        	for(int i=0; i<fitnessValues.size(); i++) {
                writer.write(i + " " + fitnessValues.get(i) + " " + objValues.get(i) + " " + fillRates.get(i));
                writer.write("\n");
        	}
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    /** Writes results of local (heuristic) LA tuning algorithm per iteration into a .txt file.
     * @param fitnessValues List of fitness values
     * @param objValues List of objective values
     * @param fillRates List of fill rates
     * @throws IOException 
     */
    public void writeLAIterationValuesLocal(String tuple, List<Double> fitnessValues, List<Double> objValues, List<Double> fillRates, int scenario) throws IOException {
        try {
            FileWriter writer = new FileWriter("ValuesPerIterationLA" + tuple + "Scenario" + scenario + ".txt", true);
            writer.write("Iteration FitnessValue ScCosts FR");
            writer.write("\n");
        	for(int i=0; i<fitnessValues.size(); i++) {
                writer.write(i + " " + fitnessValues.get(i) + " " + objValues.get(i) + " " + fillRates.get(i));
                writer.write("\n");
        	}
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Writes a .txt file with all training input data.
     */
    public void writeTrainingInputData(Inventory inventory, Material material, TreeMap<Date,Integer> demands, TreeMap<Date,Integer> demandChanges,
    		TreeMap<Date,Integer> invPositions, TreeMap<Date,Integer> leadTimes, TreeMap<Date,Double> marketPrices, TreeMap<Date,Integer> orderQuantities) {
        try {
            FileWriter writer = new FileWriter("inputData" + inventory.getInventoryId() + "," + material.getMaterialCode() + ".txt", true);
            // Write field names
            if(inventory.getStageType().equals("PLANT")) {
            	writer.write("date demand demandChange inventoryPositions leadTime marketPrice orderQuantity");
            } else if(inventory.getStageType().equals("WH")) {
            	writer.write("date demand demandChange inventoryPositions leadTime orderQuantity");
            } else {
            	writer.write("date demand demandChange inventoryPositions orderQuantity");
            }
            
            writer.write("\n");
            // Write input data
    		for(Date d : sc.getDates()) {
    			if(!d.after(sc.getSplitDate())) {
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
    		}
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
   
    /**
     * Programmatically writes a .fcl text file containing the information (DB + RB) of a KB related to a particular inventory and material.
     * A .FCL file contains all the information to run a complete fuzzy model (Fuzzification, Fuzzy Inference, Defuzzification) (see Simulation class, method: getOrderQuantity())
     * @param inventory Inventory of interest
     * @param material Material of interest
     * @param dataBase Data base --> contains definition points of all membership functions 
     * @param ruleBase Rule base --> contains all fuzzy if-then rules
     * @param varNames String names of all fuzzy variables that can be used (e.g. demand, inventoryPosition, leadTime, marketPrice, orderQuantity, etc.)
     * @param labelNames Available names of fuzzy labels (e.g. low, medium, high, etc.)
     * @param iteration Number of iteration
     * @throws IOException 
     */
    public void writeFCL(Inventory inventory, Material material, double[] dataBase, HashSet<String[]> ruleBase, List<String> varNames, 
    	List<String> labelNames, int iteration) throws IOException {
        try {
        	String fileName = "kbGeneration" + inventory.getInventoryId() + "," + material.getMaterialCode() + ".fcl";
        	// Create file
        	FileWriter writer = new FileWriter(fileName,true);
        	// Start function block
        	String fbName = "kbGeneration" + iteration;
        	writer.write("FUNCTION_BLOCK " + fbName);
        	// 2 new lines
        	writer.write("\n");
        	writer.write("\n");
        	writer.write("VAR_INPUT");
        	// New line
        	writer.write("\n");
        	// Define input variables as real values
        	for(int i=0; i<varNames.size()-1; i++) {
        		writer.write(varNames.get(i) + " : " + "REAL;");
        		writer.write("\n");
        	}
        	writer.write("END_VAR");
        	// 2 new lines
        	writer.write("\n");
        	writer.write("\n");
        	writer.write("VAR_OUTPUT");
        	// New line
        	writer.write("\n");
        	// Define output variable as a real value
       		writer.write(varNames.get(varNames.size()-1) + " : " + "REAL;");
    		writer.write("\n");
    		writer.write("END_VAR");
        	// 2 new lines
        	writer.write("\n");
        	writer.write("\n");
        	// Define membership functions of input variables
        	for(int i=0; i<varNames.size()-1; i++) {
        		int startPos = i * labelNames.size() * 3;
        		writer.write("FUZZIFY" + " " + varNames.get(i));
        		writer.write("\n");
        		for(int j=0; j<labelNames.size(); j++) {
        			// Retrieve the 3 definition points of mf
        			double defPoint1 = dataBase[startPos + j*3 + 0];
        			double defPoint2 = dataBase[startPos + j*3 + 1];
        			double defPoint3 = dataBase[startPos + j*3 + 2];
        			String mfDefinition = "TERM" + " " + labelNames.get(j) + " := " + "TRIAN " + defPoint1 + " " + defPoint2 + " " + defPoint3 + ";";
        			writer.write(mfDefinition);
        			writer.write("\n");
        		}
        		writer.write("END_FUZZIFY");
        		writer.write("\n");
        		writer.write("\n");
        	}
        	// Define membership function of output variable
    		int startPos = (varNames.size() - 1) * labelNames.size() * 3;
    		writer.write("DEFUZZIFY" + " " + varNames.get(varNames.size()-1));
    		writer.write("\n");
    		for(int i=0; i<labelNames.size(); i++) {
    			// Retrieve the 3 definition points of MF
    			double defPoint1 = dataBase[startPos + i*3 + 0];
    			double defPoint2 = dataBase[startPos + i*3 + 1];
    			double defPoint3 = dataBase[startPos + i*3 + 2];
    			String mfDefinition = "TERM" + " " + labelNames.get(i) + " := " + "TRIAN " + defPoint1 + " " + defPoint2 + " " + defPoint3 + ";";
    			writer.write(mfDefinition);
    			writer.write("\n");
    		}
    		writer.write("METHOD : COG;");
    		writer.write("\n");
    		writer.write("DEFAULT := 0;");
    		writer.write("\n");
    		writer.write("END_DEFUZZIFY");
    		writer.write("\n");
    		writer.write("\n");
    		// Define rule block
    		writer.write("RULEBLOCK No1");
    		writer.write("\n");
    		writer.write("AND : MIN;");
    		writer.write("\n");
    		writer.write("ACT : MIN;");
    		writer.write("\n");
    		writer.write("ACCU : MAX;");
    		writer.write("\n");
    		writer.write("\n");
    		int index = 1;
    		for(String[] rule : ruleBase) {
    			String ruleStr;
    			if(varNames.size() == 5) {
        			if(sc.isDemand()) {
        				ruleStr = "IF " + "demand IS " + rule[0] + " AND inventoryPosition IS " + rule[1] + " AND leadTime IS " + rule[2] +  " AND marketPrice IS " + rule[3] + /*" AND supplierReliability IS " + rule[3] + */ " THEN orderQuantity IS " + rule[4] + ";";
        			} else {
        				ruleStr = "IF " + "demandChange IS " + rule[0] + " AND inventoryPosition IS " + rule[1] + " AND leadTime IS " + rule[2] + " AND marketPrice IS " + rule[3] + /*" AND supplierReliability IS " + rule[3] + */ " THEN orderQuantity IS " + rule[4] + ";";
        			}
    			} else if(varNames.size() == 4){
        			if(sc.isDemand()) {
        				ruleStr = "IF " + "demand IS " + rule[0] + " AND inventoryPosition IS " + rule[1] + " AND leadTime IS " + rule[2] + /*" AND supplierReliability IS " + rule[3] + */ " THEN orderQuantity IS " + rule[3] + ";";
        			} else {
        				ruleStr = "IF " + "demandChange IS " + rule[0] + " AND inventoryPosition IS " + rule[1] + " AND leadTime IS " + rule[2] + /*" AND supplierReliability IS " + rule[3] + */ " THEN orderQuantity IS " + rule[3] + ";";
        			}
    			} else {
        			if(sc.isDemand()) {
        				ruleStr = "IF " + "demand IS " + rule[0] + " AND inventoryPosition IS " + rule[1] + /*" AND supplierReliability IS " + rule[3] + */ " THEN orderQuantity IS " + rule[2] + ";";
        			} else {
        				ruleStr = "IF " + "demandChange IS " + rule[0] + " AND inventoryPosition IS " + rule[1] + /*" AND supplierReliability IS " + rule[3] + */ " THEN orderQuantity IS " + rule[2] + ";";
        			}
    			}
    			writer.write("RULE " + index + " : " + ruleStr);
    			writer.write("\n");
    			index += 1;
    		}
    		writer.write("END_RULEBLOCK");
    		writer.write("\n");
    		writer.write("\n");
    		writer.write("END_FUNCTION_BLOCK");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}