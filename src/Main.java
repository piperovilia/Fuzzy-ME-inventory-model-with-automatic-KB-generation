import java.io.File;
import java.text.SimpleDateFormat;
import java.util.*;

public class Main {

	public static void main(String[] args) throws Exception {
		// SIMULATION/SENSITIVITY ANALYSIS
		// Read data
		double splitRatio = 0.7;
		boolean isDemand = true;
		int numDates = 360;
		SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd");
		Date startDate = dateFormat.parse("2021-01-01");
		DataReader reader = new DataReader(splitRatio, isDemand);
		reader.readDataSimulation(startDate, numDates);
		SC initSc = reader.getSupplyChain();
		TreeSet<Date> dates = initSc.getDates();
		
		// Define parameters
		boolean multiEchelon = true;
		int numLabels = 3;
		boolean continuousDemand = true;
		// Parameters for Simulated Annealing
		double alpha = 0.95;
		double tempGlobal = 0.25;
		double tempLocal = 0.4;
		double percentageSwitchesGlobal = 0.6;
		double percentageSwitches = 0.6; // [0,1]
		int stoppingCountGlobal = 200;
		int stoppingCount = 200;
		// Parameters for Genetic Tuning
		int globalN = 30;
		int localN = 20;
		int bitsGene = 12;
		double phi = 0.1;
		int numCycles = 100;
		int maxIterations = 100;
		double gammaFit = 1.0;
		double phiFit = 1.0;
		
		// Define parameters for data generation
		// DEMAND
		// Related to discrete demand
		List<List<Integer>> demandLists = new ArrayList<List<Integer>>();
		List<Integer> lowUncertaintyDemand;
		List<Integer> highUncertaintyDemand;
		lowUncertaintyDemand = new ArrayList<>(Arrays.asList(350, 400, 450));
		highUncertaintyDemand = new ArrayList<>(Arrays.asList(50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700));

		demandLists.add(lowUncertaintyDemand);
		demandLists.add(highUncertaintyDemand);		
		List<List<Double>> cdfMapsDemand = new ArrayList<List<Double>>();
		List<Double> lowUncertaintyCdf = new ArrayList<>(Arrays.asList(0.25, 0.75, 1.0));	
		List<Double> highUncertaintyCdf = new ArrayList<>(Arrays.asList(0.01, 0.08, 0.11, 0.17, 0.22, 0.31, 0.48, 0.68, 0.79, 0.90, 0.91, 0.95, 0.99, 1.0));
		cdfMapsDemand.add(lowUncertaintyCdf);
		cdfMapsDemand.add(highUncertaintyCdf);
		// Related to normal demand
		List<List<Double>> demandParameters = new ArrayList<List<Double>>();
		List<Double> lowUncertaintyDemandParameters = new ArrayList<>(Arrays.asList(400.0, 36.0));
		List<Double> highUncertaintyDemandParameters = new ArrayList<>(Arrays.asList(370.0, 142.0));
		demandParameters.add(lowUncertaintyDemandParameters);
		demandParameters.add(highUncertaintyDemandParameters);
		
		// LEAD TIME
		List<List<Double>> leadTimeParameters = new ArrayList<List<Double>>();
		List<Double> lowUncertaintyLtParameters = new ArrayList<>(Arrays.asList(1.5, 1.0));
		List<Double> highUncertaintyLtParameters = new ArrayList<>(Arrays.asList(3.0, 1.0));
		leadTimeParameters.add(lowUncertaintyLtParameters);
		leadTimeParameters.add(highUncertaintyLtParameters);
		
		// MARKET PRICES
		Map<Material,Double> meanPrices = new HashMap<Material,Double>();
		List<Map<Material,Double>> marketPriceStdvs = new ArrayList<Map<Material,Double>>();
		Map<Material,Double> lowPriceStdvs = new HashMap<Material,Double>();
		Map<Material,Double> highPriceStdvs = new HashMap<Material,Double>();
		for(Material rawMaterial : initSc.getMaterials()) {
			if(rawMaterial.getTypeCode().equals("ROH")) {
				if(rawMaterial.getMaterialCode().equals("CC-R01")) {
					meanPrices.put(rawMaterial, 0.99);
					lowPriceStdvs.put(rawMaterial, 0.099);
					highPriceStdvs.put(rawMaterial, 0.297);
				}
				if(rawMaterial.getMaterialCode().equals("CC-R02")) {
					meanPrices.put(rawMaterial, 0.86);
					lowPriceStdvs.put(rawMaterial, 0.086);
					highPriceStdvs.put(rawMaterial, 0.258);
				}
			}
		}
		marketPriceStdvs.add(lowPriceStdvs);
		marketPriceStdvs.add(highPriceStdvs);		
		
		int counter = 1;
		// Run all scenarios(2^3 scenarios)
		for(int i=0; i<2; i++) {
			for(int j=0; j<2; j++) {
				for(int k=0; k<2; k++) {
						System.out.println("SCENARIO " + counter + ":");
						System.out.println("");
						DataReader readerNew = new DataReader(splitRatio, isDemand);
						readerNew.readDataSimulation(startDate, numDates);
						SC sc = readerNew.getSupplyChain();
						// DATA GENERATION
						List<Integer> demandsList = demandLists.get(i);
						List<Double> cdfMapDemand = cdfMapsDemand.get(i);
						double meanDemand = demandParameters.get(i).get(0);
						double stdvDemand = demandParameters.get(i).get(1);
						double shapeWh = leadTimeParameters.get(j).get(0);
						double scaleWh = leadTimeParameters.get(j).get(1);
						double shapeP = leadTimeParameters.get(j).get(0);
						double scaleP = leadTimeParameters.get(j).get(1);
						Map<Material,Double> priceStdvs  = marketPriceStdvs.get(k);
						boolean isSimulation = true;
						sc.setIsSimulation(isSimulation);
						DataGeneration dataGenerator;
						if(continuousDemand) {
							// Normal demand
							dataGenerator = new DataGeneration(sc, readerNew, dates, continuousDemand, meanDemand, stdvDemand, null, null, shapeWh, scaleWh, shapeP, scaleP, meanPrices, priceStdvs);
						} else {
							// Discrete demand
							dataGenerator = new DataGeneration(sc, readerNew, dates, continuousDemand, 0.0, 0.0, demandsList, cdfMapDemand, shapeWh, scaleWh, shapeP, scaleP, meanPrices, priceStdvs);
						}
						dataGenerator.generateSCData();
						
						dataGenerator.printResults();
						
						// Store final results of EOQ model in .txt files
						Writer writer = new Writer(sc);
						writer.writeFinalResultsTraining(null, dataGenerator, 0, counter);
						writer.writeFinalResultsTest(null, dataGenerator, 0, counter);
						
						// EXACT MODEL
						Solver solverExact = new Solver(sc, numLabels, alpha, tempGlobal, percentageSwitchesGlobal, stoppingCountGlobal, globalN, bitsGene, phi, numCycles, maxIterations, gammaFit, phiFit, counter);
						// Run optimization model
						solverExact.trainFisModel(multiEchelon);					
						// RESET data in SC
						dataGenerator.addDataToSc();
						// Test model
						solverExact.testModel(false);
					
						// RESET data in SC
						dataGenerator.addDataToSc();
											
						// HEURISTIC					
						Solver solverHeuristic = new Solver(sc, numLabels, alpha, tempLocal, percentageSwitches, stoppingCount, localN, bitsGene, phi, numCycles, maxIterations, gammaFit, phiFit, counter);
						// Run optimization model
						solverHeuristic.trainFISModelHeuristic();										
						// RESET data in SC
						dataGenerator.addDataToSc();
						// Test model
						solverHeuristic.testModel(true);
						
						// Change name of input data files
						File directory = new File("C:\\Users\\iliad\\Desktop\\Eclipse_workspace\\KB_Generation");
						File[] files = directory.listFiles();
						for (File f : files)
						{
							if (f.getName().startsWith("inputData"))
							{
								File oldFile = new File(f.getName() + "Scenario " + counter + ".txt");
								f.renameTo(oldFile);
							}
						}
						// Increase scenario counter
						counter += 1;
				}
			}
		}
		
		/*
		// TEST SIMULATION RESULTS
		TestInstance test = new TestInstance(initSc, gammaFit, phiFit, demandLists, cdfMapsDemand, leadTimeParameters, demandParameters, meanPrices, marketPriceStdvs);
		test.runTest(100, numDates, startDate, splitRatio, isDemand, continuousDemand);
		
		// REAL WORLD APPLICATION
		// Read data
		double splitRatio2 = 0.75;
		boolean isDemand2 = true;
		DataReader reader2 = new DataReader(splitRatio2, isDemand2);
		reader2.readData();
		SC sc2 = reader2.getSupplyChain();
		TreeSet<Date> dates2 = sc2.getDates();
		
		
		// Define parameters
		int numLabels2 = 3;
		// Parameters for Simulated Annealing
		double alpha2 = 0.95;
		double temp2 = 0.4;
		double percentageSwitches2 = 0.6; // [0,1]
		int stoppingCount2 = 200;
		// Parameters for Genetic Tuning
		int N2 = 20;
		int bitsGene2 = 12;
		double phi2 = 0.1;
		int numCycles2 = 100;
		int maxIterations2 = 250;
		double gammaFit2 = 1.0;
		double phiFit2 = 1.0;
		int counter2 = 100;
		// DATA GENERATION
		boolean isSimulation2 = false;
		sc2.setIsSimulation(isSimulation2);
		DataGeneration dataGenerator2 = new DataGeneration(sc2, reader2, dates2);
		dataGenerator2.generateSCData();
		dataGenerator2.printResults();
		
		// Store final results of EOQ model in .txt files
		Writer writer2 = new Writer(sc2);
		writer2.writeFinalResultsTraining(null, dataGenerator2, 0, counter2);
		writer2.writeFinalResultsTest(null, dataGenerator2, 0, counter2);
		// Write ILs and IPs resulting from EOQ model into .txt file
		writer2.writeInventoryLevelsTraining(null, dataGenerator2, 0);
		writer2.writeInventoryLevelsTest(null, dataGenerator2, 0);
		writer2.writeInventoryPositionsTraining(null, dataGenerator2, 0);
		writer2.writeInventoryPositionsTest(null, dataGenerator2, 0);
			
		// Heuristic
		Solver solverHeuristic2 = new Solver(sc2, numLabels2, alpha2, temp2, percentageSwitches2, stoppingCount2, N2, bitsGene2, phi2, numCycles2, maxIterations2, gammaFit2, phiFit2, counter2);
		// Run optimization model
		solverHeuristic2.trainFISModelHeuristic();
		// Add input data to SC object
		dataGenerator2.addDataToSc();
		// Test model
		solverHeuristic2.testModel(true);
		
		// Change names of input data files
		File directory = new File("C:\\Users\\iliad\\Desktop\\Eclipse_workspace\\KB_Generation");
		File[] files = directory.listFiles();
		for (File f : files)
		{
			if (f.getName().startsWith("inputData"))
			{
				File oldFile = new File(f.getName() + "Scenario " + counter2 + ".txt");
				f.renameTo(oldFile);
			}
		}
		*/
	}	
}