package edu.scripps.yates.pcq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.pcq.compare.model.MyTTest;
import edu.scripps.yates.pcq.compare.model.QuantifiedSite;
import edu.scripps.yates.pcq.compare.model.QuantifiedSiteSet;
import edu.scripps.yates.pcq.compare.model.TTestMatrix;
import edu.scripps.yates.utilities.appversion.AppVersion;
import edu.scripps.yates.utilities.files.FileUtils;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.maths.PValueCorrection;
import edu.scripps.yates.utilities.maths.PValueCorrectionResult;
import edu.scripps.yates.utilities.maths.PValueCorrectionType;
import edu.scripps.yates.utilities.maths.PValuesCollection;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import edu.scripps.yates.utilities.sequence.PositionInPeptide;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import smile.math.Histogram;
import smile.math.Math;
import smile.netlib.NLMatrix;

public class QuantSiteOutputComparator {
	private static Logger log = Logger.getLogger(QuantSiteOutputComparator.class);
	private static Options options;
	private static File fileOfFiles;
	private static Map<File, String> sampleNamesByFiles = new THashMap<File, String>();
	private static AppVersion version;
	// default values
	private static final PValueCorrectionType defaultPValueCorrectionMethod = PValueCorrectionType.BY;
	private static final double defaultQValueThreshold = 0.05;
	private static final String currentFolder = System.getProperty("user.dir");
	private final List<File> inputFiles = new ArrayList<File>();
	private final Double rInf;
	private String outputFileName;
	private final PValueCorrectionType pValueCorrectionMethod;
	private final double qValueThreshold;
	private final int numberSigmas;
	private Double distributionSigma;
	private Double distributionAverage;
	private final int minNumberOfDiscoveries;

	public QuantSiteOutputComparator(List<File> inputFiles, double rInf, String outputFileName,
			PValueCorrectionType pValueCorrectionType, double qValueThreshold, int numberSigmas,
			int minNumberOfDiscoveries) {
		this.inputFiles.addAll(inputFiles);
		this.rInf = rInf;
		this.outputFileName = outputFileName;
		if (FilenameUtils.getExtension(this.outputFileName).equals("")) {
			this.outputFileName = this.outputFileName + ".tsv";
		}
		this.qValueThreshold = qValueThreshold;
		this.numberSigmas = numberSigmas;
		pValueCorrectionMethod = pValueCorrectionType;
		this.minNumberOfDiscoveries = minNumberOfDiscoveries;
	}

	public static void main(String[] args) {
		version = ProteinClusterQuant.getVersion();
		System.out.println("Running Quant Site comparator version " + version.toString());
		setupCommandLineOptions();
		final CommandLineParser parser = new BasicParser();
		QuantSiteOutputComparator quantSiteComparator = null;
		try {
			final CommandLine cmd = parser.parse(options, args);

			final List<File> inputFiles = new ArrayList<File>();

			fileOfFiles = new File(cmd.getOptionValue("f"));
			if (!fileOfFiles.exists()) {
				fileOfFiles = new File(currentFolder + File.separator + cmd.getOptionValue("f"));
				if (!fileOfFiles.exists()) {
					throw new Exception("Input file not found");
				}
			}
			final List<String> lines = Files.readAllLines(Paths.get(fileOfFiles.toURI()));
			for (final String line : lines) {
				final String sampleName = line.split("\t")[0].trim();
				final String fileFolder = line.split("\t")[1].trim();
				File folder = new File(fileFolder);
				if (folder.isFile() && folder.exists()) {
					inputFiles.add(folder);
					sampleNamesByFiles.put(folder, sampleName);
					continue;
				} else {
					folder = new File(currentFolder + File.separator + fileFolder);
					if (folder.isFile() && folder.exists()) {
						inputFiles.add(folder);
						sampleNamesByFiles.put(folder, sampleName);
						continue;
					}
				}
				if (!folder.exists()) {
					throw new Exception("Folder '" + folder.getAbsolutePath() + "' not found");
				}
				final File[] files = folder.listFiles(new FilenameFilter() {

					@Override
					public boolean accept(File dir, String name) {
						if (name.contains("peptideNodeTable")) {
							return true;
						}
						return false;
					}
				});
				if (files.length == 0) {
					throw new Exception("peptideNodeTable file not found at folder '" + folder.getAbsolutePath() + "'");
				}
				if (files.length > 1) {
					throw new Exception(
							"More than one file contains 'peptideNodeTable' tag in the name. It is not clear which one to take.");
				}
				final File file = files[0];
				if (!file.exists()) {
					throw new Exception(
							"File for sample '" + sampleName + "' is not found at: '" + file.getAbsolutePath() + "'");
				}
				inputFiles.add(file);
				sampleNamesByFiles.put(file, sampleName);
				log.info("File added for sample '" + sampleName + "' : " + file.getAbsolutePath() + "'");
			}

			Double rInf = null;
			try {
				rInf = Double.valueOf(cmd.getOptionValue("RInf"));
				if (rInf < 0) {
					throw new Exception(
							"Option 'RInf' must be a positive number. (Negative Infinities will be replaced by -RInf.)");
				}
			} catch (final NumberFormatException e) {
				throw new Exception("Option 'RInf' must be numerical");
			}
			final String outputFileName = cmd.getOptionValue("out");

			PValueCorrectionType pValueCorrectionType = defaultPValueCorrectionMethod;
			if (cmd.hasOption("pvc")) {
				try {
					pValueCorrectionType = PValueCorrectionType.valueOf(cmd.getOptionValue("pvc"));
					log.info("Using p-value correction method: " + pValueCorrectionType + " ("
							+ pValueCorrectionType.getReference() + ")");
				} catch (final Exception e) {
					final String errorMessage = "Invalid p-value correction method '" + cmd.getOptionValue("pvc")
							+ "'. Valid values are " + PValueCorrectionType.getValuesString();
					throw new Exception(errorMessage);
				}
			} else {

				log.info("pvc (p-value correction) parameter wasn't set. Using " + defaultPValueCorrectionMethod + " ("
						+ defaultPValueCorrectionMethod.getReference() + ") by default");

			}
			double qValueThreshold = defaultQValueThreshold;
			if (cmd.hasOption("qvt")) {
				try {
					qValueThreshold = Double.valueOf(cmd.getOptionValue("qvt"));
					if (qValueThreshold < 0 || qValueThreshold > 1) {
						throw new Exception();
					}
					log.info("Using q-value threshold = " + qValueThreshold);
				} catch (final Exception e) {
					final String errorMessage = "Invalid q-value threshold '" + cmd.getOptionValue("qvt")
							+ "'. A number between 0 and 1 are valid";
					throw new Exception(errorMessage);
				}
			} else {
				log.info("qvt (q-value threshold) parameter wasn't set. Using " + defaultQValueThreshold
						+ " by default");
			}
			int minNumberOfDiscoveries = 0; // by default
			if (cmd.hasOption("md")) {
				try {
					minNumberOfDiscoveries = Integer.valueOf(cmd.getOptionValue("md"));
					if (minNumberOfDiscoveries < 0) {
						throw new Exception();
					}
					log.info("Using minimum_discoveries = " + minNumberOfDiscoveries);
				} catch (final Exception e) {
					final String errorMessage = "Invalid md value '" + cmd.getOptionValue("md")
							+ "'. A positive number greater or equal to 0 is valid";
					throw new Exception(errorMessage);
				}
			} else {
				log.info("md (minimum_discoveries) parameter wasn't set. Using " + minNumberOfDiscoveries
						+ " by default. However, only sites with at least one discovery will be reported in the Excel output file.");
			}

			int numberSigmas = 2; // by default
			if (cmd.hasOption("ns")) {
				try {
					numberSigmas = Integer.valueOf(cmd.getOptionValue("ns"));
					if (numberSigmas < 0) {
						throw new Exception();
					}
					log.info("Using number_sigmas = " + numberSigmas);
				} catch (final Exception e) {
					final String errorMessage = "Invalid ns value '" + cmd.getOptionValue("ns")
							+ "'. A positive number greater or equal to 0 is valid";
					throw new Exception(errorMessage);
				}
			} else {
				log.info("ns (number_sigmas) parameter wasn't set. Using " + numberSigmas + " by default.");
			}
			quantSiteComparator = new QuantSiteOutputComparator(inputFiles, rInf, outputFileName, pValueCorrectionType,
					qValueThreshold, numberSigmas, minNumberOfDiscoveries);

		} catch (final Exception e) {
			e.printStackTrace();
			errorInParameters();
		}
		if (quantSiteComparator != null) {
			try {
				quantSiteComparator.run();
			} catch (final IOException e) {
				e.printStackTrace();
				log.error("Program exited unexpectedly.");
				System.exit(-1);
			}
		} else {
			log.error("Program exited unexpectedly.");
			System.exit(-1);
		}
		log.info("Program exited correctly.");
		System.exit(0);

	}

	private void run() throws IOException {
		QuantifiedSiteSet quantSites = null;
		for (final File file : inputFiles) {
			final String sampleName = getSampleNameByFile(file);
			if (quantSites == null) {
				quantSites = readFile(file);
				quantSites.addSampleName(sampleName);
			} else {
				final QuantifiedSiteSet quantSites2 = readFile(file);
				final List<String> sampleNames = quantSites.getSampleNames();
				quantSites = mergeQuantifiedSiteSets(quantSites, quantSites2);
				for (final String sampleName2 : sampleNames) {
					quantSites.addSampleName(sampleName2);
				}
				quantSites.addSampleName(sampleName);
			}

		}
		writeTripletsOutput(quantSites);
		writePairWisePValueMatrixes(quantSites);
	}

	private static String getSampleNameByFile(File file) {
		return sampleNamesByFiles.get(file);
	}

	private void writePairWisePValueMatrixes(QuantifiedSiteSet quantSites) throws IOException {
		final int numSamples = quantSites.getNumExperiments();

		final Map<String, TTestMatrix> matrixMap = new THashMap<String, TTestMatrix>();
		for (final QuantifiedSite quantSite : quantSites.getSortedByRatios()) {
			final TTestMatrix matrix = new TTestMatrix(numSamples, numSamples);
			matrixMap.put(quantSite.getNodeKey(), matrix);
			for (int sampleIndex1 = 0; sampleIndex1 < numSamples; sampleIndex1++) {
				for (int sampleIndex2 = sampleIndex1 + 1; sampleIndex2 < numSamples; sampleIndex2++) {
					if (quantSite.getNodeKey().equals("O75369#K1838") && sampleIndex1 == 7 && sampleIndex2 == 9) {
						log.info(quantSite);
					}
					final MyTTest ttest = performTTest(quantSite, sampleIndex1, sampleIndex2,
							getDistributionAverage(quantSites), getDistributionSigma(quantSites));
					matrix.set(sampleIndex1, sampleIndex2, ttest);
				}
			}
		}
		// correct the pValues per pair of samples, that is
		// sampleIndex1+sampleIndex2
		final ProgressCounter counter = new ProgressCounter(
				Double.valueOf(numSamples * (numSamples - 1) / 2).intValue(), ProgressPrintingType.PERCENTAGE_STEPS, 1,
				true);
		// to keep the number of significant pvalues after the pvalue correction
		final TObjectIntHashMap<String> numberOfDiscoveriesPerSite = new TObjectIntHashMap<String>();
		for (int sampleIndex1 = 0; sampleIndex1 < numSamples; sampleIndex1++) {
			for (int sampleIndex2 = sampleIndex1 + 1; sampleIndex2 < numSamples; sampleIndex2++) {
				final TObjectDoubleHashMap<String> pValues = new TObjectDoubleHashMap<String>();
				final Set<String> quantSiteKeys = quantSites.getQuantifiedSitesByKey().keySet();
				for (final String quantSite : quantSiteKeys) {
					final TTestMatrix matrix = matrixMap.get(quantSite);
					if (matrix != null) {
						final MyTTest myTTest = matrix.get(sampleIndex1, sampleIndex2);
						if (myTTest.isUseForPVAlueCorrection()) {
							final double pValue = myTTest.getPValue();
							if (!Double.isNaN(pValue)) {
								pValues.put(quantSite, pValue);
							}
						}
					}
				}
				final PValuesCollection pValueCollection = new PValuesCollection(pValues);
				final PValueCorrectionResult pAdjust = PValueCorrection.pAdjust(pValueCollection,
						pValueCorrectionMethod);
				for (final String quantSite : quantSiteKeys) {
					final TTestMatrix matrix = matrixMap.get(quantSite);
					if (matrix != null) {
						final Double adjustedPValue = pAdjust.getCorrectedPValues().getPValue(quantSite);
						if (adjustedPValue != null) {
							if (adjustedPValue < qValueThreshold) {
								if (numberOfDiscoveriesPerSite.contains(quantSite)) {
									numberOfDiscoveriesPerSite.put(quantSite,
											numberOfDiscoveriesPerSite.get(quantSite) + 1);
								} else {
									numberOfDiscoveriesPerSite.put(quantSite, 1);
								}
							}
							matrix.get(sampleIndex1, sampleIndex2).setPValue(adjustedPValue);
						}
					}
					if (numberOfDiscoveriesPerSite.get(quantSite) < minNumberOfDiscoveries) {
						matrixMap.remove(quantSite);
					}
				}
				counter.increment();
				final String printIfNecessary = counter.printIfNecessary();
				if (!"".equals(printIfNecessary)) {
					log.info("Correcting pValues... " + printIfNecessary);
				}
			}
		}
		if (matrixMap.size() < 1) {
			log.warn("There is not any matrix with a minimum number of discoveries of " + minNumberOfDiscoveries
					+ " !!!");
		} else {
			log.info("Writting output files with " + matrixMap.size() + " matrixes ...");
		}
		final File matrixSummaryFile = getMatrixSummaryFile();
		final File excelSummaryFile = getExcelMatrixSummaryFile();
		// delete excel file if exists
		if (excelSummaryFile.exists()) {
			excelSummaryFile.delete();
		}
		final TIntArrayList nums = new TIntArrayList();

		// delete previous matrixes
		final File matrixFolder = getIndividualMatrixFolder();
		if (!matrixFolder.exists()) {
			matrixFolder.mkdirs();
		}
		if (matrixFolder.listFiles().length > 0) {
			for (final File file : matrixFolder.listFiles()) {
				file.delete();
			}
		}
		// this integer matrix will store how many time each pair of samples is
		// significant
		final NLMatrix sampleComparisonMatrix = new NLMatrix(numSamples, numSamples, 0);

		final FileWriter fw = new FileWriter(matrixSummaryFile);
		boolean atLeastOneMatrix = false;
		// sort by number of discoveries and by ratios
		for (final QuantifiedSite quantifiedSite : quantSites
				.getSortedByNumDiscoveriesAndProteinsAndSites(numberOfDiscoveriesPerSite)) {
			fw.write(quantifiedSite.getNodeKey() + "\n");
			fw.write("Number of discoveries:\t" + numberOfDiscoveriesPerSite.get(quantifiedSite.getNodeKey()) + "\n");

			final TTestMatrix matrix = matrixMap.get(quantifiedSite.getNodeKey());
			if (matrix == null) {
				continue;
			}
			atLeastOneMatrix = true;
			fw.write(printMatrix(matrix, quantifiedSite) + "\n");

			// independent file with the matrix
			final int numDiscoveries = numberOfDiscoveriesPerSite.get(quantifiedSite.getNodeKey());
			nums.add(numDiscoveries);

			final File individualMatrixFile = getIndividualMatrixFile(numDiscoveries, quantifiedSite.getNodeKey());
			if (numDiscoveries > 0 && numDiscoveries >= minNumberOfDiscoveries) {
				final FileWriter individualMatrixFileWriter = new FileWriter(individualMatrixFile);
				individualMatrixFileWriter.write(printMatrix(matrix, quantifiedSite));
				individualMatrixFileWriter.close();
				for (int i = 0; i < numSamples; i++) {
					for (int j = i + 1; j < numSamples; j++) {
						final double pvalue = matrix.get(i, j).getPValue();
						if (pvalue < qValueThreshold) {
							sampleComparisonMatrix.set(i, j, sampleComparisonMatrix.get(i, j) + 1);
						}
					}
				}

				// add the individual file to the Excel file

				FileUtils.separatedValuesToXLSX(individualMatrixFile.getAbsolutePath(),
						excelSummaryFile.getAbsolutePath(), "\t", numDiscoveries + "_" + quantifiedSite.getNodeKey());

			}
		}
		final String message = "Comparison matrix: number sites in which each pair of sample has been found significant with q-value < "
				+ qValueThreshold;
		System.out.println(message);
		fw.write("\n\n" + message + "\n");
		System.out.println(printMatrix(sampleComparisonMatrix, null));
		fw.write(printMatrix(sampleComparisonMatrix, null) + "\n");
		fw.close();
		if (atLeastOneMatrix) {
			log.info("Output file written at: '" + matrixSummaryFile.getAbsolutePath() + "'");
			log.info("Output file written at: '" + excelSummaryFile.getAbsolutePath() + "'");
		} else {
			log.info(
					"No output was generated because there is not any quantified site with a minimum number of discoveries (significant pairwise comparisons with q-value<"
							+ qValueThreshold + ") of " + minNumberOfDiscoveries);
		}
		if (!nums.isEmpty() && nums.max() > 1) {
			final double[][] histogram = Histogram.histogram(nums.toArray(), nums.max());
			int numSitesWithAtLeastOneDiscovery = 0;
			for (int j = 0; j < histogram[2].length; j++) {
				System.out.println("Number of sites with " + j + " discoveries " + histogram[2][j]);
				if (j > 0) {
					numSitesWithAtLeastOneDiscovery += histogram[2][j];
				}
			}
			System.out.println("Number of sites with at least one discovery: " + numSitesWithAtLeastOneDiscovery);
		}
	}

	private File getIndividualMatrixFile(int numDiscoveries, String site) {
		final File file = new File(getIndividualMatrixFolder().getAbsoluteFile().getAbsolutePath() + File.separator
				+ numDiscoveries + "_" + site + ".txt");
		if (!file.getParentFile().exists()) {
			file.getParentFile().mkdirs();
		}
		return file;
	}

	private File getIndividualMatrixFolder() {
		return new File(getOutputFolder() + File.separator + "matrixes_" + FilenameUtils.getBaseName(outputFileName));
	}

	private File getMatrixSummaryFile() {
		return new File(getOutputFolder() + File.separator + FilenameUtils.getBaseName(outputFileName)
				+ "_QVALUES_matrixes.txt");
	}

	private File getExcelMatrixSummaryFile() {
		final File file = new File(getOutputFolder() + File.separator + FilenameUtils.getBaseName(outputFileName)
				+ "_QVALUES_matrixes.xlsx");

		return file;
	}

	private String getOutputFolder() {
		if (fileOfFiles != null) {
			return fileOfFiles.getParentFile().getAbsolutePath();
		} else {
			return inputFiles.get(0).getParentFile().getAbsolutePath();
		}
	}

	private void writeTripletsOutput(QuantifiedSiteSet mergedQuantSites) throws IOException {

		// get it sorted
		final List<QuantifiedSite> sortedQuantifiedSites = mergedQuantSites.getSortedByRatios();

		log.info("Writting output file...");
		final File outputFile = new File(getOutputFolder() + File.separator + outputFileName);
		final FileWriter fw = new FileWriter(outputFile);
		// header
		fw.write(QuantifiedSite.NODE_KEY + "\t");
		for (int index = 0; index < mergedQuantSites.getNumExperiments(); index++) {
			fw.write(QuantifiedSite.LOG2RATIO + "_" + (index + 1) + "\t");
			fw.write(QuantifiedSite.RATIOSCOREVALUE + "_" + (index + 1) + "\t");
			// fw.write(QuantifiedSite.NUMPSMS + "_" + (index + 1) + "\t");
			fw.write(QuantifiedSite.NUMMEASUREMENTS + "_" + (index + 1) + "\t");
		}
		fw.write(QuantifiedSite.SEQUENCE + "\t");
		fw.write(QuantifiedSite.POSITIONSINPEPTIDE + "\t");
		fw.write(QuantifiedSite.PROTEINS + "\t");
		fw.write(QuantifiedSite.GENES + "\t");
		fw.write("p-value" + "\t");
		fw.write("q-value (corrected by " + pValueCorrectionMethod + ")" + "\t");
		fw.write("\n");
		for (final QuantifiedSite quantSite : sortedQuantifiedSites) {
			final QuantifiedSite quantifiedSite = mergedQuantSites.getQuantifiedSitesByKey()
					.get(quantSite.getNodeKey());
			fw.write(quantifiedSite.getNodeKey() + "\t");
			for (int index = 0; index < mergedQuantSites.getNumExperiments(); index++) {
				final Double log2Ratio = replaceInfinite(quantifiedSite.getLog2Ratio(index), rInf);
				fw.write(log2Ratio + "\t");
				final int numMeasurements = quantifiedSite.getNumMeasurements(index);
				Double stdev = quantifiedSite.getRatioScoreValue(index);
				if (numMeasurements <= 1) {
					stdev = Double.NaN;
				}

				fw.write(stdev + "\t");
				fw.write(numMeasurements + "\t");
			}

			fw.write(quantifiedSite.getSequence() + "\t");
			fw.write(getSeparatedValueStringFromChars(quantifiedSite.getPositionsInPeptide(), "-") + "\t");
			fw.write(quantifiedSite.getProteins() + "\t");
			fw.write(quantifiedSite.getGenes() + "\t");

			fw.write("\n");
		}
		fw.close();
		log.info("Output file written at: '" + outputFile.getAbsolutePath() + "'");

	}

	private double replaceInfinite(Double log2Ratio, Double rInf) {

		if (log2Ratio == null || Double.isNaN(log2Ratio)) {
			return Double.NaN;
		}
		if (rInf == null) {
			return log2Ratio;
		}
		if (Double.isInfinite(log2Ratio)) {
			if (Double.POSITIVE_INFINITY == log2Ratio) {
				return rInf;
			} else if (Double.NEGATIVE_INFINITY == log2Ratio) {
				return -rInf;
			}
		}
		return log2Ratio;
	}

	private MyTTest performTTest(QuantifiedSite quantSite, int sampleIndex1, int sampleIndex2, double distributionAvg,
			double distributionSigma) {
		final double mean1 = quantSite.getLog2Ratio(sampleIndex1);
		final double mean2 = quantSite.getLog2Ratio(sampleIndex2);
		final double stdev1 = quantSite.getRatioScoreValue(sampleIndex1);
		final double var1 = Math.pow(stdev1, 2);
		final double stdev2 = quantSite.getRatioScoreValue(sampleIndex2);
		final double var2 = Math.pow(stdev2, 2);
		final int n1 = quantSite.getNumMeasurements(sampleIndex1);
		final int n2 = quantSite.getNumMeasurements(sampleIndex2);

		if (n1 < 2 || n2 < 2 || stdev1 == Double.NaN || stdev2 == Double.NaN) {
			return new MyTTest(Double.NaN, false);
		}
		if ((Double.isInfinite(mean1) && Double.isInfinite(mean2))) {
			if (mean1 == mean2) {
				return new MyTTest(1.0, false); // no significant
			} else {
				return new MyTTest(0.0, false); // significant
			}
		} else {
			if (Double.isInfinite(mean1) || Double.isInfinite(mean2)) {
				// only significant if the one that is FINITE is beyond n sigmas
				// of the distribution
				final double finite = Double.isFinite(mean1) ? mean1 : mean2;
				final double infinite = Double.isInfinite(mean1) ? mean1 : mean2;
				if (Double.POSITIVE_INFINITY == infinite) {
					if (finite < distributionAvg + numberSigmas * distributionSigma) {
						return new MyTTest(0.0, false);
					} else {
						return new MyTTest(1.0, false);
					}
				} else {
					// negative infinity
					if (finite > distributionAvg + numberSigmas * distributionSigma) {
						// significant
						return new MyTTest(0.0, false);
					} else {
						// no significant
						return new MyTTest(1.0, false);
					}
				}
			}
		}
		final edu.scripps.yates.utilities.maths.TTest pValue = edu.scripps.yates.utilities.maths.TTest.test(mean1, var1,
				n1, mean2, var2, n2, true);
		return new MyTTest(pValue, true);

	}

	private double getDistributionSigma(QuantifiedSiteSet quantSites) {
		if (distributionSigma == null) {
			log.info("Calculating the standard deviation of the whole distribution of ratios...");
			final TDoubleArrayList ratios = new TDoubleArrayList();
			for (final QuantifiedSite quantifiedSite : quantSites) {
				for (int experimentIndex = 0; experimentIndex < quantifiedSite.getNumExperiments(); experimentIndex++) {
					final Double log2Ratio = quantifiedSite.getLog2Ratio(experimentIndex);
					if (log2Ratio != null && !Double.isNaN(log2Ratio) && !Double.isInfinite(log2Ratio)) {
						ratios.add(log2Ratio);
					}
				}
			}
			distributionSigma = Maths.stddev(ratios.toArray());
			log.info("The standard deviation of the whole distribution of ratios in all experiments is: "
					+ distributionSigma);
		}
		return distributionSigma;
	}

	private double getDistributionAverage(QuantifiedSiteSet quantSites) {
		if (distributionAverage == null) {
			log.info("Calculating the mean of the whole distribution of ratios...");
			final TDoubleArrayList ratios = new TDoubleArrayList();
			for (final QuantifiedSite quantifiedSite : quantSites) {
				for (int experimentIndex = 0; experimentIndex < quantifiedSite.getNumExperiments(); experimentIndex++) {
					final Double log2Ratio = quantifiedSite.getLog2Ratio(experimentIndex);
					if (log2Ratio != null && !Double.isNaN(log2Ratio) && !Double.isInfinite(log2Ratio)) {
						ratios.add(log2Ratio);
					}
				}
			}
			distributionAverage = Maths.mean(ratios.toArray());
			log.info("The mean of the whole distribution of ratios in all experiments is: " + distributionAverage);
		}
		return distributionAverage;
	}

	public static String getSeparatedValueStringFromChars(Collection<PositionInPeptide> collection, String separator) {
		final StringBuilder sb = new StringBuilder();
		for (final Object c : collection) {
			if (!"".equals(sb.toString())) {
				sb.append(separator);
			}
			sb.append(c.toString());
		}
		return sb.toString();
	}

	private QuantifiedSiteSet mergeQuantifiedSiteSets(QuantifiedSiteSet quantSites1, QuantifiedSiteSet quantSites2) {
		log.info("Merging two quantified site sets");
		final QuantifiedSiteSet ret = new QuantifiedSiteSet();
		final Set<String> keySet1 = quantSites1.getQuantifiedSitesByKey().keySet();
		for (final String keySet : keySet1) {
			if (quantSites2.getQuantifiedSitesByKey().contains(keySet)) {
				final QuantifiedSite quantSite1 = quantSites1.getQuantifiedSitesByKey().get(keySet);
				final QuantifiedSite quantSite2 = quantSites2.getQuantifiedSitesByKey().get(keySet);
				if (Double.isNaN(quantSite1.getLog2Ratio(0))) {
					if (Double.isNaN(quantSite2.getLog2Ratio(0))) {
						continue;
					}
				}
				// use the first for the result, adding the values of the second
				quantSite1.addLog2Ratio(quantSite2.getLog2Ratio(0));
				quantSite1.addNumMeasurements(quantSite2.getNumMeasurements(0));
				quantSite1.addNumPeptides(quantSite2.getNumPeptides(0));
				quantSite1.addNumPSMs(quantSite2.getNumPSMs(0));
				quantSite1.addRatioScoreValue(quantSite2.getRatioScoreValue(0));
				quantSite1.getPositionsInPeptide().addAll(quantSite2.getPositionsInPeptide());
				ret.add(quantSite1);
			}
		}
		log.info("Result of merging is " + ret.getQuantifiedSitesByKey().size() + " quantified sites from "
				+ quantSites1.getQuantifiedSitesByKey().size() + " + " + quantSites2.getQuantifiedSitesByKey().size());
		log.info("Now having " + ret.getNumExperiments() + " experiments");
		return ret;
	}

	private QuantifiedSiteSet readFile(File inputFile) throws IOException {
		int numLine = 1;
		String line = null;
		try {
			log.info("Reading input file: '" + inputFile.getAbsolutePath() + "'");
			final BufferedReader br = new BufferedReader(new FileReader(inputFile));
			final QuantifiedSiteSet ret = new QuantifiedSiteSet();

			final TObjectIntHashMap<String> indexesByHeaders = new TObjectIntHashMap<String>();
			while ((line = br.readLine()) != null) {
				try {
					final String[] split = line.split("\t");
					if (numLine == 1) {
						for (int i = 0; i < split.length; i++) {
							indexesByHeaders.put(split[i], i);
						}
					} else {
						final QuantifiedSite quantSite = new QuantifiedSite(split, indexesByHeaders);
						ret.add(quantSite);
					}
				} finally {
					numLine++;
				}
			}
			br.close();
			log.info(ret.getQuantifiedSitesByKey().size() + " quantified sites read from input file '"
					+ inputFile.getAbsolutePath() + "'");
			return ret;
		} catch (final Exception e) {
			log.error("LINE  " + numLine + ": " + line);
			log.error("Error reading at line " + numLine + " of file '" + inputFile.getAbsolutePath() + "'");
			throw e;
		}
	}

	private static void setupCommandLineOptions() {
		// create Options object
		options = new Options();
		final Option opt1 = new Option("RInf", "replace_infinity", true,
				"[OPTIONAL] -RInf replaces +/- Infinity with a user defined (+/-) value in the output summary table file");
		opt1.setRequired(true);
		options.addOption(opt1);
		final Option opt2 = new Option("f", "input_file", true,
				"[MANDATORY] Full path to a file containing pairs (separated by TAB) of sample names and full path to the peptideNodeTable of a PCQ run to compare");
		opt2.setRequired(true);
		options.addOption(opt2);
		final Option opt3 = new Option("out", "output_file_name", true,
				"[MANDATORY] Output file name that will be created in the current folder");
		opt3.setRequired(true);
		options.addOption(opt3);
		final Option opt4 = new Option("pvc", "pvalue_correction", true,
				"[OPTIONAL] p-value correction method to apply. Valid values are: "
						+ PValueCorrectionType.getValuesString() + ". If not provided, the method will be "
						+ defaultPValueCorrectionMethod + " (Reference: " + defaultPValueCorrectionMethod.getReference()
						+ ")");
		opt4.setRequired(false);
		options.addOption(opt4);
		final Option opt5 = new Option("qvt", "qvalue_threshold", true,
				"[OPTIONAL] q-value threshold to apply to the corrected p-values. A value between 0 and 1 is permitted. If not provided, a threshold of "
						+ defaultQValueThreshold + " will be applied.");
		opt5.setRequired(false);

		final Option opt6 = new Option("md", "minimum_discoveries", true,
				"[OPTIONAL] minimum number of discoveries (significantly different between two samples) required for a quantified site to be in the output files. If not provided, there will be no minimum number, although no quant sites without any significantly different site between 2 samples will be reported in the Excel output file.");
		opt6.setRequired(false);
		options.addOption(opt6);

		final Option opt7 = new Option("ns", "number_sigmas", true,
				"[OPTIONAL] number of sigmas that will be used to decide whether an INFINITY ratio is significantly different than a FINITE ratio.\n"
						+ "If R1=POSITIVE_INFINITY and R2 < avg_distribution + ns*sigma_distribution_of_ratios, then R2 is significantly different.\n"
						+ "If R1=NEGATIVE_INFINITY and R2 > avg_distribution + ns*sigma_distribution_of_ratios, then R2 is significantly different.");
		opt7.setRequired(false);
		options.addOption(opt7);
	}

	private static void errorInParameters() {
		// automatically generate the help statement
		final HelpFormatter formatter = new HelpFormatter();

		formatter.printHelp(150, "java -jar QuantSiteoutputComparator.jar", "with the following parameters:", options,
				"\n\nContact Salvador Martinez-Bartolome at salvador@scripps.edu for more help");

		System.exit(0);
	}

	/**
	 * Returns the string representation of matrix.
	 * 
	 * @param quantifiedSite
	 * 
	 * @param full
	 *            Print the full matrix if true. Otherwise only print top left 7
	 *            x 7 submatrix.
	 */
	public String printMatrix(TTestMatrix matrix, QuantifiedSite quantifiedSite) {
		final StringBuilder sb = new StringBuilder();
		final int numRows = matrix.nrows();
		final int numCols = matrix.ncols();

		final String newline = numCols < matrix.ncols() ? "...\n" : "\n";
		// header
		if (quantifiedSite != null) {
			sb.append("Site(s)\tPosition(s)\tProtein(s)\tGene(s)\n");
			sb.append(quantifiedSite.getNodeKey() + "\t" + quantifiedSite.getPositions() + "\t"
					+ quantifiedSite.getProteins() + "\t" + quantifiedSite.getGenes() + "\n");
		}
		sb.append("\t");
		if (quantifiedSite != null) {
			sb.append("\t");
		}
		for (int i = 0; i < numCols; i++) {
			sb.append(getSampleNameByFile(inputFiles.get(i)) + "\t");
		}
		sb.append("\n");
		sb.append("\t");
		if (quantifiedSite != null) {
			sb.append("\t");
			for (int i = 0; i < numCols; i++) {
				final Double log2Ratio = quantifiedSite.getLog2Ratio(i);
				sb.append(log2Ratio + "\t");
			}
		}
		sb.append("\n");
		for (int i = 0; i < numRows; i++) {
			sb.append(getSampleNameByFile(inputFiles.get(i)) + "\t");
			if (quantifiedSite != null) {
				final Double log2Ratio = quantifiedSite.getLog2Ratio(i);
				sb.append(log2Ratio + "\t");
			}

			for (int j = 0; j < numCols; j++) {
				sb.append(matrix.get(i, j) + "\t");
			}
			sb.append(newline);
		}

		if (numRows < matrix.nrows()) {
			sb.append("  ...\n");
		}

		return sb.toString();
	}

	/**
	 * Returns the string representation of matrix.
	 * 
	 * @param quantifiedSite
	 * 
	 * @param full
	 *            Print the full matrix if true. Otherwise only print top left 7
	 *            x 7 submatrix.
	 */
	public String printMatrix(NLMatrix matrix, QuantifiedSite quantifiedSite) {
		final StringBuilder sb = new StringBuilder();
		final int numRows = matrix.nrows();
		final int numCols = matrix.ncols();

		final String newline = numCols < matrix.ncols() ? "...\n" : "\n";
		// header
		if (quantifiedSite != null) {
			sb.append("Site(s)\tPosition(s)\tProtein(s)\tGene(s)\n");
			sb.append(quantifiedSite.getNodeKey() + "\t" + quantifiedSite.getPositions() + "\t"
					+ quantifiedSite.getProteins() + "\t" + quantifiedSite.getGenes() + "\n");
		}
		sb.append("\t");
		if (quantifiedSite != null) {
			sb.append("\t");
		}
		for (int i = 0; i < numCols; i++) {
			sb.append(getSampleNameByFile(inputFiles.get(i)) + "\t");
		}
		sb.append("\n");
		sb.append("\t");
		if (quantifiedSite != null) {
			sb.append("\t");
			for (int i = 0; i < numCols; i++) {
				final Double log2Ratio = quantifiedSite.getLog2Ratio(i);
				sb.append(log2Ratio + "\t");
			}
		}
		sb.append("\n");
		for (int i = 0; i < numRows; i++) {
			sb.append(getSampleNameByFile(inputFiles.get(i)) + "\t");
			if (quantifiedSite != null) {
				final Double log2Ratio = quantifiedSite.getLog2Ratio(i);
				sb.append(log2Ratio + "\t");
			}

			for (int j = 0; j < numCols; j++) {
				sb.append(matrix.get(i, j) + "\t");
			}
			sb.append(newline);
		}

		if (numRows < matrix.nrows()) {
			sb.append("  ...\n");
		}

		return sb.toString();
	}
}
