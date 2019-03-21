package edu.scripps.yates.pcq.quantsite;

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
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.log4j.Logger;

import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.pcq.ProteinClusterQuant;
import edu.scripps.yates.pcq.compare.model.MyTTest;
import edu.scripps.yates.pcq.compare.model.QuantifiedSite;
import edu.scripps.yates.pcq.compare.model.QuantifiedSiteSet;
import edu.scripps.yates.pcq.compare.model.TTestMatrix;
import edu.scripps.yates.pcq.quantsite.tmt.TMTPairWisePCQInputParametersGenerator;
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
import gnu.trove.set.hash.THashSet;
import smile.math.Histogram;
import smile.math.Math;
import smile.netlib.NLMatrix;

public class QuantSiteOutputComparator {
	private static Logger log = Logger.getLogger(QuantSiteOutputComparator.class);
	private static Options options;
	private static File fileOfFiles;
	public static Map<File, String> sampleNamesByFiles = new THashMap<File, String>();
	private static AppVersion version;
	// default values
	public static final PValueCorrectionType defaultPValueCorrectionMethod = PValueCorrectionType.BY;
	public static final double defaultQValueThreshold = 0.05;
	private static final String SWAP = "SWAP";

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
	private THashMap<QuantifiedSite, TTestMatrix> ttestMatrixesByQuantSites;
	private TObjectIntHashMap<QuantifiedSite> numberOfDiscoveriesPerSite;
	private final boolean tmtData;
	private String outputFolder;
	private THashMap<QuantifiedSite, NLMatrix> matrixByQuantSitesForTMT;
	private final Set<File> ratioSwaps = new THashSet<File>();

	public QuantSiteOutputComparator(List<File> inputFiles, Set<File> ratioSwap, double rInf, String outputFileName,
			PValueCorrectionType pValueCorrectionType, double qValueThreshold, int numberSigmas,
			int minNumberOfDiscoveries, boolean tmtData) {
		this.inputFiles.addAll(inputFiles);
		if (ratioSwap != null) {
			this.ratioSwaps.addAll(ratioSwap);
		}
		this.rInf = rInf;
		this.outputFileName = outputFileName;
		if (FilenameUtils.getExtension(this.outputFileName).equals("")) {
			this.outputFileName = this.outputFileName + ".tsv";
		}
		this.qValueThreshold = qValueThreshold;
		this.numberSigmas = numberSigmas;
		pValueCorrectionMethod = pValueCorrectionType;
		this.minNumberOfDiscoveries = minNumberOfDiscoveries;
		this.tmtData = tmtData;
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
			final Set<File> swapFiles = new THashSet<File>();
			fileOfFiles = new File(cmd.getOptionValue("f"));
			if (!fileOfFiles.exists()) {
				fileOfFiles = new File(currentFolder + File.separator + cmd.getOptionValue("f"));
				if (!fileOfFiles.exists()) {
					throw new Exception("Input file not found");
				}
			}
			final List<String> lines = Files.readAllLines(Paths.get(fileOfFiles.toURI()));
			for (final String line : lines) {
				final String[] split = line.split("\t");
				final String sampleName = split[0].trim();
				final String fileFolder = split[1].trim();

				File folder = new File(fileFolder);
				if (folder.isFile() && folder.exists()) {
					inputFiles.add(folder);
					sampleNamesByFiles.put(folder, sampleName);
					if (split.length == 3) {
						if (split[2].trim().equals(SWAP)) {
							swapFiles.add(folder);
						}
					}
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

				final File[] files = folder.listFiles(getPeptideNodeFileFilter());
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
				if (split.length == 3) {
					if (split[2].trim().equals(SWAP)) {
						swapFiles.add(file);
					}
				}
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
			// tmt
			boolean tmtData = false;
			if (cmd.hasOption("tmt")) {
				try {
					tmtData = Boolean.valueOf(cmd.getOptionValue("tmt"));
				} catch (final NumberFormatException e) {
					throw new Exception(
							"Option 'tmt' must be boolean value (true or false). If not provided it will be 'false'.");
				}
			}

			quantSiteComparator = new QuantSiteOutputComparator(inputFiles, swapFiles, rInf, outputFileName,
					pValueCorrectionType, qValueThreshold, numberSigmas, minNumberOfDiscoveries, tmtData);

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

	public static FilenameFilter getPeptideNodeFileFilter() {
		final FilenameFilter fileFilter = new FilenameFilter() {

			@Override
			public boolean accept(File dir, String name) {
				if (name.contains("peptideNodeTable")) {
					return true;
				}
				return false;
			}
		};
		return fileFilter;
	}

	public void run() throws IOException {
		if (tmtData) {
			runTMTData();
		} else {
			runRegularData();
		}
	}

	public void runTMTData() throws IOException {
		QuantifiedSiteSet quantSites = null;
		// here each file is a comparison between 2 samples
		for (final File file : inputFiles) {
			final String samplePairName = getSampleNameByFile(file);
			if (quantSites == null) {
				quantSites = readPCQOutputFile(file, this.ratioSwaps.contains(file));
				quantSites.addSampleName(samplePairName);
			} else {
				final QuantifiedSiteSet quantSites2 = readPCQOutputFile(file, this.ratioSwaps.contains(file));
				final List<String> sampleNames = quantSites.getSampleNames();
				quantSites = mergeQuantifiedSiteSets(quantSites, quantSites2);
				for (final String sampleName2 : sampleNames) {
					quantSites.addSampleName(sampleName2);
				}
				quantSites.addSampleName(samplePairName);
			}

		}

		writePairWisePValueMatrixesForTMT(quantSites);
		writeRatioTableOutput(quantSites);
	}

	public void runRegularData() throws IOException {
		// first, get the total set of quantified sites
		final Set<String> keySets = new THashSet<String>();
		final Set<String> nonNanSites = new THashSet<String>();
		final Map<File, QuantifiedSiteSet> quantSiteSetsByFile = new THashMap<File, QuantifiedSiteSet>();
		for (final File file : inputFiles) {
			final QuantifiedSiteSet quantSitesTMP = readPCQOutputFile(file, ratioSwaps.contains(file));
			quantSiteSetsByFile.put(file, quantSitesTMP);
			keySets.addAll(quantSitesTMP.getQuantifiedSitesByKey().keySet());
			for (final QuantifiedSite quantSite : quantSitesTMP.getSortedByRatios()) {
				if (!Double.isNaN(quantSite.getLog2Ratio(0))) {
					nonNanSites.add(quantSite.getNodeKey());
				}
			}
		}
		final QuantifiedSiteSet quantSites = new QuantifiedSiteSet();
		log.info(nonNanSites.size() + "/" + keySets.size() + " quantified sites in " + inputFiles.size()
				+ " experiments");
		for (final String siteKey : keySets) {

			for (final File file : inputFiles) {
				boolean swapRatios = false;
				if (this.ratioSwaps.contains(file)) {
					swapRatios = true;
				}
				final String sampleName = getSampleNameByFile(file);
				quantSites.addSampleName(sampleName);

				final QuantifiedSite quantifiedSiteTMP = quantSiteSetsByFile.get(file).getQuantifiedSitesByKey()
						.get(siteKey);
				final QuantifiedSite quantifiedSite = quantSites.getQuantifiedSitesByKey().get(siteKey);

				if (quantifiedSite != null && quantifiedSiteTMP != null) {
					// use the first for the result, adding the values of the second
					final Double log2Ratio = quantifiedSiteTMP.getLog2Ratio(0);

					quantifiedSite.addLog2Ratio(log2Ratio);
					quantifiedSite.addNumMeasurements(quantifiedSiteTMP.getNumMeasurements(0));
					quantifiedSite.addNumPeptides(quantifiedSiteTMP.getNumPeptides(0));
					quantifiedSite.addNumPSMs(quantifiedSiteTMP.getNumPSMs(0));
					quantifiedSite.addRatioStdev(quantifiedSiteTMP.getRatioStdev(0));
					quantifiedSite.getPositionsInPeptide().addAll(quantifiedSiteTMP.getPositionsInPeptide());
					if (quantifiedSite.getProteins() == null) {
						quantifiedSite.setProteins(quantifiedSiteTMP.getProteins());
						quantifiedSite.setGenes(quantifiedSiteTMP.getGenes());
						quantifiedSite.setPositionsInPeptide(quantifiedSiteTMP.getPositionsInPeptide());
						quantifiedSite.setSequence(quantifiedSiteTMP.getSequence());
					}
				} else if (quantifiedSite != null && quantifiedSiteTMP == null) {
					quantifiedSite.addLog2Ratio(Double.NaN);
					quantifiedSite.addNumMeasurements(0);
					quantifiedSite.addNumPeptides(0);
					quantifiedSite.addNumPSMs(0);
					quantifiedSite.addRatioStdev(Double.NaN);
				} else if (quantifiedSite == null && quantifiedSiteTMP != null) {
					quantSites.add(quantifiedSiteTMP);
				} else if (quantifiedSite == null && quantifiedSiteTMP == null) {
					quantSites.add(new QuantifiedSite(siteKey));
				}

			}
		}
		calculateNumberOfDiscoveriesPerSite(quantSites);
		writeRatioTableOutput(quantSites);
		writePairWisePValueMatrixes(quantSites);
	}

	private static String getSampleNameByFile(File file) {
		return sampleNamesByFiles.get(file);
	}

	private void calculateNumberOfDiscoveriesPerSite(QuantifiedSiteSet quantSiteSet) {
		final int numSamples = quantSiteSet.getNumExperiments();

		ttestMatrixesByQuantSites = new THashMap<QuantifiedSite, TTestMatrix>();
		for (final QuantifiedSite quantSite : quantSiteSet.getSortedByRatios()) {
			final TTestMatrix matrix = new TTestMatrix(numSamples, numSamples);
			ttestMatrixesByQuantSites.put(quantSite, matrix);
			for (int sampleIndex1 = 0; sampleIndex1 < numSamples; sampleIndex1++) {
				for (int sampleIndex2 = sampleIndex1 + 1; sampleIndex2 < numSamples; sampleIndex2++) {
					final MyTTest ttest = performTTest(quantSite, sampleIndex1, sampleIndex2,
							getDistributionAverage(quantSiteSet), getDistributionSigma(quantSiteSet));
					matrix.set(sampleIndex1, sampleIndex2, ttest);
				}
			}
		}
		// correct the pValues per pair of samples, that is
		// sampleIndex1+sampleIndex2

		// to keep the number of significant pvalues after the pvalue correction
		numberOfDiscoveriesPerSite = new TObjectIntHashMap<QuantifiedSite>();
		int numSitesWithMinimumDiscoveries = 0;
		for (int sampleIndex1 = 0; sampleIndex1 < numSamples; sampleIndex1++) {
			for (int sampleIndex2 = sampleIndex1 + 1; sampleIndex2 < numSamples; sampleIndex2++) {
				final TObjectDoubleHashMap<QuantifiedSite> pValuesByQuantSite = new TObjectDoubleHashMap<QuantifiedSite>();
				final Collection<QuantifiedSite> quantSites = quantSiteSet.getQuantifiedSitesByKey().values();
				for (final QuantifiedSite quantSite : quantSites) {
					final TTestMatrix matrix = ttestMatrixesByQuantSites.get(quantSite);
					if (matrix != null) {
						final MyTTest myTTest = matrix.get(sampleIndex1, sampleIndex2);
						if (myTTest.isUseForPValueCorrection()) {
							final double pValue = myTTest.getPValue();
							if (!Double.isNaN(pValue)) {
								pValuesByQuantSite.put(quantSite, pValue);
							}
						}
					}
				}
				final PValuesCollection<QuantifiedSite> pValueCollection = new PValuesCollection<QuantifiedSite>(
						pValuesByQuantSite);
				PValueCorrectionResult<QuantifiedSite> pAdjust = null;
				if (pValuesByQuantSite.size() > 0) {
					log.info("Sample '" + getSampleNameByFile(inputFiles.get(sampleIndex1)) + "'"
							+ TMTPairWisePCQInputParametersGenerator.VS + "'"
							+ getSampleNameByFile(inputFiles.get(sampleIndex2)) + "':");
					log.info("Adjusting " + pValuesByQuantSite.size() + " p-values using method '"
							+ pValueCorrectionMethod.name() + "'");
					pAdjust = PValueCorrection.pAdjust(pValueCollection, pValueCorrectionMethod);
				} else {
					pAdjust = new PValueCorrectionResult<QuantifiedSite>();
					pAdjust.setCorrectedPValues(pValueCollection);
					pAdjust.setOriginalPValues(pValueCollection);
				}
				for (final QuantifiedSite quantSite : quantSites) {
					final TTestMatrix matrixOfTTests = ttestMatrixesByQuantSites.get(quantSite);
					if (matrixOfTTests != null) {
						Double adjustedPValue = pAdjust.getCorrectedPValues().getPValue(quantSite);
						if (adjustedPValue == null) {
							// here he have to take the sites that are
							// significant because we compared INFINITIES, but
							// in this case, they are not in the corrected
							// pvalues
							adjustedPValue = matrixOfTTests.get(sampleIndex1, sampleIndex2).getPValue();
						}
						if (adjustedPValue != null) {
							matrixOfTTests.get(sampleIndex1, sampleIndex2).setCorrectedPValue(adjustedPValue);
							if (adjustedPValue < qValueThreshold) {
								if (numberOfDiscoveriesPerSite.contains(quantSite)) {
									numberOfDiscoveriesPerSite.put(quantSite,
											numberOfDiscoveriesPerSite.get(quantSite) + 1);
								} else {
									numberOfDiscoveriesPerSite.put(quantSite, 1);
								}
							}
						}
					}
					if (numberOfDiscoveriesPerSite.get(quantSite) >= minNumberOfDiscoveries) {
						numSitesWithMinimumDiscoveries++;
					}
				}

			}
		}
		if (numSitesWithMinimumDiscoveries == 0) {
			log.warn("There is not any matrix with a minimum number of discoveries of " + minNumberOfDiscoveries
					+ " !!!");
		} else {
			log.info("Writting output files with " + numSitesWithMinimumDiscoveries
					+ " quant sites with minimum number of discoveries (" + minNumberOfDiscoveries + ")...");
		}
	}

	/**
	 * Then, for each quantified site, I build a matrix 10x10, where I store a
	 * t-test result comparing the ratios of that site in each pairwise sample
	 * comparison.
	 * 
	 * @param quantSiteSet
	 * @throws IOException
	 */
	private void writePairWisePValueMatrixes(QuantifiedSiteSet quantSiteSet) throws IOException {
		final int numSamples = quantSiteSet.getNumExperiments();

		final File matrixSummaryFile = getMatrixSummaryFile();
		final File excelSummaryFile = getExcelMatrixSummaryFile();
		// delete excel file if exists
		if (excelSummaryFile.exists()) {
			excelSummaryFile.delete();
		}
		final TIntArrayList nums = new TIntArrayList();

		// delete previous matrixes in the folder
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
		for (final QuantifiedSite quantifiedSite : quantSiteSet
				.getSortedByNumDiscoveriesAndProteinsAndSites(numberOfDiscoveriesPerSite)) {
			if (numberOfDiscoveriesPerSite.get(quantifiedSite) >= minNumberOfDiscoveries) {
				fw.write(quantifiedSite.getNodeKey() + "\n");
				fw.write("Number of discoveries:\t" + numberOfDiscoveriesPerSite.get(quantifiedSite) + "\n");

				final TTestMatrix matrixOfTTests = ttestMatrixesByQuantSites.get(quantifiedSite);
				if (matrixOfTTests == null) {
					continue;
				}

				atLeastOneMatrix = true;
				fw.write(printMatrix(matrixOfTTests, quantifiedSite, true) + "\n");

				// keep number of discoveries
				final int numDiscoveries = numberOfDiscoveriesPerSite.get(quantifiedSite);
				nums.add(numDiscoveries);

				// independent file with the matrix
				final File individualMatrixFile = File.createTempFile("node", "tmp");
				individualMatrixFile.deleteOnExit();
				if (numDiscoveries > 0 && numDiscoveries >= minNumberOfDiscoveries) {
					final FileWriter individualMatrixFileWriter = new FileWriter(individualMatrixFile);
					individualMatrixFileWriter.write(printMatrix(matrixOfTTests, quantifiedSite, true));
					individualMatrixFileWriter.close();

					for (int i = 0; i < numSamples; i++) {
						for (int j = i + 1; j < numSamples; j++) {
							final double pvalue = matrixOfTTests.get(i, j).getCorrectedPValue();
							if (pvalue < qValueThreshold) {
								sampleComparisonMatrix.set(i, j, sampleComparisonMatrix.get(i, j) + 1);
							}
						}
					}

					// add the individual file to the Excel file

//					FileUtils.separatedValuesToXLSX(individualMatrixFile.getAbsolutePath(),
//							excelSummaryFile.getAbsolutePath(), "\t",
//							numDiscoveries + "_" + quantifiedSite.getNodeKey());

				}
			}
		}
		final String message = "Comparison matrix: number sites in which each pair of sample has been found significant with q-value < "
				+ qValueThreshold;
		System.out.println(message);
		fw.write("\n\n" + message + "\n");
		System.out.println(printMatrix(sampleComparisonMatrix, null));
		fw.write(printMatrix(sampleComparisonMatrix, null) + "\n");

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
				final String m = "Number of sites with " + j + " discoveries " + histogram[2][j];
				System.out.println(m);
				fw.write(m + "\n");
				if (j > 0) {
					numSitesWithAtLeastOneDiscovery += histogram[2][j];
				}
			}
			final String l = "Number of sites with at least one discovery: " + numSitesWithAtLeastOneDiscovery;
			System.out.println(l);
			fw.write(l + "\n");
		}
		fw.close();
	}

	private int getNumSamplesFromSamplePairs(int numSamplesPairs) {
		return Double.valueOf((1 + Math.sqrt(1 + 8 * numSamplesPairs)) / 2).intValue();
	}

	private void writePairWisePValueMatrixesForTMT(QuantifiedSiteSet quantSites) throws IOException {
		final int numSamplesPairs = quantSites.getNumExperiments();
		final int numSamples = getNumSamplesFromSamplePairs(numSamplesPairs);
		// create a matrix per quantSite
		matrixByQuantSitesForTMT = new THashMap<QuantifiedSite, NLMatrix>();
		for (final QuantifiedSite quantSite : quantSites.getSortedByRatios()) {
			final NLMatrix matrix = new NLMatrix(numSamples, numSamples, Double.NaN);
			matrixByQuantSitesForTMT.put(quantSite, matrix);
		}
		// to keep the number of significant pvalues after the pvalue
		// correction
		numberOfDiscoveriesPerSite = new TObjectIntHashMap<QuantifiedSite>();

		//
		for (int sampleIndexPair = 0; sampleIndexPair < numSamplesPairs; sampleIndexPair++) {
			final String samplePairName = quantSites.getSampleNames().get(sampleIndexPair);

			// for each sample pair (each sampleIndex), get all ratios to build
			// the ratio distribution
			final TDoubleArrayList ratioCollection = new TDoubleArrayList();
			for (final QuantifiedSite quantSite : quantSites.getSortedByRatios()) {
				final Double ratio = quantSite.getLog2Ratio(sampleIndexPair);
				if (ratio != null && !Double.isNaN(ratio) && !Double.isInfinite(ratio)) {
					ratioCollection.add(ratio);
				}
			}
			// Now, we normalize the distribution to be centered at 0 (mean 0)
			final double mean = Maths.mean(ratioCollection);
			if (Double.compare(mean, 0.0) != 0) {
				log.info("Mean of all sites ratios for comparison " + samplePairName + " is " + mean
						+ ". Correcting distribution to center at 0.");
				int offset = 0;
				for (final double ratio : ratioCollection.toArray()) {
					ratioCollection.setQuick(offset++, ratio - mean);
				}
				// check that now the mean is 0
				final double newMean = Maths.mean(ratioCollection);
				log.info("After correction, the mean of all sites ratios for comparison " + samplePairName + " is "
						+ newMean);
			}
			log.info("Final ratio distribution for " + samplePairName + ": mean=" + Maths.mean(ratioCollection)
					+ " and stdev=" + Maths.stddev(ratioCollection));
			// create the gaussian fitter
			// make an histogram
			log.info("Creating histogram of distribution of ratios");
			final double[][] histogram = Histogram.histogram(ratioCollection.toArray());
			final GaussianCurveFitter gaussianfitter = GaussianCurveFitter.create();
			final WeightedObservedPoints obs = new WeightedObservedPoints();
			for (int j = 0; j < histogram[2].length; j++) {
				obs.add((histogram[0][j] + histogram[1][j]) / 2, histogram[2][j]);
			}
			final double[] bestFit = gaussianfitter.fit(obs.toList());
			// final double gaussianNorm = bestFit[0];
			final double gaussianMean = bestFit[1];
			final double gaussianSigma = bestFit[2];
			final NormalDistribution normalDistribution = new NormalDistribution(gaussianMean, gaussianSigma);
			// now, per each site, calculate the p-value against the
			// distribution
			// keep pvalues here:
			final Set<QuantifiedSite> significantSitesBecauseOfInfinity = new THashSet<QuantifiedSite>();

			final TObjectDoubleHashMap<QuantifiedSite> pvalues = new TObjectDoubleHashMap<QuantifiedSite>();
			for (final QuantifiedSite quantSite : quantSites.getSortedByRatios()) {
				double pvalue = 1.0;
				final Double ratio = quantSite.getLog2Ratio(sampleIndexPair) - mean;// shifting
																					// the
																					// ratio

				if (ratio != null && !Double.isNaN(ratio) && !Double.isInfinite(ratio)) {
					pvalue = normalDistribution.cumulativeProbability(ratio);
					// this p is the probability that a ratio is <= than ratio
					if (ratio > gaussianMean) {
						pvalue = 1 - pvalue; // this would be the probability
												// that a ratio is > than the
												// ratio
					}
					pvalues.put(quantSite, pvalue);
				} else if (Double.isInfinite(ratio)) {
					significantSitesBecauseOfInfinity.add(quantSite);
				}

			}
			// create the pvalues collection
			final PValuesCollection<QuantifiedSite> pValuesCollection = new PValuesCollection<QuantifiedSite>(pvalues);
			// correct the pvalues
			final PValueCorrectionResult<QuantifiedSite> adjustedPValues = PValueCorrection.pAdjust(pValuesCollection,
					pValueCorrectionMethod);

			for (final QuantifiedSite quantSite : quantSites.getQuantifiedSitesByKey().values()) {
				Double correctedPValue = adjustedPValues.getCorrectedPValues().getPValue(quantSite);
				if (correctedPValue == null) {
					// it is because is coming from an infinity and therefore is
					// not in the corrected pvalues because it wasn't use for
					// the pvalue correction. Anyway, the corrected pvalue is 0
					if (significantSitesBecauseOfInfinity.contains(quantSite)) {
						correctedPValue = 0.0;
					}
				}
				final int rowSampleIndex = getRowSampleIndex(samplePairName);
				final int columnSampleIndex = getColumnSampleIndex(samplePairName);
				if (correctedPValue == null) {
					// this is because the original non corrected pvalues was
					// NaN
					correctedPValue = Double.NaN;
				}
				matrixByQuantSitesForTMT.get(quantSite).set(rowSampleIndex, columnSampleIndex, correctedPValue);
				if (correctedPValue < qValueThreshold) {
					if (numberOfDiscoveriesPerSite.contains(quantSite)) {
						numberOfDiscoveriesPerSite.put(quantSite, numberOfDiscoveriesPerSite.get(quantSite) + 1);
					} else {
						numberOfDiscoveriesPerSite.put(quantSite, 1);
					}
				}
			}
		}

		// if (numberOfDiscoveriesPerSite.get(quantSite) >=
		// minNumberOfDiscoveries) {
		// numSitesWithMinimumDiscoveries++;
		// }
		//
		// if (numSitesWithMinimumDiscoveries == 0) {
		// log.warn("There is not any matrix with a minimum number of
		// discoveries of " + minNumberOfDiscoveries
		// + " !!!");
		// } else {
		// log.info("Writting output files with " +
		// numSitesWithMinimumDiscoveries
		// + " quant sites with minimum number of discoveries...");
		// }
		final File matrixSummaryFile = getMatrixSummaryFile();
		final File excelSummaryFile = getExcelMatrixSummaryFile();
		// delete excel file if exists
		if (excelSummaryFile.exists()) {
			excelSummaryFile.delete();
		}
		final TIntArrayList nums = new TIntArrayList();

		// delete previous matrixes in the folder
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
		final int total = quantSites.getQuantifiedSitesByKey().size();
		final ProgressCounter counter = new ProgressCounter(total, ProgressPrintingType.PERCENTAGE_STEPS, 0, true);
		counter.setSuffix(" writting Excel file...");
		for (final QuantifiedSite quantifiedSite : quantSites
				.getSortedByNumDiscoveriesAndProteinsAndSites(numberOfDiscoveriesPerSite)) {
			counter.increment();
			if (!"".equals(counter.printIfNecessary())) {
				log.info(counter.printIfNecessary());
			}
			fw.write(quantifiedSite.getNodeKey() + "\n");
			fw.write("Number of discoveries:\t" + numberOfDiscoveriesPerSite.get(quantifiedSite) + "\n");

			final NLMatrix matrixOfPValues = matrixByQuantSitesForTMT.get(quantifiedSite);
			if (matrixOfPValues == null) {
				continue;
			}
			if (numberOfDiscoveriesPerSite.get(quantifiedSite) >= minNumberOfDiscoveries) {
				atLeastOneMatrix = true;
				fw.write(printMatrix(matrixOfPValues, quantifiedSite) + "\n");

				// keep number of discoveries
				final int numDiscoveries = numberOfDiscoveriesPerSite.get(quantifiedSite);
				nums.add(numDiscoveries);

				// independent file with the matrix
				final File individualMatrixFile = File.createTempFile("node", "tmp");
				individualMatrixFile.deleteOnExit();
				if (numDiscoveries > 0 && numDiscoveries >= minNumberOfDiscoveries) {
					final FileWriter individualMatrixFileWriter = new FileWriter(individualMatrixFile);
					individualMatrixFileWriter.write(printMatrix(matrixOfPValues, quantifiedSite));
					individualMatrixFileWriter.close();

					for (int i = 0; i < numSamples; i++) {
						for (int j = i + 1; j < numSamples; j++) {
							final double pvalue = matrixOfPValues.get(i, j);
							if (pvalue < qValueThreshold) {
								sampleComparisonMatrix.set(i, j, sampleComparisonMatrix.get(i, j) + 1);
							}
						}
					}

					// add the individual file to the Excel file

					FileUtils.separatedValuesToXLSX(individualMatrixFile.getAbsolutePath(),
							excelSummaryFile.getAbsolutePath(), "\t",
							numDiscoveries + "_" + quantifiedSite.getNodeKey());

				}
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

	private int getColumnSampleIndex(String samplePairName) {
		if (samplePairName.contains(TMTPairWisePCQInputParametersGenerator.VS)) {
			final String labelString = samplePairName.split(TMTPairWisePCQInputParametersGenerator.VS)[1];
			final QuantificationLabel label = QuantificationLabel.valueOf(labelString);
			if (QuantificationLabel.isTMT10PLEX(label)) {
				return QuantificationLabel.getTMT10PlexLabels().indexOf(label);
			} else if (QuantificationLabel.isTMT6PLEX(label)) {
				return QuantificationLabel.getTMT6PlexLabels().indexOf(label);
			}
		}
		return -1;
	}

	private int getRowSampleIndex(String samplePairName) {
		if (samplePairName.contains(TMTPairWisePCQInputParametersGenerator.VS)) {
			final String labelString = samplePairName.split(TMTPairWisePCQInputParametersGenerator.VS)[0];
			final QuantificationLabel label = QuantificationLabel.valueOf(labelString);
			if (QuantificationLabel.isTMT10PLEX(label)) {
				return QuantificationLabel.getTMT10PlexLabels().indexOf(label);
			} else if (QuantificationLabel.isTMT6PLEX(label)) {
				return QuantificationLabel.getTMT6PlexLabels().indexOf(label);
			}
		}
		return -1;
	}

	// private File getIndividualMatrixFile(int numDiscoveries, String site) {
	// final File file = new
	// File(getIndividualMatrixFolder().getAbsoluteFile().getAbsolutePath() +
	// File.separator
	// + numDiscoveries + "_" + site + ".txt");
	// if (!file.getParentFile().exists()) {
	// file.getParentFile().mkdirs();
	// }
	// return file;
	// }

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
		if (outputFolder == null) {
			if (fileOfFiles != null) {
				outputFolder = fileOfFiles.getParentFile().getAbsolutePath();
			} else {
				outputFolder = inputFiles.get(0).getParentFile().getAbsolutePath();
			}
		}
		return outputFolder;
	}

	private void writeRatioTableOutput(QuantifiedSiteSet mergedQuantSites) throws IOException {

		log.info("Writting output file...");
		final File outputFile = new File(getOutputFolder() + File.separator + outputFileName);
		final FileWriter fw = new FileWriter(outputFile, false);
		// header
		fw.write(QuantifiedSite.NODE_KEY + "\t");
		final int numSamplesFromSamplePairs = getNumSamplesFromSamplePairs(mergedQuantSites.getNumExperiments());
		if (tmtData) {
			for (int i = 0; i < numSamplesFromSamplePairs; i++) {
				for (int j = i + 1; j < numSamplesFromSamplePairs; j++) {

					fw.write(QuantifiedSite.LOG2RATIO + "_" + (i + 1) + TMTPairWisePCQInputParametersGenerator.VS
							+ (j + 1) + "\t");
					fw.write(QuantifiedSite.STDEV + "_" + (i + 1) + TMTPairWisePCQInputParametersGenerator.VS + (j + 1)
							+ "\t");
					// fw.write(QuantifiedSite.NUMPSMS + "_" + (index + 1) +
					// "\t");
					fw.write(QuantifiedSite.NUMMEASUREMENTS + "_" + (i + 1) + TMTPairWisePCQInputParametersGenerator.VS
							+ (j + 1) + "\t");
				}
			}
		} else {
			for (int index = 0; index < mergedQuantSites.getNumExperiments(); index++) {
				final int num = index + 1;
				fw.write(QuantifiedSite.LOG2RATIO + "_" + num + "\t");
				fw.write(QuantifiedSite.STDEV + "_" + num + "\t");
				// fw.write(QuantifiedSite.NUMPSMS + "_" + (index + 1) + "\t");
				fw.write(QuantifiedSite.NUMMEASUREMENTS + "_" + num + "\t");
			}
		}
		fw.write(QuantifiedSite.SEQUENCE + "\t");
		fw.write(QuantifiedSite.QUANTPOSITIONSINPEPTIDE + "\t");
		fw.write(QuantifiedSite.PROTEINS + "\t");
		fw.write(QuantifiedSite.GENES + "\t");
		fw.write("# discoveries (q-value < " + qValueThreshold + ")\t");
		if (tmtData) {
			for (int i = 0; i < numSamplesFromSamplePairs; i++) {
				for (int j = i + 1; j < numSamplesFromSamplePairs; j++) {
					fw.write((i + 1) + TMTPairWisePCQInputParametersGenerator.VS + (j + 1) + " p-value" + "\t");
				}
			}
		} else {
			for (int i = 0; i < mergedQuantSites.getNumExperiments(); i++) {
				for (int j = i + 1; j < mergedQuantSites.getNumExperiments(); j++) {
					fw.write((i + 1) + TMTPairWisePCQInputParametersGenerator.VS + (j + 1) + " p-value" + "\t");
					fw.write((i + 1) + TMTPairWisePCQInputParametersGenerator.VS + (j + 1) + " q-value (by "
							+ pValueCorrectionMethod + ")" + "\t");
				}
			}
		}

		fw.write("\n");
		for (final QuantifiedSite quantifiedSite : mergedQuantSites
				.getSortedByNumDiscoveriesAndProteinsAndSites(numberOfDiscoveriesPerSite)) {

			fw.write(quantifiedSite.getNodeKey() + "\t");
			if (tmtData) {
				for (int i = 0; i < numSamplesFromSamplePairs; i++) {
					for (int j = i + 1; j < numSamplesFromSamplePairs; j++) {
						final int index = getPairIndexFromSamplesIndexex(i, j);
						final Double log2Ratio = replaceInfiniteWithRInfParameter(quantifiedSite.getLog2Ratio(index),
								rInf);
						fw.write(log2Ratio + "\t");
						final int numMeasurements = quantifiedSite.getNumMeasurements(index);
						Double stdev = quantifiedSite.getRatioStdev(index);
						if (numMeasurements <= 1) {
							stdev = Double.NaN;
						}
						fw.write(stdev + "\t");
						fw.write(numMeasurements + "\t");
					}
				}
			} else {
				for (int index = 0; index < mergedQuantSites.getNumExperiments(); index++) {
					final Double log2Ratio = replaceInfiniteWithRInfParameter(quantifiedSite.getLog2Ratio(index), rInf);
					fw.write(log2Ratio + "\t");
					final int numMeasurements = quantifiedSite.getNumMeasurements(index);
					Double stdev = quantifiedSite.getRatioStdev(index);
					if (numMeasurements <= 1) {
						stdev = Double.NaN;
					}

					fw.write(stdev + "\t");
					fw.write(numMeasurements + "\t");
				}
			}

			fw.write(quantifiedSite.getSequence() + "\t");
			fw.write(getSeparatedValueStringFromChars(quantifiedSite.getPositionsInPeptide(), "-") + "\t");
			fw.write(quantifiedSite.getProteins() + "\t");
			fw.write(quantifiedSite.getGenes() + "\t");
			fw.write(numberOfDiscoveriesPerSite.get(quantifiedSite) + "\t");

			if (tmtData) {
				final NLMatrix matrix = matrixByQuantSitesForTMT.get(quantifiedSite);
				for (int i = 0; i < matrix.nrows(); i++) {
					for (int j = i + 1; j < matrix.ncols(); j++) {
						if (matrix != null) {
							final double pValue = matrix.get(i, j);

							fw.write(pValue + "\t");
						}
					}
				}
			} else {
				final TTestMatrix tTestMatrix = ttestMatrixesByQuantSites.get(quantifiedSite);
				for (int i = 0; i < mergedQuantSites.getNumExperiments(); i++) {
					for (int j = i + 1; j < mergedQuantSites.getNumExperiments(); j++) {
						if (tTestMatrix != null && tTestMatrix.get(i, j) != null) {
							final MyTTest tTest = tTestMatrix.get(i, j);
							fw.write(tTest.getPValue() + "\t");
							fw.write(tTest.getCorrectedPValue() + "\t");
						} else {
							fw.write("\t\t");
						}

					}
				}
			}
			fw.write("\n");
		}
		fw.close();
		log.info("Output file written at: '" + outputFile.getAbsolutePath() + "'");

	}

	private int getPairIndexFromSamplesIndexex(int i, int j) {
		final String samplePairName = getSampleNameByFile(inputFiles.get(0));

		final String labelString = samplePairName.split(TMTPairWisePCQInputParametersGenerator.VS)[1];
		final QuantificationLabel label = QuantificationLabel.valueOf(labelString);

		int max = -1;
		if (QuantificationLabel.isTMT10PLEX(label)) {
			max = QuantificationLabel.getTMT10PlexLabels().size();
		} else if (QuantificationLabel.isTMT6PLEX(label)) {
			max = QuantificationLabel.getTMT6PlexLabels().size();
		}
		int ret = 0;
		for (int i2 = 0; i2 < max; i2++) {
			for (int j2 = i2 + 1; j2 < max; j2++) {
				if (i2 == i && j2 == j) {
					return ret;
				}
				ret++;
			}
		}
		return -1;
	}

	private double replaceInfiniteWithRInfParameter(Double log2Ratio, Double rInf) {

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
		final double stdev1 = quantSite.getRatioStdev(sampleIndex1);
		final double var1 = Math.pow(stdev1, 2);
		final double stdev2 = quantSite.getRatioStdev(sampleIndex2);
		final double var2 = Math.pow(stdev2, 2);
		final int n1 = quantSite.getNumMeasurements(sampleIndex1);
		final int n2 = quantSite.getNumMeasurements(sampleIndex2);

		if (n1 < 2 || n2 < 2 || Double.isNaN(stdev1) || Double.isNaN(stdev2)) {
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
		final Set<String> keySets = new THashSet<String>();
		keySets.addAll(quantSites1.getQuantifiedSitesByKey().keySet());
		keySets.addAll(quantSites2.getQuantifiedSitesByKey().keySet());
		for (final String keySet : keySets) {
			if (quantSites2.getQuantifiedSitesByKey().contains(keySet)) {
//				log.info("This site was already there");
			} else {
//				log.info("This site is new");
			}
			final QuantifiedSite quantSite1 = quantSites1.getQuantifiedSitesByKey().get(keySet);
			final QuantifiedSite quantSite2 = quantSites2.getQuantifiedSitesByKey().get(keySet);
			Double ratio1 = Double.NaN;
			if (quantSite1 != null) {
				ratio1 = quantSite1.getLog2Ratio(0);
			}
			Double ratio2 = Double.NaN;
			if (quantSite2 != null) {
				ratio2 = quantSite2.getLog2Ratio(0);
			}
			if (Double.isNaN(ratio1) && Double.isNaN(ratio2)) {
//				continue;
			}
			if (quantSite1 != null && quantSite2 != null) {
				// use the first for the result, adding the values of the second
				quantSite1.addLog2Ratio(quantSite2.getLog2Ratio(0));
				quantSite1.addNumMeasurements(quantSite2.getNumMeasurements(0));
				quantSite1.addNumPeptides(quantSite2.getNumPeptides(0));
				quantSite1.addNumPSMs(quantSite2.getNumPSMs(0));
				quantSite1.addRatioStdev(quantSite2.getRatioStdev(0));
				quantSite1.getPositionsInPeptide().addAll(quantSite2.getPositionsInPeptide());
				ret.add(quantSite1);
				continue;
			}
			if (quantSite1 != null && quantSite2 == null) {
				ret.add(quantSite1);
				continue;
			}
			if (quantSite2 != null && quantSite1 == null) {
				ret.add(quantSite2);
				continue;
			}
//			} else {
//
//			}
		}
		log.info("Result of merging is " + ret.getQuantifiedSitesByKey().size() + " quantified sites from "
				+ quantSites1.getQuantifiedSitesByKey().size() + " + " + quantSites2.getQuantifiedSitesByKey().size());
		log.info("Now having " + ret.getNumExperiments() + " experiments");
		return ret;
	}

	/**
	 * Reads PCQ output file
	 * 
	 * @param b
	 */
	private QuantifiedSiteSet readPCQOutputFile(File inputFile, boolean swapRatio) throws IOException {
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
						if (swapRatio) {
							final double swappedRatio = 1 / quantSite.getLog2Ratio(0);
							quantSite.setLog2Ratio(0, swappedRatio);
						}
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
		opt1.setRequired(false);
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
		options.addOption(opt5);
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
		final Option opt8 = new Option("tmt", "tmt_data", true,
				"[OPTIONAL] whether the data is TMT or not. In case of having TMT data, it will perform the analysis accordingly. Bu default this parameter is FALSE.");
		opt8.setRequired(false);
		options.addOption(opt8);
	}

	private static void errorInParameters() {
		// automatically generate the help statement
		final HelpFormatter formatter = new HelpFormatter();

		formatter.printHelp(150, "java -jar QuantSiteoutputComparator.jar", "with the following parameters:", options,
				"\n\nContact Salvador Martinez-Bartolome at salvador@scripps.edu for more help");

	}

	/**
	 * Returns the string representation of matrix.
	 * 
	 * @param quantifiedSite
	 * 
	 * @param full           Print the full matrix if true. Otherwise only print top
	 *                       left 7 x 7 submatrix.
	 */
	private String printMatrix(TTestMatrix matrix, QuantifiedSite quantifiedSite, boolean printCorrectedPValues) {
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
				if (Double.isInfinite(log2Ratio)) {
					sb.append("'");
				}
				sb.append(log2Ratio + "\t");
			}
		}
		sb.append("\n");
		for (int i = 0; i < numRows; i++) {
			sb.append(getSampleNameByFile(inputFiles.get(i)) + "\t");
			if (quantifiedSite != null) {
				final Double log2Ratio = quantifiedSite.getLog2Ratio(i);
				if (Double.isInfinite(log2Ratio)) {
					sb.append("'");
				}
				sb.append(log2Ratio + "\t");
			}

			for (int j = 0; j < numCols; j++) {
				final MyTTest ttest = matrix.get(i, j);
				if (ttest != null) {
					if (printCorrectedPValues) {
						sb.append(ttest.getCorrectedPValue() + "\t");
					} else {
						sb.append(ttest.getPValue() + "\t");
					}
				} else {
					sb.append("-\t");
				}
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
	 * @param full           Print the full matrix if true. Otherwise only print top
	 *                       left 7 x 7 submatrix.
	 */
	private String printMatrix(NLMatrix matrix, QuantifiedSite quantifiedSite) {
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
			if (tmtData) {
				sb.append(getTMTSampleName(getSampleNameByFile(inputFiles.get(i)), i) + "\t");

			} else {
				sb.append(getSampleNameByFile(inputFiles.get(i)) + "\t");
			}
		}
		sb.append("\n");
		sb.append("\t");
		if (quantifiedSite != null) {
			sb.append("\t");
			for (int i = 0; i < numCols; i++) {
				if (tmtData) {
					sb.append("-" + "\t");
				} else {
					final Double log2Ratio = quantifiedSite.getLog2Ratio(i);
					if (Double.isInfinite(log2Ratio)) {
						sb.append("'");
					}
					sb.append(log2Ratio + "\t");
				}
			}
		}
		sb.append("\n");
		for (int i = 0; i < numRows; i++) {
			if (tmtData) {
				sb.append(getTMTSampleName(getSampleNameByFile(inputFiles.get(i)), i) + "\t");
			} else {
				sb.append(getSampleNameByFile(inputFiles.get(i)) + "\t");
			}
			if (quantifiedSite != null) {
				if (tmtData) {
					sb.append("-" + "\t");
				} else {
					final Double log2Ratio = quantifiedSite.getLog2Ratio(i);
					if (Double.isInfinite(log2Ratio)) {
						sb.append("'");
					}
					sb.append(log2Ratio + "\t");
				}
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

	private String getTMTSampleName(String samplePairName, int tmtIndex) {
		if (samplePairName.contains(TMTPairWisePCQInputParametersGenerator.VS)) {
			final String labelString = samplePairName.split(TMTPairWisePCQInputParametersGenerator.VS)[1];
			final QuantificationLabel label = QuantificationLabel.valueOf(labelString);
			if (QuantificationLabel.isTMT10PLEX(label)) {
				return QuantificationLabel.getTMT10PlexLabels().get(tmtIndex).name();
			} else if (QuantificationLabel.isTMT6PLEX(label)) {
				return QuantificationLabel.getTMT6PlexLabels().get(tmtIndex).name();
			}
		}
		return null;
	}

	public void setOutputFolder(String outputFolder) {
		this.outputFolder = outputFolder;
	}
}
