package edu.scripps.yates.pcq.quantsite.tmt;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.log4j.Logger;

import edu.scripps.yates.pcq.PCQBatchRunner;
import edu.scripps.yates.pcq.ProteinClusterQuant;
import edu.scripps.yates.pcq.quantsite.QuantSiteOutputComparator;
import edu.scripps.yates.utilities.appversion.AppVersion;
import edu.scripps.yates.utilities.maths.PValueCorrectionType;

public class TMTQuantSiteComparator {
	private final static Logger log = Logger.getLogger(TMTQuantSiteComparator.class);
	private static Options options;
	private static AppVersion version;
	private static final boolean defaultCytParameterValue = false;
	private static final boolean defaultIgnoreNotFoundSignals = true;
	private final File paramFile;
	private final List<File> tmtFiles;
	private final String tmtType;
	private final double rInf;
	private final String outputFileName;
	private final PValueCorrectionType pValueCorrectionType;
	private final double qValueThreshold;
	private final int numberSigmas;
	private final int minNumberOfDiscoveries;
	private final boolean generateXGMMLFiles;
	private final boolean ignoreNotFoundSignals;

	public static void main(String[] args) {
		version = ProteinClusterQuant.getVersion();
		System.out.println("Running TMTPairWisePCQInputParametersGenerator version " + version.toString());
		setupCommandLineOptions();
		final CommandLineParser parser = new BasicParser();
		try {
			final CommandLine cmd = parser.parse(options, args);

			final File paramFile = new File(cmd.getOptionValue("pf"));
			final File inputFiles = new File(cmd.getOptionValue("ifs"));
			String tmtType = TMTPairWisePCQInputParametersGenerator.TMT10PLEX; // by
																				// default
			if (cmd.hasOption("tmt")) {
				tmtType = cmd.getOptionValue("tmt");
			}
			if (!TMTPairWisePCQInputParametersGenerator.TMT10PLEX.equals(tmtType)
					&& !TMTPairWisePCQInputParametersGenerator.TMT6PLEX.equals(tmtType)) {
				throw new IllegalArgumentException("Invalid value for tmt parameter: '" + tmtType
						+ "'. Valid values are " + TMTPairWisePCQInputParametersGenerator.TMT10PLEX + " or "
						+ TMTPairWisePCQInputParametersGenerator.TMT6PLEX + " in tmt parameter");
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

			PValueCorrectionType pValueCorrectionType = QuantSiteOutputComparator.defaultPValueCorrectionMethod;
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

				log.info("pvc (p-value correction) parameter wasn't set. Using "
						+ QuantSiteOutputComparator.defaultPValueCorrectionMethod + " ("
						+ QuantSiteOutputComparator.defaultPValueCorrectionMethod.getReference() + ") by default");

			}
			double qValueThreshold = QuantSiteOutputComparator.defaultQValueThreshold;
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
				log.info("qvt (q-value threshold) parameter wasn't set. Using "
						+ QuantSiteOutputComparator.defaultQValueThreshold + " by default");
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

			boolean generateXGMMLFiles = defaultCytParameterValue; // by default
			if (cmd.hasOption("cyt")) {
				try {
					generateXGMMLFiles = Boolean.valueOf(cmd.getOptionValue("cyt"));
					if (numberSigmas < 0) {
						throw new Exception();
					}
					log.info("Generate Cytoscape network = " + generateXGMMLFiles);
				} catch (final Exception e) {
					final String errorMessage = "Invalid cyt value '" + cmd.getOptionValue("cyt")
							+ "'. Use true or false.";
					throw new Exception(errorMessage);
				}
			} else {
				log.info("cyt parameter is not provided. Using " + generateXGMMLFiles + " by default.");
			}

			boolean ignoreNotFoundSignals = defaultIgnoreNotFoundSignals; // by
																			// default
			if (cmd.hasOption("ig")) {
				try {
					ignoreNotFoundSignals = Boolean.valueOf(cmd.getOptionValue("ig"));
					if (numberSigmas < 0) {
						throw new Exception();
					}
					log.info("Ignore not found signals parameter = " + ignoreNotFoundSignals);
				} catch (final Exception e) {
					final String errorMessage = "Invalid ig value '" + cmd.getOptionValue("ig")
							+ "'. Use true or false.";
					throw new Exception(errorMessage);
				}
			} else {
				log.info("ig parameter is not provided. Using " + ignoreNotFoundSignals + " by default.");
			}
			final List<File> tmtFiles = Files.readAllLines(Paths.get(inputFiles.toURI())).stream()
					.map(fullPath -> new File(fullPath)).collect(Collectors.toList());
			final TMTQuantSiteComparator runner = new TMTQuantSiteComparator(paramFile, tmtFiles, tmtType, rInf,
					outputFileName, pValueCorrectionType, qValueThreshold, numberSigmas, minNumberOfDiscoveries,
					generateXGMMLFiles, ignoreNotFoundSignals);
			runner.run();
			System.out.println("Program finished successfully.");
			System.exit(0);
		} catch (final Exception e) {
			errorInParameters();
			e.printStackTrace();
			System.out.println("Program finished with some error: " + e.getMessage());
			System.exit(-1);
		}
	}

	private void run() throws IOException {
		// get PCQ parameters and pairwise ratio files
		final TMTPairWisePCQInputParametersGenerator pcqInputParamtersGenerator = new TMTPairWisePCQInputParametersGenerator(
				paramFile, tmtFiles, tmtType, outputFileName, ignoreNotFoundSignals);
		final Map<String, File> pcqParameters = pcqInputParamtersGenerator.run();
		// run pcq in batch
		final PCQBatchRunner pcqBatchRunner = new PCQBatchRunner(pcqParameters, generateXGMMLFiles);
		final Map<String, File> outputFolders = pcqBatchRunner.run();
		// grab pcq outputs
		final List<File> pcqOutputQuantPerSiteFileList = getPCQOutputQuantPerSiteFile(outputFolders);
		// use them for quant site output
		final QuantSiteOutputComparator comparator = new QuantSiteOutputComparator(pcqOutputQuantPerSiteFileList, rInf,
				outputFileName, pValueCorrectionType, qValueThreshold, numberSigmas, minNumberOfDiscoveries, true);
		comparator.setOutputFolder(paramFile.getParent());
		comparator.run();

	}

	private List<File> getPCQOutputQuantPerSiteFile(Map<String, File> outputFolders) {
		final List<File> ret = new ArrayList<File>();
		for (final String expName : outputFolders.keySet()) {
			final File[] listFiles = outputFolders.get(expName)
					.listFiles(QuantSiteOutputComparator.getPeptideNodeFileFilter());
			if (listFiles.length != 1) {
				throw new IllegalArgumentException("Not able to find peptideNode output file in folder '"
						+ outputFolders.get(expName).getAbsolutePath() + "' for experiment " + expName);
			}
			ret.add(listFiles[0]);
			///////////////////////////////////////////////////////////////////////
			// not proud of this use of static: but it is ok
			// this will be used then in quantsiteoutputcomparator
			final String sampleName = expName;

			QuantSiteOutputComparator.sampleNamesByFiles.put(listFiles[0], sampleName);
			///////////////////////////////////////////////////////////////////////
		}
		return ret;
	}

	public TMTQuantSiteComparator(File paramFile, List<File> tmtFiles, String tmtType, double rInf,
			String outputFileName, PValueCorrectionType pValueCorrectionType, double qValueThreshold, int numberSigmas,
			int minNumberOfDiscoveries, boolean generateXGMMLFiles2, boolean ignoreNotFoundSignals) throws IOException {
		this.paramFile = paramFile;
		this.tmtFiles = tmtFiles;
		this.tmtType = tmtType;
		this.outputFileName = outputFileName;
		this.pValueCorrectionType = pValueCorrectionType;
		this.minNumberOfDiscoveries = minNumberOfDiscoveries;
		this.numberSigmas = numberSigmas;
		this.qValueThreshold = qValueThreshold;
		this.rInf = rInf;
		generateXGMMLFiles = generateXGMMLFiles2;
		this.ignoreNotFoundSignals = ignoreNotFoundSignals;
	}

	private static void setupCommandLineOptions() {
		// create Options object
		options = new Options();
		final Option opt1 = new Option("pf", "param_file", true, "[MANDATORY] -pf input parameter file used as base");
		opt1.setRequired(true);
		options.addOption(opt1);
		final Option opt2 = new Option("ifs", "input_files", true,
				"[MANDATORY] Full path to a file containing a list of full paths to tmt census out files to use");
		opt2.setRequired(true);
		options.addOption(opt2);
		final Option opt3 = new Option("tmt", "tmt_type", true, "[OPTIONAL] Either 10PLEX (by default) or 6PLEX");
		opt3.setRequired(false);
		options.addOption(opt3);

		// create Options object

		final Option opt4 = new Option("RInf", "replace_infinity", true,
				"[OPTIONAL] -RInf replaces +/- Infinity with a user defined (+/-) value in the output summary table file");
		opt4.setRequired(false);
		options.addOption(opt4);
		final Option opt6 = new Option("out", "output_file_name", true,
				"[MANDATORY] Output file name that will be created in the current folder");
		opt6.setRequired(true);
		options.addOption(opt6);
		final Option opt7 = new Option("pvc", "pvalue_correction", true,
				"[OPTIONAL] p-value correction method to apply. Valid values are: "
						+ PValueCorrectionType.getValuesString() + ". If not provided, the method will be "
						+ QuantSiteOutputComparator.defaultPValueCorrectionMethod + " (Reference: "
						+ QuantSiteOutputComparator.defaultPValueCorrectionMethod.getReference() + ")");
		opt7.setRequired(false);
		options.addOption(opt7);
		final Option opt8 = new Option("qvt", "qvalue_threshold", true,
				"[OPTIONAL] q-value threshold to apply to the corrected p-values. A value between 0 and 1 is permitted. If not provided, a threshold of "
						+ QuantSiteOutputComparator.defaultQValueThreshold + " will be applied.");
		opt8.setRequired(false);
		options.addOption(opt8);
		final Option opt9 = new Option("md", "minimum_discoveries", true,
				"[OPTIONAL] minimum number of discoveries (significantly different between two samples) required for a quantified site to be in the output files. If not provided, there will be no minimum number, although no quant sites without any significantly different site between 2 samples will be reported in the Excel output file.");
		opt9.setRequired(false);
		options.addOption(opt9);

		final Option opt10 = new Option("ns", "number_sigmas", true,
				"[OPTIONAL] number of sigmas that will be used to decide whether an INFINITY ratio is significantly different than a FINITE ratio.\n"
						+ "If R1=POSITIVE_INFINITY and R2 < avg_distribution + ns*sigma_distribution_of_ratios, then R2 is significantly different.\n"
						+ "If R1=NEGATIVE_INFINITY and R2 > avg_distribution + ns*sigma_distribution_of_ratios, then R2 is significantly different.");
		opt10.setRequired(false);
		options.addOption(opt10);

		final Option opt11 = new Option("cyt", "cytoscape_network", true,
				"[OPTIONAL] whether to create or not the cytoscape network for each of the pairwise comparisons performed with PCQ in order to speed-up the run.\n"
						+ "This parameter is " + defaultCytParameterValue + " by default.");
		opt11.setRequired(false);
		options.addOption(opt11);

	}

	private static void errorInParameters() {
		// automatically generate the help statement
		final HelpFormatter formatter = new HelpFormatter();

		formatter.printHelp(150, "java -jar TMTQuantSiteComparator.jar", "with the following parameters:", options,
				"\n\nContact Salvador Martinez-Bartolome at salvador@scripps.edu for more help");

	}
}
