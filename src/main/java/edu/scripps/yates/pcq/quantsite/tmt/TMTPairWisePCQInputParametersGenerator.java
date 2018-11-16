package edu.scripps.yates.pcq.quantsite.tmt;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.CensusOutParser;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.pcq.ProteinClusterQuant;
import edu.scripps.yates.pcq.params.PropertiesReader;
import edu.scripps.yates.pcq.util.AnalysisInputType;
import edu.scripps.yates.pcq.util.ExperimentFiles;
import edu.scripps.yates.utilities.appversion.AppVersion;
import edu.scripps.yates.utilities.model.enums.AmountType;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import gnu.trove.map.hash.THashMap;

/**
 * It gets a PCQ input parameters file and a set of TMT files and generates a
 * new input parameter file per pairwise comparison between the channels,
 * pointing to new input files, tab separated values type, with the ratios
 * calculated from the intensities parsed.
 * 
 * @author salvador
 *
 */
public class TMTPairWisePCQInputParametersGenerator {
	private final static Logger log = Logger.getLogger(TMTPairWisePCQInputParametersGenerator.class);
	private static AppVersion version;
	private static Options options;
	private final File baseParamFile;
	private final List<File> tmtFiles;
	private final String tmtType;
	private Map<QuantCondition, QuantificationLabel> labelsByConditions;
	private final String outputFileName;
	public final static String TMT10PLEX = "10PLEX";
	public final static String TMT6PLEX = "6PLEX";
	private static final String TMT_DATA_FILES = "tmt_pairwise_data_files";
	private static final String PCQ_PARAMETERS = "pcq_parameters";

	public TMTPairWisePCQInputParametersGenerator(File paramFile, List<File> tmtFiles, String tmtType,
			String outputFileName) {
		baseParamFile = paramFile;
		this.tmtFiles = tmtFiles;
		this.tmtType = tmtType;
		this.outputFileName = outputFileName;
	}

	public static void main(String[] args) {
		version = ProteinClusterQuant.getVersion();
		System.out.println("Running TMTPairWisePCQInputParametersGenerator version " + version.toString());
		setupCommandLineOptions();
		final CommandLineParser parser = new BasicParser();
		try {
			final CommandLine cmd = parser.parse(options, args);

			final File paramFile = new File(cmd.getOptionValue("pf"));
			final File inputFiles = new File(cmd.getOptionValue("ifs"));
			String tmtType = TMT10PLEX; // by default
			if (cmd.hasOption("tmt")) {
				tmtType = cmd.getOptionValue("tmt");
			}
			if (!TMT10PLEX.equals(tmtType) && !TMT6PLEX.equals(tmtType)) {
				throw new IllegalArgumentException("Invalid value for tmt parameter: '" + tmtType
						+ "'. Valid values are " + TMT10PLEX + " or " + TMT6PLEX + " in tmt parameter");
			}
			final String outputFileName = cmd.getOptionValue("out");
			final List<File> tmtFiles = Files.readAllLines(Paths.get(inputFiles.toURI())).stream()
					.map(fullPath -> new File(fullPath)).collect(Collectors.toList());
			final TMTPairWisePCQInputParametersGenerator script = new TMTPairWisePCQInputParametersGenerator(paramFile,
					tmtFiles, tmtType, outputFileName);
			script.run();
			System.exit(0);
		} catch (final Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
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

	}

	/**
	 * 
	 * @return
	 * @throws IOException
	 */
	public List<File> run() throws IOException {
		log.info("Running " + getClass().getCanonicalName());
		final List<File> pcqParamtersFiles = new ArrayList<File>();
		labelsByConditions = generateLabelsByConditions();

		// first create a file per TMT
		for (final File tmtFile : tmtFiles) {
			createPairWiseTMTTSVFiles(tmtFile);
		}
		// now create a parameter file per pairwise comparison of labels
		log.info(
				"Creating a PCQ parameter file for all pairwise combinations of TMT using as a base the parameter file "
						+ baseParamFile.getAbsolutePath());

		final List<QuantificationLabel> labels = getLabelList();
		for (int i = 0; i < labels.size(); i++) {
			final QuantificationLabel labelNumerator = labels.get(i);
			for (int j = i + 1; j < labels.size(); j++) {
				final QuantificationLabel labelDenominator = labels.get(j);
				final List<File> pairWiseTSVPCQInputFiles = new ArrayList<File>();
				for (final File tmtFile : tmtFiles) {
					final File pairWiseTSVPCQInputFile = getOutputFileForTSVData(tmtFile, labelNumerator,
							labelDenominator);
					log.info("Using tsv file " + pairWiseTSVPCQInputFile.getAbsolutePath());
					if (!pairWiseTSVPCQInputFile.exists()) {
						throw new IllegalArgumentException(
								"File '" + pairWiseTSVPCQInputFile.getAbsolutePath() + "' has not been created");
					}
					pairWiseTSVPCQInputFiles.add(pairWiseTSVPCQInputFile);
				}
				final File parameterFile = createPCQParameterFile(pairWiseTSVPCQInputFiles, labelNumerator,
						labelDenominator);
				pcqParamtersFiles.add(parameterFile);
				System.out.println(
						"PCQ parameter file created: " + FilenameUtils.getName(parameterFile.getAbsolutePath()));

			}
		}

		return pcqParamtersFiles;
	}

	private File createPCQParameterFile(List<File> pairWiseTSVPCQInputFiles, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator) throws IOException {
		log.info("Creating PCQ Parameter file for comparison: " + labelNumerator + " vs " + labelDenominator);
		log.info("Reading base parameter file: " + baseParamFile.getAbsolutePath());
		final List<String> paramLines = Files.readAllLines(Paths.get(baseParamFile.getAbsolutePath()));
		log.info(paramLines.size() + " lines read from parameter file");

		final File outputParamFile = getOutputFileForPCQParameters(baseParamFile, labelNumerator, labelDenominator);
		log.info("File to be created at: " + outputParamFile.getAbsolutePath());

		FileWriter fw = null;
		boolean containsInputType = false;
		boolean containsQuantChannels = false;
		boolean containsInputFilePath = false;
		boolean containsInputFiles = false;
		boolean containsOutputSuffix = false;
		final String parentPath = Paths.get(pairWiseTSVPCQInputFiles.get(0).getAbsolutePath()).getParent()
				.toAbsolutePath().toString().replace("\\", "\\\\");
		String originalInputParentPath = null;
		try {
			log.info("Creating PCQ parameters file for TMT comparison " + labelNumerator + " vs " + labelDenominator);

			fw = new FileWriter(outputParamFile);
			for (final String paramLine : paramLines) {
				try {
					if (isInputType(paramLine)) {
						fw.write("inputType = " + AnalysisInputType.SEPARATED_VALUES.name() + "\n");
						containsInputType = true;
					} else if (isQuantChannels(paramLine)) {
						fw.write("quantChannels = " + labelNumerator.name() + "/" + labelDenominator + "\n");
						containsQuantChannels = true;
					} else if (isInputFilePath(paramLine)) {
						originalInputParentPath = paramLine.split("=")[1].trim();
						fw.write("inputFilePath = " + parentPath);
						containsInputFilePath = true;
					} else if (isInputFiles(paramLine)) {
						fw.write("inputFiles = exp[");
						for (final File file : pairWiseTSVPCQInputFiles) {
							fw.write(FilenameUtils.getName(file.getAbsolutePath()) + ",");
						}
						fw.write("]");
						containsInputFiles = true;
					} else if (isOutputSuffix(paramLine)) {
						String suffix = "";
						if (paramLine.trim().contains("=")) {
							suffix = paramLine.trim().split("=")[1].trim();
						}
						suffix += labelNumerator.name() + "_vs_" + labelDenominator;
						fw.write("outputSuffix = " + suffix);
						containsOutputSuffix = true;
					} else if (isInputIDFiles(paramLine)) {
						// copy input ID files
						copyInputIDFiles(paramLine, new File(originalInputParentPath), labelNumerator,
								labelDenominator);
					} else {
						fw.write(paramLine);
					}
				} finally {
					fw.write("\n");
				}
			}
			log.info("PCQ parameters file created for TMT comparison " + labelNumerator + " vs " + labelDenominator);

			return outputParamFile;
		} finally {
			if (!containsInputType) {
				fw.write("inputType = " + AnalysisInputType.SEPARATED_VALUES.name() + "\n");
			}
			if (!containsQuantChannels) {
				fw.write("quantChannels = " + labelNumerator.name() + "/" + labelDenominator + "\n");
			}
			if (!containsInputFilePath) {
				fw.write("inputFilePath = " + parentPath + "\n");
			}
			if (!containsInputFiles) {
				fw.write("inputFiles = exp[");
				for (final File file : pairWiseTSVPCQInputFiles) {
					fw.write(FilenameUtils.getName(file.getAbsolutePath()) + ",");
				}
				fw.write("]");
			}
			if (!containsOutputSuffix) {
				final String suffix = labelNumerator.name() + "_vs_" + labelDenominator;
				fw.write("outputSuffix = " + suffix + "\n");
			}
			if (fw != null) {
				fw.close();
			}
		}

	}

	private void copyInputIDFiles(String paramLine, File originalInputFileFolder, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator) throws IOException {
		final String fileNamesString = paramLine.split("=")[1].trim();
		if (fileNamesString != null) {
			final List<ExperimentFiles> expFiles = new ArrayList<ExperimentFiles>();
			if (fileNamesString.contains("|")) {
				final String[] tmp = fileNamesString.split("\\|");
				for (int i = 0; i < tmp.length; i++) {
					final ExperimentFiles experimentFiles = PropertiesReader.parseExperimentFileNames(tmp[i].trim());
					expFiles.add(experimentFiles);
				}
			} else {
				if (!"".equals(fileNamesString)) {
					final ExperimentFiles experimentFiles = PropertiesReader.parseExperimentFileNames(fileNamesString);
					expFiles.add(experimentFiles);
				}
			}
			for (final ExperimentFiles experimentFiles : expFiles) {
				final List<String> files = experimentFiles.getRelicateFileNames();
				for (final String file : files) {
					final File tsvDataFile = getOutputFileForTSVData(originalInputFileFolder, labelNumerator,
							labelDenominator);
					final String path = Paths.get(tsvDataFile.getAbsolutePath()).getParent().toString();
					final File destFile = new File(path + File.separator + FilenameUtils.getName(file));

					final File originFile = new File(
							originalInputFileFolder.getAbsolutePath() + File.separator + FilenameUtils.getName(file));

					if (!destFile.exists() || destFile.length() != originFile.length()) {
						org.apache.commons.io.FileUtils.copyFile(originFile, destFile);
					}
				}
			}
		}
	}

	private boolean isInputIDFiles(String paramLine) {
		return !isCommentLine(paramLine) && paramLine.trim().startsWith("inputIDFiles");
	}

	private boolean isOutputSuffix(String paramLine) {
		return !isCommentLine(paramLine) && paramLine.trim().startsWith("outputSuffix");
	}

	private boolean isInputFiles(String paramLine) {
		return !isCommentLine(paramLine) && paramLine.trim().startsWith("inputFiles");
	}

	private boolean isInputFilePath(String paramLine) {
		return !isCommentLine(paramLine) && paramLine.trim().startsWith("inputFilePath");
	}

	private boolean isQuantChannels(String paramLine) {
		return !isCommentLine(paramLine) && paramLine.trim().startsWith("quantChannels");
	}

	private boolean isInputType(String paramLine) {
		return !isCommentLine(paramLine) && paramLine.trim().startsWith("inputType");
	}

	private boolean isCommentLine(String paramLine) {
		return paramLine.trim().startsWith("#");
	}

	private void createPairWiseTMTTSVFiles(File tmtFile) throws IOException {
		log.info("Creating data file (TSV format) for all pairwise combinations of TMT in file "
				+ tmtFile.getAbsolutePath());
		final CensusOutParser parser = new CensusOutParser(tmtFile, labelsByConditions, null, null);
		final List<String> psmKeys = new ArrayList<String>();
		psmKeys.addAll(parser.getPSMMap().keySet());
		Collections.sort(psmKeys);
		final List<QuantificationLabel> labels = getLabelList();
		// create a file for each pairwise comparison between channels
		for (int i = 0; i < labels.size(); i++) {
			final QuantificationLabel labelNumerator = labels.get(i);
			for (int j = i + 1; j < labels.size(); j++) {
				final QuantificationLabel labelDenominator = labels.get(j);
				log.info("Creating data file (TSV format) for " + labelNumerator + " vs " + labelDenominator);

				final File outputTSVFile = getOutputFileForTSVData(tmtFile, labelNumerator, labelDenominator);
				FileWriter fw = null;
				try {
					fw = new FileWriter(outputTSVFile);
					for (final String psmKey : psmKeys) {
						final QuantifiedPSMInterface psm = parser.getPSMMap().get(psmKey);
						final Set<QuantifiedProteinInterface> quantifiedProteins = psm.getQuantifiedProteins();
						// one line, the same, per protein
						for (final QuantifiedProteinInterface quantifiedProteinInterface : quantifiedProteins) {

							fw.write(psm.getKey() + "\t");
							fw.write(psm.getFullSequence() + "\t");
							double intensityNumerator = Double.NaN;
							double intensityDenominator = Double.NaN;
							for (final Amount amount : psm.getAmounts()) {
								if (amount.getAmountType() == AmountType.NORMALIZED_INTENSITY) {
									final QuantificationLabel quantificationLabel = labelsByConditions
											.get(amount.getCondition());
									if (quantificationLabel == labelNumerator) {
										intensityNumerator = amount.getValue();
									} else if (quantificationLabel == labelDenominator) {
										intensityDenominator = amount.getValue();
									}
								}
							}
							// ratio
							if (Double.isNaN(intensityNumerator) || Double.isNaN(intensityDenominator)) {
								fw.write(Double.NaN + "\t");
							} else {
								final double ratio = intensityNumerator / intensityDenominator;
								fw.write(String.valueOf(ratio) + "\t");
							}
							// weight
							double max = 0.0;
							for (final Amount amount : psm.getAmounts()) {
								if (amount.getAmountType() == AmountType.NORMALIZED_INTENSITY) {
									max = Math.max(max, amount.getValue());
								}
							}
							fw.write(max + "\t");
							// protein
							fw.write(quantifiedProteinInterface.getAccession());
							fw.write("\n");
						}
					}
					log.info("Data file (TSV format) for " + labelNumerator + " vs " + labelDenominator + " created");
				} finally {
					if (fw != null) {
						fw.close();
					}
				}
			}
		}
	}

	private List<QuantificationLabel> getLabelList() {
		final List<QuantificationLabel> labels = new ArrayList<QuantificationLabel>();
		if (isTMT6Plex()) {
			labels.addAll(QuantificationLabel.getTMT6PlexLabels());
		} else if (isTMT10Plex()) {
			labels.addAll(QuantificationLabel.getTMT10PlexLabels());
		} else {
			throw new IllegalArgumentException("Non valid value for tmt parameter: '" + tmtType + "'. Valid values are "
					+ TMT10PLEX + " and " + TMT6PLEX);
		}
		labels.sort(new Comparator<QuantificationLabel>() {

			@Override
			public int compare(QuantificationLabel o1, QuantificationLabel o2) {
				return Integer.compare(o1.ordinal(), o2.ordinal());
			}
		});
		return labels;
	}

	private File getOutputFileForPCQParameters(File file, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator) {
		return getOutputFile(file, PCQ_PARAMETERS, labelNumerator, labelDenominator, "properties");
	}

	private File getOutputFileForTSVData(File file, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator) {
		final File outputFile = getOutputFile(file, TMT_DATA_FILES, labelNumerator, labelDenominator, "tsv");
		return outputFile;
	}

	private File getOutputFile(File file, String folder, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator, String extension) {
		log.info("Returning file for creation, using: " + file.getAbsolutePath() + ", folder: " + folder + ", "
				+ labelNumerator + " vs " + labelDenominator + ", extension: " + extension);
		String basefolder = file.getAbsolutePath();
		if (file.isFile()) {
			basefolder = Paths.get(file.getAbsolutePath()).getParent().toAbsolutePath().toString();
		}
		log.info("base folder: " + basefolder);
		String fileFullPath = basefolder + File.separator;
		if (folder != null && !"".equals(folder)) {
			fileFullPath += folder + File.separator;
		}
		log.info("fileFullPath: " + fileFullPath);
		fileFullPath += FilenameUtils.getBaseName(file.getAbsolutePath()) + "_" + labelNumerator.name() + "_vs_"
				+ labelDenominator.name() + "." + extension;
		log.info("fileFullPath: " + fileFullPath);
		final File ret = new File(fileFullPath);
		log.info("file " + ret.getAbsolutePath());
		if (!ret.getParentFile().exists()) {
			ret.getParentFile().mkdirs();
		}
		return ret;
	}

	private boolean isTMT10Plex() {
		return TMT10PLEX.equals(tmtType);
	}

	private boolean isTMT6Plex() {
		return TMT6PLEX.equals(tmtType);
	}

	private Map<QuantCondition, QuantificationLabel> generateLabelsByConditions() {
		final Map<QuantCondition, QuantificationLabel> labelsByConditions = new THashMap<QuantCondition, QuantificationLabel>();
		Set<QuantificationLabel> labels = null;
		if (TMT10PLEX.equals(tmtType)) {
			labels = QuantificationLabel.getTMT10PlexLabels();
		} else if (TMT6PLEX.equals(tmtType)) {
			labels = QuantificationLabel.getTMT6PlexLabels();
		}
		for (final QuantificationLabel label : labels) {
			switch (label) {
			case TMT_10PLEX_126_127726:
				labelsByConditions.put(getCondition("AD1"), label);
				break;
			case TMT_10PLEX_127_124761:
				labelsByConditions.put(getCondition("AD2"), label);
				break;
			case TMT_10PLEX_127_131081:
				labelsByConditions.put(getCondition("AD3"), label);
				break;
			case TMT_10PLEX_128_128116:
				labelsByConditions.put(getCondition("AD4"), label);
				break;
			case TMT_10PLEX_128_134436:
				labelsByConditions.put(getCondition("AD5"), label);
				break;
			case TMT_10PLEX_129_131471:
				labelsByConditions.put(getCondition("C1"), label);
				break;
			case TMT_10PLEX_129_13779:
				labelsByConditions.put(getCondition("C2"), label);
				break;
			case TMT_10PLEX_130_134825:
				labelsByConditions.put(getCondition("C3"), label);
				break;
			case TMT_10PLEX_130_141145:
				labelsByConditions.put(getCondition("C4"), label);
				break;
			case TMT_10PLEX_131_13818:
				labelsByConditions.put(getCondition("C5"), label);
				break;
			case TMT_6PLEX_126:
				labelsByConditions.put(getCondition("AD1"), label);
				break;
			case TMT_6PLEX_127:
				labelsByConditions.put(getCondition("AD2"), label);
				break;
			case TMT_6PLEX_128:
				labelsByConditions.put(getCondition("AD3"), label);
				break;
			case TMT_6PLEX_129:
				labelsByConditions.put(getCondition("C1"), label);
				break;
			case TMT_6PLEX_130:
				labelsByConditions.put(getCondition("C2"), label);
				break;
			case TMT_6PLEX_131:
				labelsByConditions.put(getCondition("C3"), label);
				break;
			default:
				break;
			}
		}
		return labelsByConditions;
	}

	private QuantCondition getCondition(String condName) {
		return new QuantCondition(condName);
	}
}
