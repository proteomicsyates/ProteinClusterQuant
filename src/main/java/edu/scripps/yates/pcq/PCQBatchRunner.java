package edu.scripps.yates.pcq;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.pcq.params.PropertiesReader;
import edu.scripps.yates.pcq.params.ProteinClusterQuantParameters;
import edu.scripps.yates.utilities.appversion.AppVersion;
import gnu.trove.map.hash.THashMap;

public class PCQBatchRunner {
	private final static Logger log = Logger.getLogger(PCQBatchRunner.class);
	private final Map<String, File> pcqParametersFilesMap = new THashMap<String, File>();
	private boolean generateXGMMLFiles = true;

	public static void main(String[] args) {
		final AppVersion version = ProteinClusterQuant.getVersion();
		System.out.println("Running ProteinClusterQuant (PCQ) version " + version.toString() + " in batch mode");
		ProteinClusterQuant.setupCommandLineOptions(true);
		final CommandLineParser parser = new BasicParser();
		try {
			final CommandLine cmd = parser.parse(ProteinClusterQuant.options, args);
			if (cmd.getOptionValue("f") == null) {
				throw new Exception("Provide input parameter file with 'f' option");
			}
			final File inputFile = new File(cmd.getOptionValue("f"));
			String propertiesFilePath = System.getProperty("user.dir") + File.separator
					+ ProteinClusterQuant.SETUP_PROPERTIES;
			if (inputFile != null) {
				propertiesFilePath = inputFile.getAbsolutePath();

			}
			if (cmd.getOptionValue("bf") == null) {
				throw new Exception("Provide batch file with 'bf' option");
			}
			final File batchFile = new File(cmd.getOptionValue("bf"));
			log.info("Using setup.properties file at: " + propertiesFilePath);
			log.info("Using batch file at: " + batchFile.getAbsolutePath());
			final File setupPropertiesFile = new File(propertiesFilePath);

			PCQBatchRunner.run(batchFile, setupPropertiesFile, true);
			System.out.println("Everything finished correctly.");
			System.exit(0);
		} catch (final Exception e) {
			e.printStackTrace();
			ProteinClusterQuant.errorInParameters(true);
		}
		System.err.println("Some error occurred.");
		System.exit(-1);
	}

	/**
	 * In this case, the batch file is a file containing lines with inputFiles
	 * 
	 * @param batchFile
	 * @param generateXGMMLFiles
	 * @throws IOException
	 */
	public static void run(File batchFile, File paramaterFile, boolean generateXGMMLFiles) throws IOException {

		final List<String> lines = Files.readAllLines(Paths.get(batchFile.toURI()));
		String inputFilesLine = null;
		String inputIDFilesLine = null;
		String suffixLine = null;
		int numLine = 0;
		for (final String line : lines) {
			numLine++;
			final List<String> sublines = new ArrayList<String>();
			if (line.contains(";")) {
				final String[] split = line.split(";");
				for (final String subline : split) {
					sublines.add(subline);
				}
			} else {
				sublines.add(line);
			}
			for (final String subline : sublines) {

				if (subline.contains("=")) {
					final String[] split = subline.split("=");
					final String property = split[0].trim();
					final String value = split[1].trim();
					if (property.equals("inputFiles")) {
						if (inputFilesLine != null) {
							runPCQ(paramaterFile, inputFilesLine, inputIDFilesLine, suffixLine, generateXGMMLFiles);
							// reset params
							inputFilesLine = value;
							inputIDFilesLine = null;
							suffixLine = null;
						} else {
							inputFilesLine = value;
						}
					} else if (property.equals("inputIDFiles")) {
						inputIDFilesLine = value;
					} else if (property.equals("outputSuffix")) {
						suffixLine = value;
					} else {
						log.info("line " + numLine + " ignored ('" + subline + "')");
					}
				} else {
					log.info("line " + numLine + " ignored ('" + subline + "')");
				}
			}
		}
		if (inputFilesLine != null) {
			runPCQ(paramaterFile, inputFilesLine, inputIDFilesLine, suffixLine, generateXGMMLFiles);

		}
	}

	private static void runPCQ(File parameterFile, String inputFilesLine, String inputIDFilesLine, String suffixLine,
			boolean generateXGMMLFiles2) throws IOException {
		log.info("Running PCQ with:\ninputFiles=" + inputFilesLine + "\ninputIDFiles=" + inputIDFilesLine
				+ "\noutputSuffix=" + suffixLine + "\ngenerateXGMML=" + generateXGMMLFiles2);
		// this sets the pcq parameters
		PropertiesReader.readProperties(parameterFile, true);
		final ProteinClusterQuantParameters params = ProteinClusterQuantParameters.getInstance();
		params.clearInputIdentificationFiles();
		params.clearInputQuantificationFiles();
		if (inputFilesLine != null) {
			PropertiesReader.parseInputFiles(inputFilesLine, params);
		}
		if (inputIDFilesLine != null) {
			PropertiesReader.parseInputIDFiles(inputIDFilesLine, params);
		}
		if (suffixLine != null) {
			params.setOutputSuffix(suffixLine.trim());
		}
		final ProteinClusterQuant pcq = new ProteinClusterQuant(params, parameterFile);
		params.setAnalysisRun(true);
		pcq.setCreateXGMMLFile(generateXGMMLFiles2);
		pcq.run();
		final File peptideNodeTableFile = pcq.getFinalPeptideNodeTableFile();
		// append this full path to a file that can be used by PCQ comparator
		final File output = new File(params.getOutputFileFolder().getParentFile().getParentFile().getAbsolutePath()
				+ File.separator + "QuantSiteOutputComparator_input_file.txt");
		final FileWriter fw = new FileWriter(output, true);
		final String peptideNodeTableFileLine = suffixLine + "\t" + peptideNodeTableFile.getAbsolutePath();
		if (output.length() > 0l) {
			fw.write("\n");
		}
		fw.write(peptideNodeTableFileLine);
		fw.close();
		log.info("Added to " + output.getAbsolutePath() + ": '" + peptideNodeTableFileLine + "'");
	}

	/**
	 * Constructor using generateXGMMLFiles as TRUE by default
	 * 
	 * @param pcqParameters
	 */
	public PCQBatchRunner(List<File> pcqParameters) {
		this(pcqParameters, true);
	}

	public PCQBatchRunner(List<File> pcqParameters, boolean generateXGMMLFiles) {
		this.generateXGMMLFiles = generateXGMMLFiles;
		final String expName = "exp";
		int num = 1;
		for (final File file : pcqParameters) {
			try {
				// to validate file
				PropertiesReader.readProperties(file, false);
			} catch (final IOException e) {
				e.printStackTrace();
				log.error("Error reading pcq parameters file '" + file.getAbsolutePath() + "': " + e.getMessage());
			}
			pcqParametersFilesMap.put(expName + num, file);
			num++;
		}
	}

	public PCQBatchRunner(Map<String, File> pcqParameters, boolean generateXGMMLFiles) {
		this.generateXGMMLFiles = generateXGMMLFiles;
		for (final String expName : pcqParameters.keySet()) {
			final File pcqParameterFile = pcqParameters.get(expName);
			// to validate file
			try {
				PropertiesReader.readProperties(pcqParameterFile, false);
			} catch (final IOException e) {
				e.printStackTrace();
				log.error("Error reading pcq parameters file '" + pcqParameterFile.getAbsolutePath() + "': "
						+ e.getMessage());
			}

			pcqParametersFilesMap.put(expName, pcqParameterFile);
		}
		log.info(pcqParametersFilesMap.size() + " pcq parameters files received.");

	}

	public Map<String, File> runFromBatchFile() throws IOException {
		final Map<String, File> outputFolders = new THashMap<String, File>();
		int count = 1;
		for (final String expName : pcqParametersFilesMap.keySet()) {
			final File setupPropertiesFile = pcqParametersFilesMap.get(expName);
			log.info("Running PCQ with file " + FilenameUtils.getName(setupPropertiesFile.getAbsolutePath()) + " ("
					+ expName + ") " + count++ + "/" + pcqParametersFilesMap.size());
			final ProteinClusterQuant pcq = new ProteinClusterQuant(setupPropertiesFile, true);
			pcq.setCreateXGMMLFile(generateXGMMLFiles);
			pcq.run();
			final File outputFolder = pcq.getParams().getOutputFileFolder();
			outputFolders.put(expName, outputFolder);
		}
		log.info("All PCQ runs are done with no problems");
		return outputFolders;
	}

	public Map<String, File> run() throws IOException {
		final Map<String, File> outputFolders = new THashMap<String, File>();
		int count = 1;
		for (final String expName : pcqParametersFilesMap.keySet()) {
			final File setupPropertiesFile = pcqParametersFilesMap.get(expName);
			log.info("Running PCQ with file " + FilenameUtils.getName(setupPropertiesFile.getAbsolutePath()) + " ("
					+ expName + ") " + count++ + "/" + pcqParametersFilesMap.size());
			final ProteinClusterQuant pcq = new ProteinClusterQuant(setupPropertiesFile, true);
			pcq.setCreateXGMMLFile(generateXGMMLFiles);
			pcq.run();
			final File outputFolder = pcq.getParams().getOutputFileFolder();
			outputFolders.put(expName, outputFolder);
		}
		log.info("All PCQ runs are done with no problems");
		return outputFolders;
	}

	public boolean isGenerateXGMMLFiles() {
		return generateXGMMLFiles;
	}

	public void setGenerateXGMMLFiles(boolean generateXGMMLFiles) {
		this.generateXGMMLFiles = generateXGMMLFiles;
	}
}
