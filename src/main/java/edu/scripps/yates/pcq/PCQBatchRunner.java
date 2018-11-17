package edu.scripps.yates.pcq;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.pcq.params.PropertiesReader;
import gnu.trove.map.hash.THashMap;

public class PCQBatchRunner {
	private final static Logger log = Logger.getLogger(PCQBatchRunner.class);
	private final Map<String, File> pcqParametersMap = new THashMap<String, File>();

	public PCQBatchRunner(List<File> pcqParameters) {
		final String expName = "exp";
		int num = 1;
		for (final File file : pcqParameters) {
			try {
				// to validate file
				PropertiesReader.readProperties(file);
			} catch (final IOException e) {
				e.printStackTrace();
				log.error("Error reading pcq parameters file '" + file.getAbsolutePath() + "': " + e.getMessage());
			}
			pcqParametersMap.put(expName + num, file);
			num++;
		}
	}

	public PCQBatchRunner(Map<String, File> pcqParameters) {
		for (final String expName : pcqParameters.keySet()) {
			final File pcqParameterFile = pcqParameters.get(expName);
			// to validate file
			try {
				PropertiesReader.readProperties(pcqParameterFile);
			} catch (final IOException e) {
				e.printStackTrace();
				log.error("Error reading pcq parameters file '" + pcqParameterFile.getAbsolutePath() + "': "
						+ e.getMessage());
			}

			pcqParametersMap.put(expName, pcqParameterFile);
		}
		log.info(pcqParametersMap.size() + " pcq parameters files received.");

	}

	public Map<String, File> run() throws IOException {
		final Map<String, File> outputFolders = new THashMap<String, File>();
		int count = 1;
		for (final String expName : pcqParametersMap.keySet()) {
			final File setupPropertiesFile = pcqParametersMap.get(expName);
			log.info("Running PCQ with file " + FilenameUtils.getName(setupPropertiesFile.getAbsolutePath()) + " ("
					+ expName + ") " + count++ + "/" + pcqParametersMap.size());
			final ProteinClusterQuant pcq = new ProteinClusterQuant(setupPropertiesFile, true);
			pcq.run();
			final File outputFolder = pcq.getParams().getOutputFileFolder();
			outputFolders.put(expName, outputFolder);
		}
		log.info("All PCQ runs are done with no problems");
		return outputFolders;
	}
}
