package edu.scripps.yates.pcq.sanxot;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import org.apache.commons.exec.CommandLine;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.census.analysis.QuantAnalysis;
import edu.scripps.yates.census.analysis.QuantAnalysis.ANALYSIS_LEVEL_OUTCOME;
import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.analysis.QuantExperiment;
import edu.scripps.yates.census.analysis.QuantParameters;
import edu.scripps.yates.census.analysis.QuantReplicate;
import edu.scripps.yates.census.analysis.QuantificationType;
import edu.scripps.yates.census.analysis.SanXotInterfaze;
import edu.scripps.yates.census.analysis.wrappers.IntegrationResultWrapper;
import edu.scripps.yates.census.analysis.wrappers.OutlierRemovalResultWrapper;
import edu.scripps.yates.census.analysis.wrappers.SanXotAnalysisResult;
import edu.scripps.yates.census.read.model.interfaces.QuantParser;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.pcq.ProteinClusterQuant;
import edu.scripps.yates.pcq.params.ProteinClusterQuantParameters;
import edu.scripps.yates.pcq.util.ExperimentFiles;
import edu.scripps.yates.pcq.util.PCQUtils;
import edu.scripps.yates.utilities.exec.ProcessExecutor;
import edu.scripps.yates.utilities.exec.ProcessExecutorHandler;

public class SanxotRunner {
	private final QuantAnalysis quantAnalysis;
	private final static Logger log = Logger.getLogger(SanxotRunner.class);
	public static final long DEFAULT_TIMEOUT = 1000 * 60 * 20;// 20 min

	public SanxotRunner(ProteinClusterQuant pcq, QuantificationType quantType, File workingFolder,
			QuantCondition condition1, QuantCondition condition2, File fastaFile, QuantParameters quantParameters,
			Set<String> peptideInclusionList) throws FileNotFoundException {
		final ProteinClusterQuantParameters params = ProteinClusterQuantParameters.getInstance();

		quantAnalysis = new QuantAnalysis(quantType, workingFolder, condition1, condition2,
				ANALYSIS_LEVEL_OUTCOME.PEPTIDE);
		// set the ratio name
		quantParameters.setRatioName(PCQUtils.getRatioNameByAnalysisType());
		quantAnalysis.setQuantParameters(quantParameters);
		for (final ExperimentFiles experimentFiles : params.getInputQuantificationFileNames()) {
			final QuantExperiment quantExperiment = new QuantExperiment(experimentFiles.getExperimentName());
			for (final String replicateFileName : experimentFiles.getRelicateFileNames()) {
				final Map<QuantCondition, QuantificationLabel> labelsByConditions = pcq
						.getLabelsByConditions(replicateFileName);
				final QuantParser quantParser = getQuantParser(replicateFileName, labelsByConditions,
						peptideInclusionList);
				final QuantReplicate replicate = new QuantReplicate(replicateFileName, quantParser, labelsByConditions);
				quantExperiment.addReplicate(replicate);
			}
			quantAnalysis.addQuantExperiment(quantExperiment);
		}
		quantAnalysis.setKeepExperimentsSeparated(true);

		// QuantAnalysis analysis = new QuantAnalysis(quantType, workingFolder,
		// condition1, condition2,
		// ANALYSIS_LEVEL_OUTCOME.PEPTIDE);
		// analysis.runSanXot();
		//
		// final SanXotAnalysisResult result = analysis.getResult();
		// final Map<String, SanxotQuantResult> highLevelRatios =
		// result.getLastIntegrationResults()
		// .getLowerLevelRatiosAndNormalizedWeigths();
		// for (String peptide : highLevelRatios.keySet()) {
		// System.out.println(peptide + "\t" +
		// highLevelRatios.get(peptide).getLog2ratio() + "\t"
		// + highLevelRatios.get(peptide).getWeight());
		// }
	}

	private QuantParser getQuantParser(String replicateFileName,
			Map<QuantCondition, QuantificationLabel> labelsByConditions, Set<String> peptideInclusionList)
			throws FileNotFoundException {
		final QuantParser parser = PCQUtils.getQuantParser(ProteinClusterQuantParameters.getInstance(),
				labelsByConditions, replicateFileName, true, peptideInclusionList, false);
		return parser;
	}

	public SanXotAnalysisResult run() throws IOException {
		// run sanxot until peptide level and asign a
		quantAnalysis.runSanXot();
		final SanXotAnalysisResult result = quantAnalysis.getResult();

		return result;
	}

	/**
	 * @return the quantAnalysis
	 */
	public QuantAnalysis getQuantAnalysis() {
		return quantAnalysis;
	}

	private static Long runCommand(CommandLine commandLine, long timeout)
			throws IOException, InterruptedException, ExecutionException {
		final String commandString = commandLine.toString();
		log.info("Running: " + commandString);

		final ProcessExecutorHandler handler = new ProcessExecutorHandler() {

			@Override
			public void onStandardOutput(String msg) {
				log.debug("OUTPUT:" + msg);

			}

			@Override
			public void onStandardError(String msg) {
				log.error("ERROR:" + msg);

			}
		};
		final Future<Long> runProcess = ProcessExecutor.runProcess(commandLine, handler, timeout);
		while (!runProcess.isDone() && !runProcess.isCancelled()) {
			Thread.sleep(1000);
		}
		final Long processExitCode = runProcess.get();
		log.info("Process exitValue: " + processExitCode);
		return processExitCode;
	}

	public static IntegrationResultWrapper integrate(File relatFile, File dataFile, File infoFile, String prefix,
			Double forzedVariance, File workingFolder, boolean checkRelationshipValidity,
			QuantParameters quantParameters) throws IOException, InterruptedException, ExecutionException {

		if (checkRelationshipValidity && !SanXotInterfaze.checkDataValidity(relatFile, dataFile)) {
			throw new IllegalArgumentException("Combination of data file and relat file is not valid: "
					+ FilenameUtils.getName(relatFile.getAbsolutePath()) + " and "
					+ FilenameUtils.getName(dataFile.getAbsolutePath()));
		}

		final String msg = "Integrating data ...";
		log.info(msg);
		if (relatFile == null) {
			log.info("Not using relationship file. Using -C option to correct protein loading error");
		} else {
			log.info("Using relationship file " + FilenameUtils.getName(relatFile.getAbsolutePath()));
		}
		log.info("Using data file  " + FilenameUtils.getName(dataFile.getAbsolutePath()));
		CommandLine integratingCommandLine = SanXotInterfaze.getIntegrationCommandLine(relatFile, dataFile, infoFile,
				prefix, forzedVariance, workingFolder, quantParameters);

		final Long exitCode = runCommand(integratingCommandLine, DEFAULT_TIMEOUT);
		if (exitCode.longValue() != 0) {
			if (exitCode.longValue() == ProcessExecutor.TIMEOUT_ERROR_CODE) {
				final String message = "The process cound't finish before the timeout of " + DEFAULT_TIMEOUT + " ms";
				log.warn(message);
				log.info("Trying to fix the problem by forzing variance to 0 (Using -f v0)");
				integratingCommandLine = SanXotInterfaze.getIntegrationCommandLine(relatFile, dataFile, infoFile,
						prefix, 0.0, workingFolder, quantParameters);
				final Long newExitCode = runCommand(integratingCommandLine, DEFAULT_TIMEOUT);
				if (newExitCode.longValue() != 0) {
					if (newExitCode.longValue() == ProcessExecutor.TIMEOUT_ERROR_CODE) {
						log.warn(message);
						throw new IllegalArgumentException(message);
					}
					throw new IllegalArgumentException("Some error happen while integration process");
				}
			} else {
				throw new IllegalArgumentException("Some error happen while integration process");
			}
		}
		final IntegrationResultWrapper integrationResults = new IntegrationResultWrapper(workingFolder, prefix, -1, -1);

		log.info("Integration performed. Integration file at:"
				+ integrationResults.getHigherLevelDataFile().getAbsolutePath());
		log.info("Integration Variance=" + integrationResults.getIntegrationVariance());
		return integrationResults;
	}

	public static OutlierRemovalResultWrapper removeOutliers(File relatFile, File dataFile, File infoFile,
			String prefix, File workingFolder, QuantParameters quantParameters)
			throws IOException, InterruptedException, ExecutionException {
		final String msg = "Removing outliers...";
		log.info(msg);
		final CommandLine removeOutlierCommandLine = SanXotInterfaze.getRemoveOutliersCommandLine(relatFile, prefix,
				dataFile, infoFile, workingFolder, quantParameters);

		final Long exitCode = runCommand(removeOutlierCommandLine, DEFAULT_TIMEOUT);
		if (exitCode.longValue() != 0)
			throw new IllegalArgumentException("Some error happen while outlier removal process");
		final OutlierRemovalResultWrapper outliersRemovalResults = new OutlierRemovalResultWrapper(workingFolder,
				prefix);

		log.info("Outlier removal performed. New relation file at:"
				+ outliersRemovalResults.getRelatFile().getAbsolutePath());

		return outliersRemovalResults;

	}
}
