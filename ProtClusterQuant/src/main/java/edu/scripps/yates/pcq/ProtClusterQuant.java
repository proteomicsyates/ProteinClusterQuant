package edu.scripps.yates.pcq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;

import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.annotations.uniprot.xml.Entry;
import edu.scripps.yates.census.analysis.QuantAnalysis;
import edu.scripps.yates.census.analysis.QuantAnalysis.ANALYSIS_LEVEL_OUTCOME;
import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.analysis.QuantificationType;
import edu.scripps.yates.census.analysis.wrappers.IntegrationResultWrapper;
import edu.scripps.yates.census.analysis.wrappers.OutlierRemovalResultWrapper;
import edu.scripps.yates.census.analysis.wrappers.SanXotAnalysisResult;
import edu.scripps.yates.census.analysis.wrappers.SanxotQuantResult;
import edu.scripps.yates.census.read.AbstractQuantParser;
import edu.scripps.yates.census.read.model.CensusRatio;
import edu.scripps.yates.census.read.model.RatioScore;
import edu.scripps.yates.census.read.model.interfaces.QuantParser;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.pcq.filter.PCQFilter;
import edu.scripps.yates.pcq.model.PCQPeptideNode;
import edu.scripps.yates.pcq.model.PCQProteinNode;
import edu.scripps.yates.pcq.model.ProteinCluster;
import edu.scripps.yates.pcq.model.ProteinPair;
import edu.scripps.yates.pcq.sanxot.SanxotRunner;
import edu.scripps.yates.pcq.util.ExperimentFiles;
import edu.scripps.yates.pcq.util.InputType;
import edu.scripps.yates.pcq.util.PropertiesReader;
import edu.scripps.yates.pcq.util.ProteinClusterQuantParameters;
import edu.scripps.yates.pcq.util.Utils;
import edu.scripps.yates.pcq.xgmml.XgmmlExporter;
import edu.scripps.yates.utilities.alignment.nwalign.NWResult;
import edu.scripps.yates.utilities.dates.DatesUtil;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.model.enums.AggregationLevel;
import edu.scripps.yates.utilities.model.enums.CombinationType;
import edu.scripps.yates.utilities.proteomicsmodel.Score;

public class ProtClusterQuant {
	private final static Logger log = Logger.getLogger(ProtClusterQuant.class);
	private static final String SETUP_PROPERTIES = "setup.properties";
	private static final String sep = "\t";
	private final QuantCondition cond1 = new QuantCondition("cond1");
	private final QuantCondition cond2 = new QuantCondition("cond2");
	// private static final String SETUP_PROPERTIES = "setup.properties";
	private Map<String, Set<NWResult>> peptideAlignmentsBySequence;
	private final File setupPropertiesFile;
	private QuantParser parser;
	private HashSet<String> replicateNames;
	private final ProteinClusterQuantParameters params;
	private Map<String, Entry> annotatedProteins;

	public ProtClusterQuant(File setupPropertiesFile) {
		this.setupPropertiesFile = setupPropertiesFile;
		params = ProteinClusterQuantParameters.getInstance();
	}

	public static void main(String args[]) throws IOException {
		String propertiesFilePath = System.getProperty("user.dir") + File.separator + SETUP_PROPERTIES;
		if (args.length > 0) {
			propertiesFilePath = args[0];

		}
		log.info("Using setup.properties file at: " + propertiesFilePath);
		final File setupPropertiesFile = new File(propertiesFilePath);
		PropertiesReader.readProperties(setupPropertiesFile);
		ProtClusterQuant clusterCreator = new ProtClusterQuant(setupPropertiesFile);
		clusterCreator.run();
		System.exit(0);
	}

	public void run() throws IOException {

		Set<ProteinCluster> clusterSet = new HashSet<ProteinCluster>();

		try {
			List<Map<QuantCondition, QuantificationLabel>> labelsByConditionsList = getLabelsByconditionsList();

			parser = Utils.getQuantParser(params, labelsByConditionsList, params.getInputFileNamesArray());

			log.info("Reading input files...");
			Map<String, QuantifiedPeptideInterface> pepMap = parser.getPeptideMap();

			// List to hold all peptides
			List<QuantifiedPeptideInterface> peptideList = new ArrayList<QuantifiedPeptideInterface>();
			peptideList.addAll(parser.getPeptideMap().values());

			if (params.isMakeAlignments()) {
				makePeptideAlignments(peptideList);
			}
			int numClusters = 0;
			int i = 0;
			while (i < 1) {
				clusterSet = createClusters(pepMap);

				// filtering clusters
				List<PCQFilter> filters = params.getFilters();
				if (!filters.isEmpty()) {
					log.info("Filtering " + clusterSet.size() + " clusters");
					int numFilteredClusters = 0;
					final Iterator<ProteinCluster> iterator = clusterSet.iterator();
					while (iterator.hasNext()) {
						final ProteinCluster cluster = iterator.next();
						for (PCQFilter pcqFilter : filters) {
							pcqFilter.filter(cluster);
						}
						// check if the cluster still contains peptide nodes.
						if (cluster.getPeptideNodes().isEmpty()) {
							iterator.remove();
							numFilteredClusters++;
						}
					}
					log.info(numFilteredClusters + " clusters were discarded");
				}

				if (!filters.isEmpty() && numClusters != clusterSet.size()) {
					i = 0;
				} else {
					i++;
				}
				numClusters = clusterSet.size();
				// get peptide map
				pepMap = Utils.getPeptideMapFromClusters(clusterSet);
			}
			log.info("Final number of clusters after iterations: \t" + clusterSet.size());

			log.info("Identifying protein pairs in " + clusterSet.size() + " clusters...");
			int numProteinPairs = 0;
			int loopIndex = 0;
			int previousPercent = 0;
			for (ProteinCluster cluster : clusterSet) {
				int percent = (int) ((loopIndex++ * 100.0f) / clusterSet.size());
				if (previousPercent != percent) {
					previousPercent = percent;
					log.info(percent + "% clusters processed (" + loopIndex + "/" + clusterSet.size() + ")");
				}

				// create protein pairs in the cluster
				cluster.createPairs(peptideAlignmentsBySequence);
				// count protein pairs
				numProteinPairs += cluster.getProteinPairs().size();
			}
			log.info(numProteinPairs + " protein pairs identified.");

			Map<String, SanxotQuantResult> ratioStatsByPeptideNodeKey = null;
			if (params.isPerformRatioIntegration()) {
				// calculating consensus ratios up to peptide_exp_rep level
				SanXotAnalysisResult peptideRepSanxotResult = null;
				peptideRepSanxotResult = calculatePeptideReplicateRatios();

				// // if we have filtered, then remove
				// if (!params.getFilters().isEmpty())

				// make a custom sanxot analysis from peptide_rep to
				// peptide_node_rep
				final SanXotAnalysisResult peptideNodeRepSanxotResult = calculatePeptideNodeReplicateRatios(
						peptideRepSanxotResult, clusterSet);

				// make a custom sanxot analysis from peptide_node_rep to
				// peptide_node
				final SanXotAnalysisResult peptideNodeSanxotResult = calculatePeptideNodeRatios(
						peptideNodeRepSanxotResult, clusterSet);
				ratioStatsByPeptideNodeKey = peptideNodeSanxotResult.getLastIntegrationResults().getOutStatsRatios();

				// set the calculated values as final peptide node ratio values
				setIntegrationResultsIntoPeptideNodes(clusterSet, ratioStatsByPeptideNodeKey);
			}
			// print statistics and output files
			classiffyAndPrintStatistics(clusterSet);

			// print integration file
			if (ratioStatsByPeptideNodeKey != null) {
				printIntegrationFinalFile(clusterSet, ratioStatsByPeptideNodeKey);
			}

			// print PSEA QUANT files
			printPSEAQuantFiles(clusterSet);

			// export to XGMML format
			exportToXGMML(clusterSet);

			// rename TEMP output folder to output folder
			moveResultsToFinalFolder();

			log.info("DONE.");
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			log.error(e.getMessage());
		} catch (InterruptedException e) {
			e.printStackTrace();
			log.error(e.getMessage());
		} catch (ExecutionException e) {
			e.printStackTrace();
			log.error(e.getMessage());
		}

	}

	/**
	 * Writes a file with a line per each peptide node, sorted by FDR.
	 *
	 * @param clusterSet
	 * @param ratioStatsByPeptideNodeKey
	 */
	private void printIntegrationFinalFile(Set<ProteinCluster> clusterSet,
			Map<String, SanxotQuantResult> ratioStatsByPeptideNodeKey) {

		Map<String, ProteinCluster> proteinClustersByNodeID = new HashMap<String, ProteinCluster>();
		for (ProteinCluster cluster : clusterSet) {
			final Set<String> peptideNodeKeys = cluster.getPeptideNodeKeys();
			for (String peptideNodeID : peptideNodeKeys) {
				proteinClustersByNodeID.put(peptideNodeID, cluster);
			}
		}

		Map<String, Entry> annotatedProteins = getAnnotatedProteins();

		FileWriter outputIntegrationFinalFile = null;
		log.info("Printing Integration final file");
		try {
			final ProteinClusterQuantParameters params = ProteinClusterQuantParameters.getInstance();
			final File outputFileFolder = params.getTemporalOutputFolder();
			final String outputPrefix = params.getOutputPrefix();
			final String outputSuffix = params.getOutputSuffix();
			outputIntegrationFinalFile = new FileWriter(outputFileFolder.getAbsolutePath() + File.separator
					+ outputPrefix + "_FDR_" + outputSuffix + ".txt");

			outputIntegrationFinalFile.write(getPeptideNodeHeaderLine() + "\n");
			// sort peptide Nodes by FDR
			List<String> peptideNodeIDsSortedByFDR = getPeptideNodeIDsSortedByFDR(ratioStatsByPeptideNodeKey,
					proteinClustersByNodeID.keySet());
			for (String peptideNodeID : peptideNodeIDsSortedByFDR) {
				ProteinCluster cluster = proteinClustersByNodeID.get(peptideNodeID);
				final String geneNameString = Utils.getGeneNameString(annotatedProteins, cluster, null,
						params.isPrintOnlyFirstGene());

				String peptideNodeLine = getPeptideNodeLine(peptideNodeID, cluster, geneNameString,
						ratioStatsByPeptideNodeKey.get(peptideNodeID));

				outputIntegrationFinalFile.write(peptideNodeLine + "\n");
				outputIntegrationFinalFile.flush();
			}

			log.info("FDR file writen");

		} catch (IOException e) {
			log.error(e.getMessage());
		} finally {
			try {
				if (outputIntegrationFinalFile != null)
					outputIntegrationFinalFile.close();
			} catch (IOException e) {
				log.error(e.getMessage());
			}
		}
	}

	private String getPeptideNodeHeaderLine() {
		StringBuilder sb = new StringBuilder();

		// peptide Node ID
		sb.append("peptideNodeID");
		sb.append(sep);
		// protein cluster ID
		sb.append("clusterID");
		sb.append(sep);
		// geneNameString
		sb.append("genes");
		sb.append(sep);
		// unique node
		sb.append("uniquePeptideNode");
		sb.append(sep);
		// num peptide sequences
		sb.append("numPeptideSequences");
		sb.append(sep);
		// num psms
		sb.append("numPSMs");
		sb.append(sep);
		// log2 ratio
		sb.append("log2Ratio");
		sb.append(sep);
		// zRatio
		sb.append("normZRatio");
		sb.append(sep);
		// FDR
		sb.append("FDR");
		sb.append(sep);
		// weight
		sb.append("weight");
		sb.append(sep);
		// variance
		sb.append("variance");
		sb.append(sep);

		return sb.toString();
	}

	private String getPeptideNodeLine(String peptideNodeID, ProteinCluster cluster, String geneNameString,
			SanxotQuantResult sanxotQuantResult) {
		StringBuilder sb = new StringBuilder();

		PCQPeptideNode peptideNode = null;
		final Set<PCQPeptideNode> peptideNodes = cluster.getPeptideNodes();
		for (PCQPeptideNode pcqPeptideNode : peptideNodes) {
			if (pcqPeptideNode.getKey().equals(peptideNodeID)) {
				peptideNode = pcqPeptideNode;
			}
		}
		// peptide Node ID
		sb.append(peptideNodeID);
		sb.append(sep);
		// protein cluster ID
		sb.append(cluster.getClusterID());
		sb.append(sep);
		// geneNameString
		sb.append(geneNameString);
		sb.append(sep);
		// unique node
		sb.append(peptideNode.getPCQProteinNodes().size() == 1);
		sb.append(sep);
		// num peptide sequences
		sb.append(peptideNode.getIndividualPeptides().size());
		sb.append(sep);
		// num psms
		sb.append(peptideNode.getQuantifiedPSMs().size());
		sb.append(sep);
		// log2 ratio
		if (sanxotQuantResult != null) {
			sb.append(sanxotQuantResult.getLog2ratio());
		}
		sb.append(sep);
		// zRatio
		if (sanxotQuantResult != null) {
			sb.append(sanxotQuantResult.getzValue());
		}
		sb.append(sep);
		// FDR
		if (sanxotQuantResult != null) {
			sb.append(sanxotQuantResult.getFdr());
		}
		sb.append(sep);
		// weight
		if (sanxotQuantResult != null) {
			sb.append(sanxotQuantResult.getWeight());
		}
		sb.append(sep);
		// variance
		if (sanxotQuantResult != null) {
			sb.append(1 / sanxotQuantResult.getWeight());
		}
		sb.append(sep);
		return sb.toString();
	}

	private List<String> getPeptideNodeIDsSortedByFDR(final Map<String, SanxotQuantResult> ratioStatsByPeptideNodeKey,
			Set<String> peptideNodeIDs) {
		List<String> ret = new ArrayList<String>();
		ret.addAll(peptideNodeIDs);
		Collections.sort(ret, new Comparator<String>() {

			@Override
			public int compare(String peptideNodeID1, String peptideNodeID2) {
				double fdr1 = 1;
				double fdr2 = 1;
				if (ratioStatsByPeptideNodeKey.containsKey(peptideNodeID1)) {
					fdr1 = ratioStatsByPeptideNodeKey.get(peptideNodeID1).getFdr();
				}
				if (ratioStatsByPeptideNodeKey.containsKey(peptideNodeID2)) {
					fdr2 = ratioStatsByPeptideNodeKey.get(peptideNodeID2).getFdr();
				}
				return Double.compare(fdr1, fdr2);
			}

		});
		return ret;
	}

	private void setIntegrationResultsIntoPeptideNodes(Set<ProteinCluster> clusterSet,
			Map<String, SanxotQuantResult> ratioStatsByPeptideNodeKey) {
		log.info("Assigning calculated peptide node ratios and FDRs to network...");
		int numNodesWithConsensusRatio = 0;
		for (ProteinCluster proteinCluster : clusterSet) {
			final Set<PCQPeptideNode> peptideNodes = proteinCluster.getPeptideNodes();
			for (PCQPeptideNode pcqPeptideNode : peptideNodes) {
				final String peptideNodeKey = pcqPeptideNode.getKey();
				if (ratioStatsByPeptideNodeKey.containsKey(peptideNodeKey)) {
					final SanxotQuantResult sanxotQuantResult = ratioStatsByPeptideNodeKey.get(peptideNodeKey);

					final Map<QuantCondition, QuantificationLabel> labelsByConditions = getConsensusQuantificationLabelsByConditions();
					final Map<QuantificationLabel, QuantCondition> conditionsByLabels = swapMap(labelsByConditions);
					CensusRatio ratio = new CensusRatio(sanxotQuantResult.getNonLog2ratio(), false, conditionsByLabels,
							labelsByConditions.get(cond1), labelsByConditions.get(cond2), AggregationLevel.PROTEINGROUP,
							"Peptide node ratio");
					ratio.setCombinationType(CombinationType.WEIGHTED_AVERAGE);
					RatioScore fdrScore = new RatioScore(String.valueOf(sanxotQuantResult.getFdr()), "FDR",
							"PSM-level quantification confidence metric", "FDR of significantly changing ratio");
					ratio.setRatioScore(fdrScore);
					pcqPeptideNode.addConsensusRatio(ratio);
					// weight as confidence value
					pcqPeptideNode.setConfidenceValue(sanxotQuantResult.getWeight());
					numNodesWithConsensusRatio++;
				} else {
					// there are peptide with no ratios (singletons) and it is
					// fine
				}
			}
		}
		log.info("Information assigned to network successfully for " + numNodesWithConsensusRatio + " peptide nodes");
	}

	/**
	 * This is just used in the calculated ratios at peptide node level
	 *
	 * @return
	 */
	private Map<QuantCondition, QuantificationLabel> getConsensusQuantificationLabelsByConditions() {
		for (ExperimentFiles experimentFiles : ProteinClusterQuantParameters.getInstance().getInputFileNames()) {
			for (String replicateFileName : experimentFiles.getRelicateFileNames()) {
				final Map<QuantCondition, QuantificationLabel> labelsByConditions = getLabelsByConditions(
						replicateFileName);
				if (labelsByConditions != null) {
					return labelsByConditions;
				}
			}
		}
		return null;
	}

	private Map<QuantificationLabel, QuantCondition> swapMap(Map<QuantCondition, QuantificationLabel> map) {
		Map<QuantificationLabel, QuantCondition> ret = new HashMap<QuantificationLabel, QuantCondition>();
		for (QuantCondition condition : map.keySet()) {
			final QuantificationLabel quantificationLabel = map.get(condition);
			ret.put(quantificationLabel, condition);
		}
		return ret;
	}

	private SanXotAnalysisResult calculatePeptideNodeRatios(SanXotAnalysisResult peptideNodeRepSanxotResult,
			Set<ProteinCluster> clusterSet) throws IOException, InterruptedException, ExecutionException {
		final HashMap<String, IntegrationResultWrapper> experimentIntegrationResults = peptideNodeRepSanxotResult
				.getExperimentIntegrationResults();
		File workingFolder = QuantAnalysis.createWorkingFolder(params.getTemporalOutputFolder(),
				ANALYSIS_LEVEL_OUTCOME.PEPTIDE);
		File relationshipFile = writeRelationshipFileFromPeptideExpNodeToPeptideNode(clusterSet);
		List<File> files = new ArrayList<File>();
		for (String experimentName : experimentIntegrationResults.keySet()) {
			files.add(experimentIntegrationResults.get(experimentName).getHigherLevelDataFile());
		}
		SanXotAnalysisResult ret = new SanXotAnalysisResult(null);
		IntegrationResultWrapper integrationResult = null;
		if (files.size() == 1) {
			// if there is only one file, use the file from previous integration
			integrationResult = peptideNodeRepSanxotResult.getLastIntegrationResults();

		} else {
			File mergedFile = new File(
					workingFolder.getAbsolutePath() + File.separator + "PeptideNodeExp2PeptideNode.xls");
			edu.scripps.yates.utilities.files.FileUtils.mergeFiles(files, mergedFile, true);
			integrationResult = SanxotRunner.integrate(relationshipFile, mergedFile, null, "PeptideNodeExp2PeptideNode",
					null, workingFolder, true);

			if (isRemoveOutliers()) {
				File infoFile = integrationResult.getInfoFile();
				String prefix = "outliers_removed";
				final OutlierRemovalResultWrapper removeOutliers = SanxotRunner.removeOutliers(relationshipFile,
						mergedFile, infoFile, prefix, getParams().getOutliersRemovalFDR(), workingFolder);
				relationshipFile = removeOutliers.getRelatFile();
				integrationResult = SanxotRunner.integrate(relationshipFile, mergedFile, null,
						"PeptideNodeExp2PeptideNode", null, workingFolder, false);
			}
			ret.addIntegrationResult(integrationResult);
		}
		// make the last integration in order to get the FDR
		IntegrationResultWrapper finalIntegrationResult = SanxotRunner.integrate(null,
				integrationResult.getHigherLevelDataFile(), null, "PeptideNode2All", null, workingFolder, true);
		ret.addIntegrationResult(finalIntegrationResult);
		return ret;
	}

	private boolean isRemoveOutliers() {
		return getParams().getOutliersRemovalFDR() != null;
	}

	private File writeRelationshipFileFromPeptideExpNodeToPeptideNode(Set<ProteinCluster> clusterSet)
			throws IOException {
		final File workingFolder = QuantAnalysis.createWorkingFolder(params.getTemporalOutputFolder(),
				ANALYSIS_LEVEL_OUTCOME.PEPTIDE);
		File outputFile = new File(
				workingFolder.getAbsolutePath() + File.separator + "Relat_PeptideExpNode2PeptideNodeRep.csv");
		// Set<String> replicateNames =
		// getReplicateNamesFromClusters(clusterSet);
		FileWriter fw = new FileWriter(outputFile);
		BufferedWriter bw = new BufferedWriter(fw);
		try {
			fw.write("#Peptide_node\tPeptide_node_exp\n");

			final Map<String, List<String>> replicateNamesByExperimentNameMap = params
					.getReplicateNamesByExperimentNameMap();

			for (ProteinCluster proteinCluster : clusterSet) {
				for (PCQPeptideNode peptideNode : proteinCluster.getPeptideNodes()) {
					final String peptideNodeKey = peptideNode.getKey();
					String experimentKey = "";
					for (String experimentName : replicateNamesByExperimentNameMap.keySet()) {
						if (replicateNamesByExperimentNameMap.size() > 1) {
							experimentKey = experimentName;
						}
						// peptide_node_rep

						// per each peptide in a replicate in the cluster,
						// write a line like peptideNode1_rep1 \t
						// peptide1_rep1

						// peptide_rep

						StringBuilder sb = new StringBuilder();
						sb.append(peptideNodeKey).append("\t");
						sb.append(peptideNodeKey);
						if (!"".equals(experimentKey)) {
							sb.append("_").append(experimentKey);
						}
						sb.append("\n");
						bw.write(sb.toString());

					}
				}
			}
		} catch (Exception e) {
			if (bw != null) {
				bw.close();
			}
			throw e;
		} finally {
			if (bw != null) {
				bw.close();
			}
		}
		return outputFile;
	}

	/**
	 * Make a custom sanxot analysis from peptide_rep to peptide_node_rep, that
	 * is from peptides in each replicate to peptide nodes in each replicate.
	 *
	 * @param peptideRepSanxotResult
	 * @param clusterSet
	 * @return
	 * @throws ExecutionException
	 * @throws InterruptedException
	 */
	private SanXotAnalysisResult calculatePeptideNodeReplicateRatios(SanXotAnalysisResult peptideRepSanxotResult,
			Set<ProteinCluster> clusterSet) throws IOException, InterruptedException, ExecutionException {
		File relationshipFile = writeRelationshipFileFromPeptideExpRepToPeptideNodeExpRep(clusterSet);

		File workingFolder = QuantAnalysis.createWorkingFolder(params.getTemporalOutputFolder(),
				ANALYSIS_LEVEL_OUTCOME.PEPTIDE);
		SanXotAnalysisResult ret = new SanXotAnalysisResult(null);
		final Map<String, List<String>> replicateNamesByExperimentNameMap = params
				.getReplicateNamesByExperimentNameMap();

		String experimentKey = "";
		for (String experimentName : replicateNamesByExperimentNameMap.keySet()) {
			if (replicateNamesByExperimentNameMap.size() > 1) {
				experimentKey = experimentName;
			}
			Map<String, IntegrationResultWrapper> replicateIntegrations = peptideRepSanxotResult
					.getReplicateIntegrationResultsByExperiment().get(experimentName);
			List<File> dataFiles = new ArrayList<File>();
			String replicateKey = "";
			IntegrationResultWrapper replicateIntegrationResult = null;
			for (String replicateName : replicateIntegrations.keySet()) {
				if (replicateIntegrations.size() > 1) {
					replicateKey = replicateName;
				}
				File infoFile = null;
				final IntegrationResultWrapper integrationResultForReplicate = replicateIntegrations.get(replicateName);
				final String prefix = "PeptideExpRep2PeptideNodeExpRep_" + experimentKey + replicateKey;
				boolean checkRelationshipValidity = params.getFilters().isEmpty();
				replicateIntegrationResult = SanxotRunner.integrate(relationshipFile,
						integrationResultForReplicate.getHigherLevelDataFile(), infoFile, prefix, null, workingFolder,
						checkRelationshipValidity);

				if (isRemoveOutliers()) {
					infoFile = replicateIntegrationResult.getInfoFile();
					String outliersPrefix = "outliers_removed_" + prefix;
					final OutlierRemovalResultWrapper removeOutliers = SanxotRunner.removeOutliers(relationshipFile,
							integrationResultForReplicate.getHigherLevelDataFile(), infoFile, outliersPrefix,
							getParams().getOutliersRemovalFDR(), workingFolder);
					relationshipFile = removeOutliers.getRelatFile();
					replicateIntegrationResult = SanxotRunner.integrate(relationshipFile,
							integrationResultForReplicate.getHigherLevelDataFile(), null, prefix, null, workingFolder,
							false);
				}

				dataFiles.add(replicateIntegrationResult.getHigherLevelDataFile());
				ret.addReplicateExperimentIntegrationResult(replicateIntegrationResult, experimentName, replicateName);

			}
			if (dataFiles.size() > 1) {
				File relationshipFile2 = writeRelationshipFileFromPeptideNodeExpRepToPeptideNodeExp(clusterSet);
				File mergedFile = new File(workingFolder.getAbsolutePath() + File.separator
						+ "PeptideRep2PeptideNodeExpRep_" + experimentKey + ".xls");
				// append the results files
				edu.scripps.yates.utilities.files.FileUtils.mergeFiles(dataFiles, mergedFile, true);
				String prefixExperiment = "PeptideNodeExpRep2PeptideNodeExp_" + experimentKey;
				IntegrationResultWrapper experimentIntegrationResult = SanxotRunner.integrate(relationshipFile2,
						mergedFile, null, prefixExperiment, null, workingFolder, true);
				if (isRemoveOutliers()) {
					File infoFile = experimentIntegrationResult.getInfoFile();
					String outliersPrefix = "outliers_removed_" + prefixExperiment;
					final OutlierRemovalResultWrapper removeOutliers = SanxotRunner.removeOutliers(relationshipFile2,
							mergedFile, infoFile, outliersPrefix, getParams().getOutliersRemovalFDR(), workingFolder);
					relationshipFile2 = removeOutliers.getRelatFile();
					experimentIntegrationResult = SanxotRunner.integrate(relationshipFile2, mergedFile, null,
							prefixExperiment, null, workingFolder, false);
				}
				ret.addExperimentIntegrationResult(experimentIntegrationResult, experimentName);
			} else {
				ret.addExperimentIntegrationResult(replicateIntegrationResult, experimentName);
			}

		}

		return ret;
	}

	private File writeRelationshipFileFromPeptideNodeExpRepToPeptideNodeExp(Set<ProteinCluster> clusterSet)
			throws IOException {

		final File workingFolder = QuantAnalysis.createWorkingFolder(params.getTemporalOutputFolder(),
				ANALYSIS_LEVEL_OUTCOME.PEPTIDE);
		File outputFile = new File(
				workingFolder.getAbsolutePath() + File.separator + "Relat_PeptideNodeExpRep2PeptideNodeExp.csv");
		// Set<String> replicateNames =
		// getReplicateNamesFromClusters(clusterSet);
		FileWriter fw = new FileWriter(outputFile);
		BufferedWriter bw = new BufferedWriter(fw);
		try {
			fw.write("#Peptide_node_exp\tPeptide_node_exp_rep\n");
			final Map<String, List<String>> replicateNamesByExperimentNameMap = params
					.getReplicateNamesByExperimentNameMap();
			for (ProteinCluster proteinCluster : clusterSet) {
				for (PCQPeptideNode peptideNode : proteinCluster.getPeptideNodes()) {
					final String peptideNodeKey = peptideNode.getKey();
					String experimentKey = "";
					for (String experimentName : replicateNamesByExperimentNameMap.keySet()) {
						if (replicateNamesByExperimentNameMap.size() > 1) {
							experimentKey = experimentName;
						}
						// peptide_node_rep
						StringBuilder peptideNodeExp = new StringBuilder();
						peptideNodeExp.append(peptideNodeKey);
						if (!"".equals(experimentKey)) {
							peptideNodeExp.append("_").append(experimentKey);
						}
						String replicateKey = "";
						for (String replicateName : replicateNamesByExperimentNameMap.get(experimentName)) {
							if (replicateNamesByExperimentNameMap.get(experimentName).size() > 1) {
								replicateKey = replicateName;
							}
							// per each peptide in a replicate in the cluster,
							// write a line like peptideNode1_rep1 \t
							// peptide1_rep1

							// peptide_rep

							StringBuilder sb = new StringBuilder();
							sb.append(peptideNodeExp).append("\t");
							sb.append(peptideNodeKey);
							if (!"".equals(replicateKey)) {
								sb.append("_").append(replicateKey);
							}
							if (!"".equals(experimentKey)) {
								sb.append("_").append(experimentKey);
							}
							sb.append("\n");
							bw.write(sb.toString());

						}
					}
				}
			}
		} catch (Exception e) {
			if (bw != null) {
				bw.close();
			}
			throw e;
		} finally {
			if (bw != null) {
				bw.close();
			}
		}
		return outputFile;
	}

	/**
	 * peptide_rep to peptide_node_rep<br>
	 * Example in a peptide node "peptideNode1" where there are two peptides
	 * (peptide1, peptide2) in two replicates (rep1, rep2):<br>
	 * peptideNode1_rep1 peptide1_rep1<br>
	 * peptideNode1_rep1 peptide2_rep1<br>
	 * peptideNode1_rep2 peptide1_rep2<br>
	 * peptideNode1_rep2 peptide2_rep2<br>
	 * peptideNode1_rep2 peptide3_rep2<br>
	 *
	 * @param clusterSet
	 * @return
	 * @throws IOException
	 */
	private File writeRelationshipFileFromPeptideExpRepToPeptideNodeExpRep(Set<ProteinCluster> clusterSet)
			throws IOException {
		final File workingFolder = QuantAnalysis.createWorkingFolder(params.getTemporalOutputFolder(),
				ANALYSIS_LEVEL_OUTCOME.PEPTIDE);
		File outputFile = new File(
				workingFolder.getAbsolutePath() + File.separator + "Relat_PeptideExpRep2PeptideNodeExpRep.csv");
		// Set<String> replicateNames =
		// getReplicateNamesFromClusters(clusterSet);
		FileWriter fw = new FileWriter(outputFile);
		BufferedWriter bw = new BufferedWriter(fw);
		try {
			fw.write("#Peptide_node_exp_rep\tPeptide_exp_rep\n");
			final Map<String, List<String>> replicateNamesByExperimentNameMap = params
					.getReplicateNamesByExperimentNameMap();

			for (ProteinCluster proteinCluster : clusterSet) {
				for (PCQPeptideNode peptideNode : proteinCluster.getPeptideNodes()) {
					final String peptideNodeKey = peptideNode.getKey();
					String experimentKey = "";
					for (String experimentName : replicateNamesByExperimentNameMap.keySet()) {
						if (replicateNamesByExperimentNameMap.size() > 1) {
							experimentKey = experimentName;
						}
						String replicateKey = "";
						for (String replicateName : replicateNamesByExperimentNameMap.get(experimentName)) {
							if (replicateNamesByExperimentNameMap.get(experimentName).size() > 1) {
								replicateKey = replicateName;
							}
							// per each peptide in a replicate in the cluster,
							// write
							// a line like peptideNode1_rep1 \t peptide1_rep1

							// peptide_node_rep

							// peptide_rep
							Set<QuantifiedPeptideInterface> peptidesInReplicate = peptideNode
									.getQuantifiedPeptidesInReplicate(replicateName);
							for (QuantifiedPeptideInterface quantifiedPeptideInterface : peptidesInReplicate) {
								if (quantifiedPeptideInterface.getSequence().contains("AALCAVHVIR")) {
									log.info(quantifiedPeptideInterface);
								}
								// a line per peptide
								bw.write(peptideNodeKey);
								if (!"".equals(replicateKey)) {
									bw.write("_");
									bw.write(replicateKey);
								}
								if (!"".equals(experimentKey)) {
									bw.write("_");
									bw.write(experimentKey);
								}
								bw.write("\t");
								bw.write(quantifiedPeptideInterface.getSequence());
								if (!"".equals(replicateKey)) {
									bw.write("_");
									bw.write(replicateKey);
								}
								if (!"".equals(experimentKey)) {
									bw.write("_");
									bw.write(experimentKey);
								}
								bw.write("\n");

							}
						}
					}
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			if (bw != null) {
				bw.close();
			}
			throw e;
		} finally {
			if (bw != null) {
				bw.close();
			}
		}
		return outputFile;
	}

	private Set<String> getReplicateNamesFromClusters(Set<ProteinCluster> clusterSet) {
		if (replicateNames == null) {
			replicateNames = new HashSet<String>();
			for (ProteinCluster cluster : clusterSet) {
				final Set<QuantifiedPeptideInterface> peptideSet = cluster.getPeptideSet();
				replicateNames.addAll(getReplicateNamesFromPeptides(peptideSet));

			}
		}
		return replicateNames;
	}

	/**
	 * make a custom sanxot analysis from peptide_node_rep to peptide_node
	 *
	 * @return
	 * @throws IOException
	 */
	private SanXotAnalysisResult calculatePeptideReplicateRatios() throws IOException {
		log.info("Running SanXot algorithm for integrating quantitative ratios");
		long t1 = System.currentTimeMillis();
		SanxotRunner sanxotRunner = new SanxotRunner(this, getQuantType(), getParams().getTemporalOutputFolder(), cond1,
				cond2, getParams().getFastaFile(), getParams().getOutliersRemovalFDR());
		SanXotAnalysisResult result = sanxotRunner.run();
		final String time = DatesUtil.getDescriptiveTimeFromMillisecs(System.currentTimeMillis() - t1);
		log.info("Sanxot run finished in " + time);
		return result;

	}

	private QuantificationType getQuantType() {
		final InputType inputType = getParams().getInputType();
		if (inputType == InputType.CENSUS_CHRO) {
			return QuantificationType.ISOTOPOLOGUES;
		} else if (inputType == InputType.CENSUS_OUT) {
			return QuantificationType.SILAC;
		} else
			// if (inputType == InputType.SEPARATED_VALUES) {
			return QuantificationType.UNKNOWN;
		// }
	}

	private void moveResultsToFinalFolder() throws IOException {
		log.info("Moving results to final output folder...");
		FileUtils.copyDirectory(getParams().getTemporalOutputFolder(), getParams().getOutputFileFolder());
		log.info("Deleting TEMP folder...");
		FileUtils.deleteDirectory(getParams().getTemporalOutputFolder());
		log.info("TEMP folder deleted.");
	}

	private void exportToXGMML(Set<ProteinCluster> clusterSet) {
		Map<String, Entry> annotatedProteins = getAnnotatedProteins();
		XgmmlExporter exporter = new XgmmlExporter();
		exporter.exportToXGMMLUsingNodes(clusterSet, annotatedProteins, cond1, cond2);
	}

	private void makePeptideAlignments(List<QuantifiedPeptideInterface> peptideList) throws IOException {
		ProteinClusterQuantParameters params = ProteinClusterQuantParameters.getInstance();
		final File outputFileFolder = params.getTemporalOutputFolder();
		final String outputPrefix = params.getOutputPrefix();
		final String outputSuffix = params.getOutputSuffix();
		log.info("Aligning peptides...");
		FileWriter output = null;
		try {
			output = new FileWriter(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix + "_pep_"
					+ outputSuffix + ".txt");
			long t = System.currentTimeMillis();
			// store the peptide alignments
			peptideAlignmentsBySequence = Utils.alignPeptides(peptideList, cond1, cond2, output);
			log.debug("Alignments done in " + Utils.format((t - System.currentTimeMillis()) * 1.0 / 1000.0) + "sg.");
		} finally {
			if (output != null)
				output.close();
		}

	}

	public Map<QuantCondition, QuantificationLabel> getLabelsByConditions(String replicateName) {

		List<Map<QuantCondition, QuantificationLabel>> labelsByConditionsList = getLabelsByconditionsList();
		final List<ExperimentFiles> inputFileNames = ProteinClusterQuantParameters.getInstance().getInputFileNames();
		int j = 0;
		for (ExperimentFiles experimentFiles : inputFileNames) {
			for (String replicateFileName : experimentFiles.getRelicateFileNames()) {
				if (replicateFileName.equals(replicateName)) {
					return labelsByConditionsList.get(j);
				}
				j++;
			}
		}
		return null;
	}

	private List<Map<QuantCondition, QuantificationLabel>> getLabelsByconditionsList() {
		List<Map<QuantCondition, QuantificationLabel>> labelsByConditionsList = new ArrayList<Map<QuantCondition, QuantificationLabel>>();

		// per L/H xml file
		Map<QuantCondition, QuantificationLabel> labelsByConditions = new HashMap<QuantCondition, QuantificationLabel>();
		// TODO
		labelsByConditions.put(cond1, QuantificationLabel.LIGHT);
		labelsByConditions.put(cond2, QuantificationLabel.HEAVY);

		int max = getParams().getInputFileNamesArray().length;
		if (getParams().getExperimentNames().size() < 2 && getParams().isLabelSwap()) {
			throw new IllegalArgumentException(
					"Having only one experiment is not compatible with label swap option enabled");
		}
		if (getParams().isLabelSwap()) {
			max = max / 2;
		}

		for (int i = 0; i < max; i++) {
			labelsByConditionsList.add(labelsByConditions);
		}
		// labelsByConditionsList.add(labelsByConditions);
		// labelsByConditionsList.add(labelsByConditions);
		// labelsByConditionsList.add(labelsByConditions);

		// per H/L xml file
		if (getParams().isLabelSwap()) {
			Map<QuantCondition, QuantificationLabel> labelsByConditions2 = new HashMap<QuantCondition, QuantificationLabel>();
			labelsByConditions2.put(cond1, QuantificationLabel.HEAVY);
			labelsByConditions2.put(cond2, QuantificationLabel.LIGHT);
			for (int i = 0; i < max; i++) {
				labelsByConditionsList.add(labelsByConditions2);
			}
		}
		return labelsByConditionsList;
	}

	private ProteinClusterQuantParameters getParams() {
		return ProteinClusterQuantParameters.getInstance();
	}

	private Set<ProteinCluster> createClusters(Map<String, QuantifiedPeptideInterface> peptideMap) {
		Map<String, ProteinCluster> clustersByPeptideSequence = new HashMap<String, ProteinCluster>();
		Set<ProteinCluster> clusterSet = new HashSet<ProteinCluster>();
		log.info("Starting clustering " + peptideMap.size() + " peptides...");
		long t0 = System.currentTimeMillis();
		for (QuantifiedPeptideInterface peptide1 : peptideMap.values()) {
			String sequence1 = peptide1.getSequence();

			ProteinCluster cluster = null;

			// if peptide has a cluster associated
			if (clustersByPeptideSequence.containsKey(sequence1)) {
				// grab it
				cluster = clustersByPeptideSequence.get(sequence1);
			} else {
				cluster = new ProteinCluster();
				// adds cluster to cluster set
				clusterSet.add(cluster);
				// Map<String, protCluster> <- cluster (key = sequence)
				clustersByPeptideSequence.put(sequence1, cluster);
			}
			if (sequence1.equals("AALCAVHVIR")) {
				log.info(cluster.getPeptideSet());
			}
			// put peptide in cluster
			cluster.addPeptide(peptide1);
			// in case of having peptide alignments done
			if (peptideAlignmentsBySequence != null && peptideAlignmentsBySequence.containsKey(sequence1)) {
				Set<NWResult> gAS = peptideAlignmentsBySequence.get(sequence1);
				for (NWResult AlignResult : gAS) {
					if (AlignResult.getFinalAlignmentScore() >= ProteinClusterQuantParameters.getInstance()
							.getFinalAlignmentScore()
							&& AlignResult.getSequenceIdentity() >= ProteinClusterQuantParameters.getInstance()
									.getSequenceIdentity()
							&& AlignResult.getMaxConsecutiveIdenticalAlignment() >= ProteinClusterQuantParameters
									.getInstance().getMaxConsecutiveIdenticalAlignment()) {

						String seq2 = AlignResult.getSeq1();

						// since there are two sequences when comparing,
						// this makes sure we have both, not a two of
						// the same sequence
						if (seq2.equals(sequence1)) {
							seq2 = AlignResult.getSeq2();
						}

						QuantifiedPeptideInterface pep2 = peptideMap.get(seq2);

						// checking to see if peptide 2 is already in a
						// cluster
						if (pep2 != null) {
							if (clustersByPeptideSequence.containsKey(pep2.getSequence())) {
								ProteinCluster cluster2 = clustersByPeptideSequence.get(pep2.getSequence());
								if (!cluster.equals(cluster2)) {
									// merges the clusters with similar peptides
									cluster = Utils.mergeClusters(cluster, cluster2);
									clusterSet.remove(cluster2);
									for (QuantifiedPeptideInterface quantifiedPeptide : cluster.getPeptideSet()) {
										clustersByPeptideSequence.put(quantifiedPeptide.getSequence(), cluster);
									}
								}

								// Add (Quant pep2, hisProtein) to Cluster)
								cluster.addPeptide(pep2);
								clustersByPeptideSequence.put(pep2.getSequence(), cluster);
							}
							// add alignment to the cluster
							cluster.addAlignment(AlignResult);
						}
					}
				}
			}
			// get proteins of the peptide and add them to the cluster
			Set<QuantifiedProteinInterface> proteinSet = peptide1.getQuantifiedProteins();
			for (QuantifiedProteinInterface protein : proteinSet) {
				// put protein in cluster
				cluster.addProtein(protein);

				// peptide 2 <- protein
				for (QuantifiedPeptideInterface peptide2 : protein.getQuantifiedPeptides()) {
					// checking to see if peptide 2 is already in a
					// cluster
					if (clustersByPeptideSequence.containsKey(peptide2.getSequence())) {
						ProteinCluster cluster2 = clustersByPeptideSequence.get(peptide2.getSequence());
						if (!cluster.equals(cluster2)) {
							// merge the clusters
							cluster = Utils.mergeClusters(cluster, cluster2);
							clusterSet.remove(cluster2);
							for (QuantifiedPeptideInterface quantifiedPeptide : cluster.getPeptideSet()) {
								clustersByPeptideSequence.put(quantifiedPeptide.getSequence(), cluster);
							}
						}
					}

					// add new peptide 2 to cluster
					cluster.addPeptide(peptide2);

					// Map <- peptide, cluster
					clustersByPeptideSequence.put(peptide2.getSequence(), cluster);
				}
			}
		}
		double time = System.currentTimeMillis() - t0;

		log.info("Clustering done in " + DatesUtil.getDescriptiveTimeFromMillisecs(time));

		log.info("Collapsing clusters...");
		t0 = System.currentTimeMillis();
		for (ProteinCluster proteinCluster : clusterSet) {
			proteinCluster.reorganize();
		}
		time = System.currentTimeMillis() - t0;
		log.info("Cluster collapsing done in " + DatesUtil.getDescriptiveTimeFromMillisecs(time));

		return clusterSet;
	}

	private void printPSEAQuantFiles(Set<ProteinCluster> clusterSet) {

		Map<String, Entry> annotatedProteins = getAnnotatedProteins();

		if (getParams().getInputFileNamesArray() != null && !"".equals(getParams().getInputFileNamesArray()[0])) {
			log.info("Replicate names are empty. PSEA-Quant files will contain a column for each original RAW file");

			FileWriter outputPSEAQuant = null;
			FileWriter outputPSEAQuant_inv = null;
			log.info("Printing PSEA-Quant input files");
			try {
				final ProteinClusterQuantParameters params = ProteinClusterQuantParameters.getInstance();
				final File outputFileFolder = params.getTemporalOutputFolder();
				final String outputPrefix = params.getOutputPrefix();
				final String outputSuffix = params.getOutputSuffix();
				outputPSEAQuant = new FileWriter(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix
						+ "_pseaQuant_" + outputSuffix + ".txt");
				outputPSEAQuant_inv = new FileWriter(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix
						+ "_pseaQuant_inv_" + outputSuffix + ".txt");
				final Mean meanCalculator = new Mean();
				List<String> replicateNameList = new ArrayList<String>();
				Set<String> replicateNames = new HashSet<String>();
				for (ProteinCluster cluster : clusterSet) {
					final Set<QuantifiedPeptideInterface> peptideSet = cluster.getPeptideSet();
					replicateNames.addAll(getReplicateNamesFromPeptides(peptideSet));

				}

				replicateNameList.addAll(replicateNames);
				Collections.sort(replicateNameList);
				int numClustersSkippedNotHavingGeneName = 0;
				int numClustersSkippedNotHavingGoodRatios = 0;

				// do not include simulans and virilis in PSEA Quant output
				Set<String> validTaxonomies = new HashSet<String>();
				validTaxonomies.add("melanogaster");
				validTaxonomies.add("human");
				validTaxonomies.add("homo");
				validTaxonomies.add("musculus");
				validTaxonomies.add("cerevisiae");
				for (ProteinCluster cluster : clusterSet) {

					final String geneNameString = Utils.getGeneNameString(annotatedProteins, cluster, validTaxonomies,
							params.isPrintOnlyFirstGene());
					if (geneNameString.contains("Dsim")) {
						log.info(Utils.getGeneNameString(annotatedProteins, cluster, validTaxonomies,
								params.isPrintOnlyFirstGene()));
					}
					if ("".equals(geneNameString)) {
						numClustersSkippedNotHavingGeneName++;
						continue;
					}
					List<Double> ratioMeans = new ArrayList<Double>();
					// TODO
					// according to the manuscript, it is the average of the
					// ratios
					// of the peptide nodes.

					for (String replicateName : replicateNameList) {
						final Set<ProteinPair> proteinPairs = cluster.getProteinPairs();
						Set<Double> ratioValues = new HashSet<Double>();
						if (proteinPairs.isEmpty()) {
							final QuantRatio pepRatio = Utils.getConsensusRatio(
									cluster.getPeptideSet(params.isExcludeUniquePeptides()), cond1, cond2,
									replicateName);
							if (pepRatio != null) {
								final Double countRatio = pepRatio.getNonLogRatio(cond1, cond2);
								if (countRatio != null && !Double.isNaN(countRatio) && !Double.isInfinite(countRatio)) {
									ratioValues.add(countRatio);
								}
							}
						} else {
							for (ProteinPair proteinPair : proteinPairs) {
								final Map<String, Set<QuantifiedPeptideInterface>> sharedPeptidesMap_S12 = Utils
										.getSharedPeptidesMap(proteinPair.getProtein1(), proteinPair.getProtein2(),
												false);
								for (Set<QuantifiedPeptideInterface> peptidesS12 : sharedPeptidesMap_S12.values()) {
									List<Set<QuantifiedPeptideInterface>> pep12DoubleList = new ArrayList<Set<QuantifiedPeptideInterface>>();
									if (ProteinClusterQuantParameters.getInstance()
											.isCollapseIndistinguishablePeptides()) {
										pep12DoubleList.add(peptidesS12);
									} else {
										for (QuantifiedPeptideInterface peptide : peptidesS12) {
											Set<QuantifiedPeptideInterface> set = new HashSet<QuantifiedPeptideInterface>();
											set.add(peptide);
											pep12DoubleList.add(set);
										}
									}
									for (Set<QuantifiedPeptideInterface> sharedPeptides_S12 : pep12DoubleList) {
										final QuantRatio pepRatio = Utils.getConsensusRatio(sharedPeptides_S12, cond1,
												cond2, replicateName);
										if (pepRatio != null) {
											final Double countRatio = pepRatio.getNonLogRatio(cond1, cond2);
											if (countRatio != null && !Double.isNaN(countRatio)
													&& !Double.isInfinite(countRatio)) {
												ratioValues.add(countRatio);
											}
										}
									}
								}
							}
						}
						Double ratioMean = null;
						if (ratioValues.isEmpty()) {
							numClustersSkippedNotHavingGoodRatios++;
						} else {
							double[] ratios = new double[ratioValues.size()];
							int index = 0;
							for (double ratioValue : ratioValues) {
								ratios[index] = ratioValue;
							}
							ratioMean = meanCalculator.evaluate(ratios);
						}
						ratioMeans.add(ratioMean);
					}
					boolean validCluster = false;
					for (Double ratioMean : ratioMeans) {
						if (ratioMean != null) {
							validCluster = true;
							break;
						}
					}
					if (!validCluster) {
						numClustersSkippedNotHavingGoodRatios++;
						continue;
					}
					outputPSEAQuant.write(geneNameString + "\t");
					outputPSEAQuant_inv.write(geneNameString + "\t");
					for (Double ratioMean : ratioMeans) {
						if (ratioMean != null) {
							outputPSEAQuant.write(String.valueOf(ratioMean));
							outputPSEAQuant_inv.write(String.valueOf(1.0 / ratioMean));
						}
						outputPSEAQuant.write("\t");
						outputPSEAQuant_inv.write("\t");
					}
					outputPSEAQuant.write("\n");
					outputPSEAQuant.flush();
					outputPSEAQuant_inv.write("\n");
					outputPSEAQuant_inv.flush();

				}

				log.info("Clusters skipped due to not having appropiate gene name: "
						+ numClustersSkippedNotHavingGeneName);
				log.info("Clusters skipped due to not having at least one valid ratio: "
						+ numClustersSkippedNotHavingGoodRatios);
				log.info("PSEA-Quant files writen");

			} catch (

			IOException e)

			{
				log.error(e.getMessage());
			} finally

			{
				try {
					if (outputPSEAQuant != null)
						outputPSEAQuant.close();
					if (outputPSEAQuant_inv != null)
						outputPSEAQuant_inv.close();
				} catch (IOException e) {
					log.error(e.getMessage());
				}
			}
		} else {
			log.info("PSEA-Quant intput files will NOT be created");
		}
	}

	private Map<String, Entry> getAnnotatedProteins() {
		if (annotatedProteins == null) {
			if (getParams().getUniprotReleasesFolder() != null) {
				log.info("Getting UniprotKB annotations for " + parser.getProteinMap().size() + " proteins");
				UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(
						getParams().getUniprotReleasesFolder(), true);

				// exporting to xgmml
				annotatedProteins = uplr.getAnnotatedProteins(getParams().getUniprotVersion(),
						parser.getUniprotAccSet());
				log.info(annotatedProteins.size() + " annotations retrieved out of " + parser.getProteinMap().size()
						+ " proteins");
			}
		}
		return annotatedProteins;
	}

	private void classiffyAndPrintStatistics(Set<ProteinCluster> clusterSet) {

		FileWriter outputCases = null;
		FileWriter outputGroups = null;
		FileWriter outputProteins = null;
		FileWriter outputNames = null;
		FileWriter outputSummary = null;
		FileWriter outputPairs = null;
		File outputSummaryFile = null;
		Map<Classification1Case, Integer> classification1Counters = new HashMap<Classification1Case, Integer>();
		for (Classification1Case case1 : Classification1Case.values()) {
			classification1Counters.put(case1, 0);
		}
		Map<Classification2Case, Integer> classification2Counters = new HashMap<Classification2Case, Integer>();
		for (Classification2Case case2 : Classification2Case.values()) {
			classification2Counters.put(case2, 0);
		}

		int numWrongPeptides = 0;
		int lightPosInfinity = 0;
		int lightNegativeInfinity = 0;
		int heavyPositiveInfinity = 0;
		int heavyNegativeInfinity = 0;
		int peptidesSharedbyBothSpecies = 0;
		int peptidesSharedbyBothSpeciesAndInfinityRatio = 0;

		int numProteinPairs = 0;
		int numProt = 0;
		int numProtNodes = 0;
		int numPep = 0;
		int numPepNodes = 0;
		int numClustWithMoreThanOneIndividualProtein = 0;
		int numClustWithMoreThanOneProteinNode = 0;
		int numSignificantPeptideNodes = 0;
		int numSignificantClusters = 0;

		List<ProteinPairPValue> ranking = new ArrayList<ProteinPairPValue>();
		final ProteinClusterQuantParameters params = ProteinClusterQuantParameters.getInstance();
		try {
			final String outputPrefix = params.getOutputPrefix();
			final String outputSuffix = params.getOutputSuffix();
			final File outputFileFolder = params.getTemporalOutputFolder();

			if (params.isGenerateMiscellaneousFiles()) {
				outputCases = new FileWriter(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix
						+ "_cases_" + outputSuffix + ".txt");
				outputGroups = new FileWriter(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix
						+ "_groups_" + outputSuffix + ".txt");
				outputProteins = new FileWriter(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix
						+ "_prot_" + outputSuffix + ".txt");
				outputNames = new FileWriter(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix
						+ "_names_" + outputSuffix + ".txt");

				outputPairs = new FileWriter(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix
						+ "_proteinPairs_" + outputSuffix + ".txt");

				outputGroups.append("Prot" + "\t" + "Edge" + "\t" + "Pep" + "\n");
				outputNames.append("OriginalName" + "\t" + "NewName" + "\t" + "Description" + "\t" + "Species" + "\t"
						+ "ProtOrPep" + "\n");
				outputCases.append("Protein" + "\t" + "Case" + "\n");
			}
			outputSummaryFile = new File(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix
					+ "_summary_" + outputSuffix + ".txt");
			outputSummary = new FileWriter(outputSummaryFile);

			log.info("Classifying cases...");
			// to calculate the mean and std of ratios
			List<Double> peptideNodeConsensusRatios = new ArrayList<Double>();
			for (ProteinCluster cluster : clusterSet) {
				boolean containsSignificantPeptideNodes = false;

				Set<PCQProteinNode> proteinNodeSet = cluster.getProteinNodes();
				Set<PCQPeptideNode> peptideNodeSet = cluster.getPeptideNodes();
				if (cluster.getProteinSet().size() > 1) {
					numClustWithMoreThanOneIndividualProtein++;
				}
				if (proteinNodeSet.size() > 1) {
					numClustWithMoreThanOneProteinNode++;
				}
				numProtNodes += proteinNodeSet.size();
				numProt += cluster.getNumIndividualProteins();
				numPep += cluster.getPeptideSet().size();
				numPepNodes += cluster.getPeptideNodeKeys().size();
				if (params.isGenerateMiscellaneousFiles()) {
					for (PCQProteinNode proteinNode : proteinNodeSet) {
						Set<PCQPeptideNode> peptideNodeSet2 = proteinNode.getPeptideNodes();
						// prints out protein and peptides
						for (PCQPeptideNode peptideNode : peptideNodeSet2) {

							outputProteins.append(proteinNode.getAccession() + "\t");
							outputProteins.append(proteinNode.getDescription() + "\t");
							outputProteins.append(peptideNode.getSequence() + "\t");
							Double ratio = Utils.getRatioValue(peptideNode.getConsensusRatio(cond1, cond2), cond1,
									cond2);

							if (ratio == Double.NEGATIVE_INFINITY) {
								outputProteins.append("NEG_INF" + "\t");
							} else if (ratio == Double.POSITIVE_INFINITY) {
								outputProteins.append("POS_INF" + "\t");
							} else {

								outputProteins.append(ratio + "\t");
							}
							outputProteins.append("\t");
							outputProteins.append("1" + "\t");
							if (params.isGetTax()) {
								outputProteins.append("0" + "\t");
								outputProteins.append(proteinNode.getTaxonomy() + "\n");
							} else {
								outputProteins.append("0" + "\n");
							}
							outputProteins.flush();

						}
					}
				}

				for (PCQPeptideNode peptideNode : peptideNodeSet) {

					Double ratioValue = Utils.getRatioValue(peptideNode.getConsensusRatio(cond1, cond2), cond1, cond2);

					if (Double.isNaN(ratioValue)) {
						continue;
					}

					if (!Double.isInfinite(ratioValue)) {
						// to make the average and the stdev
						peptideNodeConsensusRatios.add(ratioValue);
					} else {
						// infinite nodes are considered as significantly
						// regulated
						numSignificantPeptideNodes++;
						containsSignificantPeptideNodes = true;
					}
					if (params.getSignificantFDRThreshold() != null) {

						final QuantRatio consensusRatio = peptideNode.getConsensusRatio(cond1, cond2);
						if (consensusRatio != null) {
							final Score score = consensusRatio.getAssociatedConfidenceScore();
							if (score != null) {
								try {
									if (score.getValue() != null) {
										double fdr = Double.valueOf(score.getValue());
										if (params.getSignificantFDRThreshold() >= fdr) {
											numSignificantPeptideNodes++;
											containsSignificantPeptideNodes = true;
										}
									}
								} catch (NumberFormatException e) {

								}
							}
						}
					}

					// see if Dm has neg or pos infinity associated
					if (peptideNode.getTaxonomies().contains(params.getLightSpecies())
							&& peptideNode.getTaxonomies().contains(params.getHeavySpecies())) {
						peptidesSharedbyBothSpecies++;
						if (Double.isInfinite(ratioValue)) {
							peptidesSharedbyBothSpeciesAndInfinityRatio++;
						}
					} else if (peptideNode.getTaxonomies().contains(params.getLightSpecies())) {

						if (Double.compare(Double.POSITIVE_INFINITY, ratioValue) == 0) {
							lightPosInfinity++;
						}
						if (Double.compare(Double.NEGATIVE_INFINITY, ratioValue) == 0) {
							lightNegativeInfinity++;
							numWrongPeptides++;
						}

					}

					// see if Dv has neg or pos infinity assocaited
					if (peptideNode.getTaxonomies().contains(params.getHeavySpecies())) {
						if (Double.compare(Double.POSITIVE_INFINITY, ratioValue) == 0) {
							heavyPositiveInfinity++;
							numWrongPeptides++;
						}
						if (Double.compare(Double.NEGATIVE_INFINITY, ratioValue) == 0) {
							heavyNegativeInfinity++;
						}
					}
				}
				if (containsSignificantPeptideNodes) {
					numSignificantClusters++;
				}
			}
			// iterate over cluster set (for each cluster)
			for (ProteinCluster cluster : clusterSet) {

				Set<ProteinPair> proteinPairs = cluster.getProteinPairs();

				for (ProteinPair proteinPair : proteinPairs) {
					proteinPair.setOutputCases(outputCases);
					proteinPair.setOutputGroups(outputGroups);
					proteinPair.setOutputNames(outputNames);

					numProteinPairs++;

					if (ProteinClusterQuantParameters.getInstance().isApplyClassificationsByProteinPair()) {
						proteinPair.compareAverages(cond1, cond2);
					}
					if (params.isGenerateMiscellaneousFiles()) {
						// print proteinPair summary
						for (String summaryLine : proteinPair.getSummaryLines(cond1, cond2)) {
							outputPairs.write(summaryLine + "\n");
						}
					}
					ProteinPairPValue first = proteinPair.getFirstCase();
					ProteinPairPValue second = proteinPair.getSecondCase();

					if (first != null && second != null) {
						ranking.add(first);
						ranking.add(second);
					}

					Collection<Classification1Case> pairCases1 = proteinPair.getClassification1Case().values();
					for (Classification1Case pairCase1 : pairCases1) {
						int count = classification1Counters.get(pairCase1) + 1;
						classification1Counters.put(pairCase1, count);
					}
					Collection<Classification2Case> pairCases2 = proteinPair.getClassification2Cases().values();
					for (Classification2Case pairCase2 : pairCases2) {
						int count = classification2Counters.get(pairCase2) + 1;
						classification2Counters.put(pairCase2, count);

					}
				}

				if (cluster.isInconsistenceWithQTest(cond1, cond2)) {
					// numInconsistenceClusters++;
				}

			}

			Collections.sort(ranking, getComparatorForProteinPairPValues());
			log.info("Printing statistics:");
			for (ProteinPairPValue proteinPairPValue : ranking) {
				log.info(proteinPairPValue.getpValue() + "\t");
				log.info(proteinPairPValue.getProteinpair().getNameProt1() + "\t");
				log.info(proteinPairPValue.getProteinpair().getNameProt2());
			}

			// log.info("Number Inconsistent with Cluster: " +
			// numInconsistenceClusters);
			StringBuilder stats = new StringBuilder();

			stats.append("\n-----------------\n");
			stats.append("Protein pair classification 1:\t\n");

			for (Classification1Case case1 : Classification1Case.values()) {
				stats.append("Number of Case " + case1.getCaseID() + " " + case1.getExplanation() + ":\t"
						+ classification1Counters.get(case1) + "\n");
			}

			stats.append("-----------------\n");
			stats.append("\n");

			stats.append("-----------------\n");
			stats.append("Protein pair classification 2:\t\n");

			for (Classification2Case case2 : Classification2Case.values()) {
				stats.append("Number of Case " + case2.getCaseID() + " " + case2.getExplanation() + ":\t"
						+ classification2Counters.get(case2) + "\n");
			}
			stats.append("-----------------\n");
			stats.append("\n");

			if (params.isIonsPerPeptideNodeThresholdOn()) {
				stats.append("Peptides discarded due to minimum number of ions ("
						+ params.getIonsPerPeptideNodeThreshold() + "):\t" + ProteinCluster.discardedPeptides + "\n\n");
			} else {
				stats.append("No peptides discarded due to minimum number of ions. Filter was disabled\n\n");
			}

			stats.append("Peptides shared by '" + params.getLightSpecies() + "' and '" + params.getHeavySpecies()
					+ "' species:\t" + peptidesSharedbyBothSpecies + "\n");
			stats.append("Peptides shared by '" + params.getLightSpecies() + "' and '" + params.getHeavySpecies()
					+ "' species having +/-INF ratio  (wrong):\t" + peptidesSharedbyBothSpeciesAndInfinityRatio + "\n");

			stats.append(
					"'" + params.getLightSpecies() + "' specific peptides with POS_INF:\t " + lightPosInfinity + "\n");
			stats.append("'" + params.getLightSpecies() + "' specific peptides with NEG_INF (wrong):\t "
					+ lightNegativeInfinity + "\n");
			stats.append("'" + params.getHeavySpecies() + "' specific peptides with POS_INF (wrong):\t "
					+ heavyPositiveInfinity + "\n");
			stats.append("'" + params.getHeavySpecies() + "' specific peptides with NEG_INF:\t " + heavyNegativeInfinity
					+ "\n");
			stats.append("Total number of peptides with WRONG ratio assignment:\t" + numWrongPeptides + "\n");
			stats.append("\n");

			stats.append("Number of Clusters:\t " + clusterSet.size() + "\n");
			stats.append("Number of Clusters with only one Protein:\t "
					+ (clusterSet.size() - numClustWithMoreThanOneIndividualProtein) + "\n");
			stats.append("Number of Clusters with only one Protein Node:\t "
					+ (clusterSet.size() - numClustWithMoreThanOneProteinNode) + "\n");
			stats.append("Number of Clusters with more than one Protein:\t " + numClustWithMoreThanOneIndividualProtein
					+ "\n");
			stats.append("Number of Clusters with more than one Protein Node:\t " + numClustWithMoreThanOneProteinNode
					+ "\n");
			stats.append("Number of ProteinPairs:\t " + numProteinPairs + "\n");
			stats.append("Number of Proteins:\t " + numProt + "\n");
			stats.append("Number of Proteins Nodes:\t " + numProtNodes + "\n");
			stats.append("Number of Peptides:\t " + numPep + "\n");
			stats.append("Number of Peptides Nodes:\t " + numPepNodes + "\n");
			// if (params.getSignificantFDRThreshold() != null) {
			stats.append("Number of significantly changing peptide nodes (FDR<" + params.getSignificantFDRThreshold()
					+ "):\t" + numSignificantPeptideNodes + "\n");
			stats.append("Number of significantly changing protein clusters (FDR<" + params.getSignificantFDRThreshold()
					+ "):\t" + numSignificantClusters + "\n");
			// }

			stats.append(
					"Number of Peptides NOT found in DB: " + AbstractQuantParser.peptidesMissingInDB.size() + "\n");
			final double mean = Maths.mean(peptideNodeConsensusRatios.toArray(new Double[0]));
			final double stdev = Maths.stddev(peptideNodeConsensusRatios.toArray(new Double[0]));
			stats.append("Average of peptide node ratios: " + mean + "\n");
			stats.append("Standard deviation of peptide node ratios: " + stdev + "\n");
			// print to console
			log.info(stats.toString());
			// print to file
			outputSummary.write(stats.toString() + "\n");
		} catch (

		IOException e) {
			log.error(e.getMessage());
		} finally {
			try {
				if (outputCases != null)
					outputCases.close();
				if (outputProteins != null)
					outputProteins.close();
				if (outputGroups != null)
					outputGroups.close();
				if (outputNames != null)
					outputNames.close();
				if (outputSummary != null)
					outputSummary.close();
				if (outputPairs != null)
					outputPairs.close();
				// append parameters file into the summary file
				appendFiles(outputSummaryFile, setupPropertiesFile);
			} catch (IOException e) {
				log.error(e.getMessage());
			}

		}
		log.info("Statistics printed");
	}

	private void appendFiles(File file1, File file2) {

		try {
			// File to write
			File file3;
			file3 = File.createTempFile("TEMP", "tmp");
			file3.deleteOnExit();
			FileWriter filewriter = new FileWriter(file3, true);

			// Read the file as string
			String file1Str = FileUtils.readFileToString(file1);
			String file2Str = FileUtils.readFileToString(file2);

			// Write the file
			filewriter.write(file1Str);
			filewriter.write("\n\n\nINPUT PARAMETERS:\n\n");
			filewriter.write(file2Str); // true for append
			filewriter.close();
			FileUtils.copyFile(file3, file1);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	/**
	 * Given a set of peptides, returns the set of psms detected in a given
	 * replicate
	 *
	 * @param replicateName
	 * @param peptideSet
	 * @return
	 */
	private Set<QuantifiedPSMInterface> getPSMsFromReplicate(String replicateName,
			Set<QuantifiedPeptideInterface> peptideSet) {
		Set<QuantifiedPSMInterface> ret = new HashSet<QuantifiedPSMInterface>();
		for (QuantifiedPeptideInterface quantifiedPeptide : peptideSet) {
			final Set<QuantifiedPSMInterface> quantifiedPSMs = quantifiedPeptide.getQuantifiedPSMs();
			for (QuantifiedPSMInterface quantifiedPSM : quantifiedPSMs) {
				final String fileName = quantifiedPSM.getRawFileName();
				if (fileName.contains(replicateName)) {
					ret.add(quantifiedPSM);
				} else {
					if (replicateName.equals("E1")) {
						ret.add(quantifiedPSM);
					}
				}
			}
		}
		return ret;
	}

	private Set<String> getReplicateNamesFromPeptides(Set<QuantifiedPeptideInterface> peptideSet) {
		Set<String> ret = new HashSet<String>();
		for (QuantifiedPeptideInterface quantifiedPeptide : peptideSet) {
			final Set<String> fileNames = quantifiedPeptide.getFileNames();
			for (String string : fileNames) {
				boolean someAdded = false;
				final String[] replicateIdentifiers = ProteinClusterQuantParameters.getInstance()
						.getInputFileNamesArray();
				if (replicateIdentifiers != null) {
					for (String replicateIdentifier : replicateIdentifiers) {
						if (string.contains(replicateIdentifier)) {
							ret.add(replicateIdentifier);
							someAdded = true;
						}
					}
				}
				if (!someAdded) {
					ret.add(string);
				}
				// FOR SERIES:
				// if (string.contains("_A1_")) {
				// ret.add("A1");
				// } else if (string.contains("_A2_")) {
				// ret.add("A2");
				// } else if (string.contains("_A3_")) {
				// ret.add("A3");
				// } else {
				// ret.add("A1");
				// }

				// for DmDv
				// ret.add(string);
			}
		}
		return ret;
	}

	/**
	 * Gets a comparator to use for sorting {@link ProteinPairPValue} by the
	 * pvalue
	 *
	 * @return
	 */
	private Comparator<ProteinPairPValue> getComparatorForProteinPairPValues() {
		Comparator<ProteinPairPValue> comparator = new Comparator<ProteinPairPValue>() {

			@Override
			public int compare(ProteinPairPValue o1, ProteinPairPValue o2) {
				return Double.compare(o1.getpValue(), o2.getpValue());
			}
		};
		return comparator;
	}

}
