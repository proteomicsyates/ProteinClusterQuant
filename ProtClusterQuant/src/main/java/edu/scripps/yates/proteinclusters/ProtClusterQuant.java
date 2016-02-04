package edu.scripps.yates.proteinclusters;

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

import javax.xml.bind.JAXBException;

import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.annotations.uniprot.xml.Entry;
import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.CensusChroParser;
import edu.scripps.yates.census.read.model.IsobaricQuantifiedPeptide;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.proteinclusters.util.ColorManager;
import edu.scripps.yates.proteinclusters.util.PropertiesReader;
import edu.scripps.yates.proteinclusters.util.ProteinClusterQuantParameters;
import edu.scripps.yates.proteinclusters.util.Utils;
import edu.scripps.yates.proteinclusters.xgmml.XgmmlExporter;
import edu.scripps.yates.utilities.alignment.nwalign.NWResult;

public class ProtClusterQuant {
	private final static Logger log = Logger.getLogger(ProtClusterQuant.class);
	private static final String SETUP_PROPERTIES = "setup.properties";
	private final QuantCondition cond1 = new QuantCondition("cond1");
	private final QuantCondition cond2 = new QuantCondition("cond2");
	// private static final String SETUP_PROPERTIES = "setup.properties";
	// CASIMIR
	private Map<String, Set<NWResult>> peptideAlignments;
	private final File setupPropertiesFile;

	public ProtClusterQuant(File setupPropertiesFile) {
		this.setupPropertiesFile = setupPropertiesFile;
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
		Map<String, ProteinCluster> clusterMap = new HashMap<String, ProteinCluster>();

		final ProteinClusterQuantParameters params = ProteinClusterQuantParameters.getInstance();
		try {

			List<Map<QuantCondition, QuantificationLabel>> labelsByConditionsList = new ArrayList<Map<QuantCondition, QuantificationLabel>>();

			// per L/H xml file
			Map<QuantCondition, QuantificationLabel> labelsByConditions = new HashMap<QuantCondition, QuantificationLabel>();
			labelsByConditions.put(cond1, QuantificationLabel.LIGHT);
			labelsByConditions.put(cond2, QuantificationLabel.HEAVY);
			log.info("Input parameters:\n" + params);
			int max = params.getInputFileNames().length;
			if (params.isLabelSwap()) {
				max = max / 2;
			}
			for (int i = 0; i < max; i++) {
				labelsByConditionsList.add(labelsByConditions);
			}
			// labelsByConditionsList.add(labelsByConditions);
			// labelsByConditionsList.add(labelsByConditions);
			// labelsByConditionsList.add(labelsByConditions);

			// per H/L xml file
			if (params.isLabelSwap()) {
				Map<QuantCondition, QuantificationLabel> labelsByConditions2 = new HashMap<QuantCondition, QuantificationLabel>();
				labelsByConditions2.put(cond1, QuantificationLabel.HEAVY);
				labelsByConditions2.put(cond2, QuantificationLabel.LIGHT);
				for (int i = 0; i < max; i++) {
					labelsByConditionsList.add(labelsByConditions2);
				}
			}
			log.info("Reading input files...");
			CensusChroParser parser = Utils.getCensusChroParser(params.getFastaFile(), params.getInputFileFolder(),
					params.getInputFileNames(), labelsByConditionsList, params.getEnzymeArray(),
					params.getMissedCleavages(), params.getUniprotReleasesFolder(), params.getDecoyRegexp());

			// if (params.isIonsPerPeptideThresholdOn()) {
			// log.info("Filtering peptides by number of ions
			// data...");
			// log.info("Requiring " +
			// params.getIonsPerPeptideThreshold() + " ions per peptide...");
			// if (params.isIonsPerPeptideThresholdOn()) {
			// final Collection<String> peptideSequences =
			// parser.getPeptideMap().keySet();
			// Set<String> sequencesToRemove = new HashSet<String>();
			// for (String peptideSequence : peptideSequences) {
			// QuantifiedPeptide quantifiedPeptide =
			// parser.getPeptideMap().get(peptideSequence);
			// if (Utils.getIonCount(quantifiedPeptide) <
			// params.getIonsPerPeptideThreshold()) {
			// numPeptidesDiscarded++;
			// // remove this peptide from its proteins
			// final Set<QuantifiedProteinInterface> quantifiedProteins =
			// quantifiedPeptide.getQuantifiedProteins();
			// for (QuantifiedProteinInterface quantifiedProtein :
			// quantifiedProteins) {
			// quantifiedProtein.getQuantifiedPeptides().remove(quantifiedPeptide);
			// }
			// // get all psms of that peptide
			// final Set<QuantifiedPSM> quantifiedPSMs =
			// quantifiedPeptide.getQuantifiedPSMs();
			// for (QuantifiedPSM quantifiedPSM : quantifiedPSMs) {
			// // remove this psms from its proteins
			// final Iterator<QuantifiedProteinInterface> quantifiedProteins2 =
			// quantifiedPSM
			// .getQuantifiedProteins().iterator();
			// while (quantifiedProteins2.hasNext()) {
			// final QuantifiedProteinInterface quantifiedProtein =
			// quantifiedProteins2.next();
			// quantifiedProtein.getQuantifiedPSMs().remove(quantifiedPSM);
			// }
			// // remove this psm from teh parser
			// parser.getPSMMap().remove(quantifiedPSM.getPSMIdentifier());
			//
			// }
			// sequencesToRemove.add(peptideSequence);
			// }
			// }
			// for (String sequenceToRemove : sequencesToRemove) {
			// parser.getPeptideMap().remove(sequenceToRemove);
			// }
			// }
			//
			// log.info(numPeptidesDiscarded + " peptides discarded");
			// }

			// UniProt_D_simulans_and_melanogaster_11-01-2014.fasta
			// UniProt_drosophila_melanogaster_and_virilis_11-01-2014_reversed.fasta
			// UniProt_Human_DF508CFTR_06-03-2010_reversed.fasta
			// UniProt_Human_02-09-2013_reversed.fasta

			Map<String, QuantifiedPeptideInterface> pepMap = parser.getPeptideMap();
			// gets map from peptide map
			// Map<String, QuantifiedPeptide> pepMap =
			// parser.getPeptideMap();

			// List to hold all peptides
			List<QuantifiedPeptideInterface> peptideList = new ArrayList<QuantifiedPeptideInterface>();
			peptideList.addAll(parser.getPeptideMap().values());

			peptideAlignments = new HashMap<String, Set<NWResult>>();
			final File outputFileFolder = params.getTemporalOutputFolder();
			final String outputPrefix = params.getOutputPrefix();
			final String outputSuffix = params.getOutputSuffix();
			if (params.isMakeAlignments()) {
				FileWriter output = null;
				try {
					output = new FileWriter(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix + "_pep_"
							+ outputSuffix + ".txt");
					log.info("Aligning peptides...");
					long t = System.currentTimeMillis();
					peptideAlignments = Utils.alignPeptides(peptideList, cond1, cond2, output);
					log.info("Alignments done in " + Utils.format((t - System.currentTimeMillis()) * 1.0 / 1000.0)
							+ "sg.");
				} finally {
					if (output != null)
						output.close();
				}
			}
			int numClusters = 0;
			int i = 0;
			while (i < 2) {
				clusterSet = cluster(pepMap);
				log.info("Identifying protein pairs in " + clusterSet.size() + " clusters...");
				int numProteinPairs = 0;
				for (ProteinCluster cluster : clusterSet) {
					if (params.isCollapseIndistinguishableProteins()) {
						cluster.createPairsCollapsingIndistinguisibleProteins(peptideAlignments);
					} else {
						cluster.createPairs(peptideAlignments);
					}
					numProteinPairs += cluster.getProteinPairs().size();
				}
				log.info(numProteinPairs + " protein pairs identified.");

				if (params.isIonsPerPeptideThresholdOn()) {
					log.info("Filtering peptide nodes containing less than " + params.getIonsPerPeptideThreshold()
							+ " ions");
					final Iterator<ProteinCluster> iterator = clusterSet.iterator();
					while (iterator.hasNext()) {
						ProteinCluster cluster = iterator.next();
						// filter out the nodes that doesn't pass the threshold
						// of minimum number of ions
						if (params.isIonsPerPeptideThresholdOn()) {
							cluster.applyIonsThreshold(params.getIonsPerPeptideThreshold());
						}
						// check if the cluster still contains proteins and
						// peptides.
						if (cluster.getPeptideSet().isEmpty()) {
							iterator.remove();
						}
					}
				}
				pepMap = Utils.getPeptideMapFromClusters(clusterSet);
				if (numClusters != clusterSet.size()) {
					i = 0;
				} else {
					i++;
				}
				numClusters = clusterSet.size();
			}
			log.info("Final number of clusters after iterations: \t" + clusterSet.size());

			// int numInconsistenceClusters = 0;

			// print statistics and output files
			classiffyAndPrintStatistics(clusterSet);

			log.info("Getting UniprotKB annotations for " + parser.getProteinMap().size() + " proteins");
			// get Uniprot annotations for retrieving geneNames and taxonomies
			UniprotProteinLocalRetriever uplr = null;
			if (params.getUniprotReleasesFolder() != null) {
				uplr = new UniprotProteinLocalRetriever(params.getUniprotReleasesFolder(), true);
			}
			// exporting to xgmml
			final Map<String, Entry> annotatedProteins = uplr.getAnnotatedProteins(null,
					parser.getProteinMap().keySet());
			log.info(annotatedProteins.size() + " annotations retrieved out of " + parser.getProteinMap().size()
					+ " proteins");
			// print PSEA QUANT files
			if (params.getReplicateIdentifiers() != null) {
				if ("".equals(params.getReplicateIdentifiers()[0])) {
					log.info(
							"Replicate names are empty. PSEA-Quant files will contain a column for each original RAW file");
				}
				printPSEAQuantFiles(clusterSet, annotatedProteins);
			} else {
				log.info("PSEA-Quant intput files will NOT be created");
			}

			exportToXGMML(outputFileFolder, outputPrefix, outputSuffix, cond1, cond2, params.getColorManager(),
					clusterSet, annotatedProteins);

			// rename TEMP output folder to output folder
			log.info("Moving results to final output folder...");
			FileUtils.copyDirectory(params.getTemporalOutputFolder(), params.getOutputFileFolder());
			log.info("Deleting TEMP folder...");
			FileUtils.deleteDirectory(params.getTemporalOutputFolder());
			log.info("DONE.");
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			log.error(e.getMessage());
		}

	}

	private Set<ProteinCluster> cluster(Map<String, QuantifiedPeptideInterface> peptideMap) {
		Map<String, ProteinCluster> clusterMap = new HashMap<String, ProteinCluster>();
		Set<ProteinCluster> clusterSet = new HashSet<ProteinCluster>();
		log.info("Starting clustering " + peptideMap.size() + " peptides...");
		long t0 = System.currentTimeMillis();
		for (QuantifiedPeptideInterface peptide1 : peptideMap.values()) {
			String sequence1 = peptide1.getSequence();
			ProteinCluster cluster = null;

			// if peptide has a cluster associated
			if (clusterMap.containsKey(sequence1)) {
				// grab it
				cluster = clusterMap.get(sequence1);
			} else {
				cluster = new ProteinCluster();
				// adds cluster to cluster set
				clusterSet.add(cluster);
				// Map<String, protCluster> <- cluster (key = sequence)
				clusterMap.put(sequence1, cluster);
			}

			// put peptides in cluster
			cluster.addPeptides(peptide1);
			if (peptideAlignments.containsKey(sequence1)) {
				Set<NWResult> gAS = peptideAlignments.get(sequence1);
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
							if (clusterMap.containsKey(pep2.getSequence())) {
								ProteinCluster cluster2 = clusterMap.get(pep2.getSequence());
								if (!cluster.equals(cluster2)) {
									// merges the clusters with similar peptides
									cluster = Utils.mergeClusters(cluster, cluster2);
									clusterSet.remove(cluster2);
									for (QuantifiedPeptideInterface quantifiedPeptide : cluster.getPeptideSet()) {
										clusterMap.put(quantifiedPeptide.getSequence(), cluster);
									}
								}

								// Add (Quant pep2, hisProtein) to Cluster)
								cluster.addPeptides(pep2);
								clusterMap.put(pep2.getSequence(), cluster);
							}
							// add alignment to the cluster
							cluster.addAlignment(AlignResult);
						}
					}
				}
			}
			// get proteins of the peptide
			Set<QuantifiedProteinInterface> proteinSet = peptide1.getQuantifiedProteins();
			for (QuantifiedProteinInterface protein : proteinSet) {
				// put protein in cluster
				cluster.addProteins(protein);

				// peptide 2 <- protein
				for (QuantifiedPeptideInterface peptide2 : protein.getQuantifiedPeptides()) {
					// checking to see if peptide 2 is already in a
					// cluster
					if (clusterMap.containsKey(peptide2.getSequence())) {
						ProteinCluster cluster2 = clusterMap.get(peptide2.getSequence());
						if (!cluster.equals(cluster2)) {
							// merges the clusters with similar peptides
							cluster = Utils.mergeClusters(cluster, cluster2);
							clusterSet.remove(cluster2);
							for (QuantifiedPeptideInterface quantifiedPeptide : cluster.getPeptideSet()) {
								clusterMap.put(quantifiedPeptide.getSequence(), cluster);
							}
						}
					}

					// add new peptides to cluster
					cluster.addPeptides(peptide2);

					// Map <- peptide, cluster
					clusterMap.put(peptide2.getSequence(), cluster);
				}
			}
		}
		double time = System.currentTimeMillis() - t0;
		String units = "sg";
		if (time < 1000) {
			units = "ms";
		} else {
			time = time * 1.0 / 1000;
		}
		log.info("Clustering done in " + Utils.format(time) + " " + units);
		return clusterSet;
	}

	private void printPSEAQuantFiles(Set<ProteinCluster> clusterSet, Map<String, Entry> annotatedProteins) {
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
				replicateNames.addAll(getReplicateNames(peptideSet));

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
				// according to the manuscript, it is the average of the ratios
				// of the peptide nodes.

				for (String replicateName : replicateNameList) {
					final Set<ProteinPair> proteinPairs = cluster.getProteinPairs();
					Set<Double> ratioValues = new HashSet<Double>();
					if (proteinPairs.isEmpty()) {
						final Ratio pepRatio = Utils.getPepRatio(
								cluster.getPeptideSet(
										ProteinClusterQuantParameters.getInstance().isExcludeUniquePeptides()),
								cond1, cond2, replicateName);
						if (pepRatio != null) {
							final Double countRatio = pepRatio.getCountRatio(cond1, cond2);
							if (countRatio != null && !Double.isNaN(countRatio) && !Double.isInfinite(countRatio)) {
								ratioValues.add(countRatio);
							}
						}
					} else {
						for (ProteinPair proteinPair : proteinPairs) {
							final Map<String, Set<QuantifiedPeptideInterface>> sharedPeptidesMap_S12 = Utils
									.getSharedPeptidesMap(proteinPair.getProtein1(), proteinPair.getProtein2(), false);
							for (Set<QuantifiedPeptideInterface> peptidesS12 : sharedPeptidesMap_S12.values()) {
								List<Set<QuantifiedPeptideInterface>> pep12DoubleList = new ArrayList<Set<QuantifiedPeptideInterface>>();
								if (ProteinClusterQuantParameters.getInstance().isCollapseIndistinguishablePeptides()) {
									pep12DoubleList.add(peptidesS12);
								} else {
									for (QuantifiedPeptideInterface peptide : peptidesS12) {
										Set<QuantifiedPeptideInterface> set = new HashSet<QuantifiedPeptideInterface>();
										set.add(peptide);
										pep12DoubleList.add(set);
									}
								}
								for (Set<QuantifiedPeptideInterface> sharedPeptides_S12 : pep12DoubleList) {
									final Ratio pepRatio = Utils.getPepRatio(sharedPeptides_S12, cond1, cond2,
											replicateName);
									if (pepRatio != null) {
										final Double countRatio = pepRatio.getCountRatio(cond1, cond2);
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

			log.info("Clusters skipped due to not having appropiate gene name: " + numClustersSkippedNotHavingGeneName);
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

	}

	private void classiffyAndPrintStatistics(Set<ProteinCluster> clusterSet) {

		FileWriter Output = null;
		FileWriter Output5 = null;
		FileWriter Output3 = null;
		FileWriter OutputNames = null;
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
		int numClustMoreOne = 0;
		List<ProteinPairPValue> ranking = new ArrayList<ProteinPairPValue>();
		final ProteinClusterQuantParameters params = ProteinClusterQuantParameters.getInstance();
		try {
			final String outputPrefix = params.getOutputPrefix();
			final String outputSuffix = params.getOutputSuffix();
			final File outputFileFolder = params.getTemporalOutputFolder();
			Output = new FileWriter(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix + "_cases_"
					+ outputSuffix + ".txt");
			Output5 = new FileWriter(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix + "_groups_"
					+ outputSuffix + ".txt");
			Output3 = new FileWriter(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix + "_prot_"
					+ outputSuffix + ".txt");
			OutputNames = new FileWriter(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix + "_names_"
					+ outputSuffix + ".txt");
			outputSummaryFile = new File(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix
					+ "_summary_" + outputSuffix + ".txt");
			outputPairs = new FileWriter(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix
					+ "_proteinPairs_" + outputSuffix + ".txt");
			outputSummary = new FileWriter(outputSummaryFile);

			Output5.append("Prot" + "\t" + "Edge" + "\t" + "Pep" + "\n");
			OutputNames.append("OriginalName" + "\t" + "NewName" + "\t" + "Description" + "\t" + "Species" + "\t"
					+ "ProtOrPep" + "\n");
			Output.append("Protein" + "\t" + "Case" + "\n");

			log.info("Classifying cases...");

			for (ProteinCluster cluster : clusterSet) {
				Set<QuantifiedProteinInterface> proteinSet = cluster.getProteinSet();
				Set<QuantifiedPeptideInterface> peptideSet = cluster.getPeptideSet();
				if (proteinSet.size() > 1) {
					numClustMoreOne++;
				}
				numProtNodes += cluster.getProteinSet().size();
				numProt += cluster.getNumIndividualProteins();
				numPep += cluster.getPeptideSet().size();
				numPepNodes += cluster.getPeptideNodeKeys().size();
				for (QuantifiedProteinInterface quantifiedProtein : proteinSet) {

					Set<QuantifiedPeptideInterface> peptideSet2 = quantifiedProtein.getQuantifiedPeptides();
					// prints out protein and peptides
					for (QuantifiedPeptideInterface peptide : peptideSet2) {
						Output3.append(quantifiedProtein.getAccession() + "\t");
						Output3.append(quantifiedProtein.getDescription() + "\t");
						Output3.append(peptide.getSequence() + "\t");
						Double ratio = null;
						if (peptide instanceof IsobaricQuantifiedPeptide) {
							ratio = ((IsobaricQuantifiedPeptide) peptide).getCountRatio(cond1, cond2);
						}
						if (ratio == Double.NEGATIVE_INFINITY) {
							Output3.append("NEG_INF" + "\t");
						} else if (ratio == Double.POSITIVE_INFINITY) {
							Output3.append("POS_INF" + "\t");
						} else {
							Output3.append(ratio + "\t");
						}
						Output3.append("\t");
						Output3.append("1" + "\t");
						if (params.isGetTax()) {
							Output3.append("0" + "\t");
							Output3.append(quantifiedProtein.getTaxonomy() + "\n");
						} else {
							Output3.append("0" + "\n");
						}
						Output3.flush();

					}

				}

				for (QuantifiedPeptideInterface quantifiedPeptide : peptideSet) {

					Double ratio = null;
					if (quantifiedPeptide instanceof IsobaricQuantifiedPeptide) {
						ratio = ((IsobaricQuantifiedPeptide) quantifiedPeptide).getCountRatio(cond1, cond2);
					}
					if (Double.isNaN(ratio)) {
						continue;
					}

					// see if Dm has neg or pos infinity associated
					if (quantifiedPeptide.getTaxonomies().contains(params.getLightSpecies())
							&& quantifiedPeptide.getTaxonomies().contains(params.getHeavySpecies())) {
						peptidesSharedbyBothSpecies++;
						if (Double.isInfinite(ratio)) {
							peptidesSharedbyBothSpeciesAndInfinityRatio++;
						}
					} else if (quantifiedPeptide.getTaxonomies().contains(params.getLightSpecies())) {

						if (Double.compare(Double.POSITIVE_INFINITY, ratio) == 0) {
							lightPosInfinity++;
						}
						if (Double.compare(Double.NEGATIVE_INFINITY, ratio) == 0) {
							lightNegativeInfinity++;
							numWrongPeptides++;
						}

					}

					// see if Dv has neg or pos infinity assocaited
					if (quantifiedPeptide.getTaxonomies().contains(params.getHeavySpecies())) {
						if (Double.compare(Double.POSITIVE_INFINITY, ratio) == 0) {
							heavyPositiveInfinity++;
							numWrongPeptides++;
						}
						if (Double.compare(Double.NEGATIVE_INFINITY, ratio) == 0) {
							heavyNegativeInfinity++;
						}
					}
				}

			}
			// iterate over cluster set (for each cluster)
			for (ProteinCluster cluster : clusterSet) {

				Set<ProteinPair> proteinPairs = cluster.getProteinPairs();

				for (ProteinPair proteinPair : proteinPairs) {
					proteinPair.setOutput(Output);
					proteinPair.setOutput5(Output5);
					proteinPair.setOutputNames(OutputNames);

					numProteinPairs++;
					proteinPair.compareAverages(cond1, cond2);

					// print proteinPair summary
					for (String summaryLine : proteinPair.getSummaryLines(cond1, cond2)) {
						outputPairs.write(summaryLine + "\n");
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

			if (params.isIonsPerPeptideThresholdOn()) {
				stats.append("Peptides discarded due to minimum number of ions (" + params.getIonsPerPeptideThreshold()
						+ "):\t" + ProteinCluster.discardedPeptides + "\n\n");
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
			stats.append("Number of Clusters with only one Protein:\t " + (clusterSet.size() - numClustMoreOne) + "\n");
			stats.append("Number of Clusters with more than one Protein:\t " + numClustMoreOne + "\n");
			stats.append("Number of ProteinPairs:\t " + numProteinPairs + "\n");
			stats.append("Number of Proteins:\t " + numProt + "\n");
			stats.append("Number of Proteins Nodes:\t " + numProtNodes + "\n");
			stats.append("Number of Peptides:\t " + numPep + "\n");
			stats.append("Number of Peptides Nodes:\t " + numPepNodes + "\n");

			// print to console
			log.info(stats.toString());
			// print to file
			outputSummary.write(stats.toString() + "\n");
		} catch (IOException e) {
			log.error(e.getMessage());
		} finally {
			try {
				if (Output != null)
					Output.close();
				if (Output3 != null)
					Output3.close();
				if (Output5 != null)
					Output5.close();
				if (OutputNames != null)
					OutputNames.close();
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

	private void exportToXGMML(File outputFileFolder, String outputPrefix, String outputSuffix, QuantCondition cond1,
			QuantCondition cond2, ColorManager colorManager, Collection<ProteinCluster> clusterCollection,
			Map<String, Entry> annotatedProteins) {
		final XgmmlExporter xgmmlExporter = new XgmmlExporter();
		xgmmlExporter.setAnnotatedProteins(annotatedProteins);
		log.info("Creating XGMML files for Cytoscape...");
		try {
			log.info("Creating XGMML for the entire network...");
			// export the total network
			File xgmmlOutPutFile = new File(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix
					+ "_cytoscape_" + outputSuffix + ".xgmml");
			xgmmlExporter.exportToGmmlFromProteinClusters(xgmmlOutPutFile, outputPrefix + "_" + outputSuffix,
					clusterCollection, cond1, cond2, colorManager);
			// classification 1
			final Map<Classification1Case, Set<ProteinCluster>> proteinPairsByClassification1 = getProteinClustersByClassification1(
					clusterCollection);
			final Classification1Case[] cases1 = Classification1Case.values();
			for (Classification1Case case1 : cases1) {
				if (proteinPairsByClassification1.containsKey(case1)) {
					log.info("Creating XGMML for case " + case1.getCaseID() + "(" + case1.getExplanation() + ")...");

					File xgmmlOutPutFile2 = new File(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix
							+ "_cytoscape_" + case1.getCaseID() + "-" + case1.name() + "_" + outputSuffix + ".xgmml");
					xgmmlExporter.exportToGmmlFromProteinClusters(xgmmlOutPutFile2,
							outputPrefix + "_" + case1.getCaseID() + "-" + case1.name() + "_" + outputSuffix,
							proteinPairsByClassification1.get(case1), cond1, cond2, colorManager);
				}
			}
			// classification 2
			final Map<Classification2Case, Set<ProteinCluster>> proteinPairsByClassification2 = getProteinClustersByClassification2(
					clusterCollection);
			final Classification2Case[] cases2 = Classification2Case.values();
			for (Classification2Case case2 : cases2) {
				if (proteinPairsByClassification2.containsKey(case2)) {
					log.info("Creating XGMML for case " + case2.getCaseID() + "(" + case2.getExplanation() + ")...");
					File xgmmlOutPutFile2 = new File(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix
							+ "_cytoscape_" + case2.getCaseID() + "-" + case2.name() + "_" + outputSuffix + ".xgmml");
					xgmmlExporter.exportToGmmlFromProteinClusters(xgmmlOutPutFile2,
							outputPrefix + "_" + case2.getCaseID() + "-" + case2.name() + "_" + outputSuffix,
							proteinPairsByClassification2.get(case2), cond1, cond2, colorManager);
				}
			}

		} catch (JAXBException e) {
			e.printStackTrace();
			log.error(e.getMessage());
		}

	}

	private Map<Classification1Case, Set<ProteinCluster>> getProteinClustersByClassification1(
			Collection<ProteinCluster> clusters) {
		Map<Classification1Case, Set<ProteinCluster>> ret = new HashMap<Classification1Case, Set<ProteinCluster>>();
		for (ProteinCluster proteinCluster : clusters) {
			final Set<ProteinPair> proteinPairs = proteinCluster.getProteinPairs();
			for (ProteinPair proteinPair : proteinPairs) {
				final Collection<Classification1Case> classification1PairCases = proteinPair.getClassification1Case()
						.values();
				for (Classification1Case classification1PairCase : classification1PairCases) {
					if (ret.containsKey(classification1PairCase)) {
						ret.get(classification1PairCase).add(proteinCluster);
					} else {
						Set<ProteinCluster> set = new HashSet<ProteinCluster>();
						set.add(proteinCluster);
						ret.put(classification1PairCase, set);
					}
				}
			}
		}
		return ret;
	}

	private Map<Classification2Case, Set<ProteinCluster>> getProteinClustersByClassification2(
			Collection<ProteinCluster> clusters) {
		Map<Classification2Case, Set<ProteinCluster>> ret = new HashMap<Classification2Case, Set<ProteinCluster>>();
		for (ProteinCluster proteinCluster : clusters) {
			final Set<ProteinPair> proteinPairs = proteinCluster.getProteinPairs();
			for (ProteinPair proteinPair : proteinPairs) {
				final Collection<Classification2Case> classification2PairCases = proteinPair.getClassification2Cases()
						.values();
				for (Classification2Case classification2PairCase : classification2PairCases) {
					if (ret.containsKey(classification2PairCase)) {
						ret.get(classification2PairCase).add(proteinCluster);
					} else {
						Set<ProteinCluster> set = new HashSet<ProteinCluster>();
						set.add(proteinCluster);
						ret.put(classification2PairCase, set);
					}
				}
			}
		}
		return ret;
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
				final String fileName = quantifiedPSM.getFileName();
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

	private Set<String> getReplicateNames(Set<QuantifiedPeptideInterface> peptideSet) {
		Set<String> ret = new HashSet<String>();
		for (QuantifiedPeptideInterface quantifiedPeptide : peptideSet) {
			final Set<String> fileNames = quantifiedPeptide.getRawFileNames();
			for (String string : fileNames) {
				boolean someAdded = false;
				final String[] replicateIdentifiers = ProteinClusterQuantParameters.getInstance()
						.getReplicateIdentifiers();
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

			public int compare(ProteinPairPValue o1, ProteinPairPValue o2) {
				return Double.compare(o1.getpValue(), o2.getpValue());
			}
		};
		return comparator;
	}

}
