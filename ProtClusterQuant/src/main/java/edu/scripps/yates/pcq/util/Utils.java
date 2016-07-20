package edu.scripps.yates.proteinclusters.util;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.annotations.uniprot.xml.Entry;
import edu.scripps.yates.annotations.uniprot.xml.GeneNameType;
import edu.scripps.yates.annotations.uniprot.xml.GeneType;
import edu.scripps.yates.annotations.uniprot.xml.OrganismNameType;
import edu.scripps.yates.annotations.uniprot.xml.OrganismType;
import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.CensusChroParser;
import edu.scripps.yates.census.read.model.Ion;
import edu.scripps.yates.census.read.model.IonSerie.IonSerieType;
import edu.scripps.yates.census.read.model.IsobaricQuantifiedPSM;
import edu.scripps.yates.census.read.model.IsobaricQuantifiedPeptide;
import edu.scripps.yates.census.read.model.QuantifiedProteinFromCensusOut;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.dbindex.DBIndexInterface;
import edu.scripps.yates.dbindex.io.DBIndexSearchParams;
import edu.scripps.yates.dbindex.io.DBIndexSearchParamsImpl;
import edu.scripps.yates.proteinclusters.ProteinCluster;
import edu.scripps.yates.proteinclusters.Ratio;
import edu.scripps.yates.utilities.alignment.nwalign.NWAlign;
import edu.scripps.yates.utilities.alignment.nwalign.NWResult;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;

public class Utils {

	public static final String PROTEIN_DESCRIPTION_SEPARATOR = "####";
	public static DecimalFormat df = new DecimalFormat("#.#");
	private final static Logger log = Logger.getLogger(Utils.class);
	public static final String PROTEIN_ACC_SEPARATOR = " ";
	public static final double factor = 1.2;

	public static Map<String, Set<NWResult>> alignPeptides(List<QuantifiedPeptideInterface> peptideList,
			QuantCondition cond1, QuantCondition cond2, FileWriter Output2) throws IOException {
		Set<NWResult> gAS = new HashSet<NWResult>();
		Map<String, Set<NWResult>> gAM = new HashMap<String, Set<NWResult>>();

		try {

			for (int i = 0; i < peptideList.size(); i++) {
				QuantifiedPeptideInterface pep1 = peptideList.get(i);

				if (i % 25 == 0) {
					log.info("Aligning Peptides " + i + "/" + peptideList.size());
				}

				for (int j = i + 1; j < peptideList.size(); j++) {
					QuantifiedPeptideInterface pep2 = peptideList.get(j);
					NWResult pepResults = NWAlign.needlemanWunsch(pep1.getSequence(), pep2.getSequence(), -11, -1);

					if ((pepResults.getFinalAlignmentScore() >= ProteinClusterQuantParameters.getInstance()
							.getFinalAlignmentScore()
							&& pepResults.getSequenceIdentity() >= ProteinClusterQuantParameters.getInstance()
									.getSequenceIdentity())
							&& pepResults.getMaxConsecutiveIdenticalAlignment() >= ProteinClusterQuantParameters
									.getInstance().getMaxConsecutiveIdenticalAlignment()) {
						// pepResults call to function to put result in map
						// pass map, result, call on both pep1
						gAM = putResultInMap(gAM, pepResults, gAS, pep1);
						gAM = putResultInMap(gAM, pepResults, gAS, pep2);

						// Set<QuantifiedProteinInterface> protSet =
						// pep1.getQuantifiedProteins();
						//
						// if ((pepResults.getFinalAlignmentScore() >= 30 &&
						// pepResults.getSequenceIdentity() >= 0.8)
						// && pepResults.getMaxConsecutiveIdenticalAlignment()
						// >= 6)
						// {
						// for (QuantifiedProteinInterface quantifiedProtein :
						// protSet) {
						// print p1 and pep1
						/*
						 * System.out.print(quantifiedProtein.getAccession() +
						 * "\t"); System.out.print(pep1.getSequence() + "\t");
						 * System.out.print("1" + "\t"); System.out.print("0" +
						 * "\t"); log.info("1");
						 */

						// print pep1 and pep2
						Output2.append(pep1.getSequence() + "\t");
						Output2.append(pep2.getSequence() + "\t");
						Double ratio1 = null;
						if (pep1 instanceof IsobaricQuantifiedPeptide) {
							ratio1 = ((IsobaricQuantifiedPeptide) pep1).getCountRatio(cond1, cond2);
						}
						if (ratio1 == Double.NEGATIVE_INFINITY) {
							Output2.append("NEG_INF" + "\t");
						} else if (ratio1 == Double.POSITIVE_INFINITY) {
							Output2.append("POS_INF" + "\t");
						} else {
							Output2.append(ratio1 + "\t");
						}
						Double ratio2 = null;
						if (pep2 instanceof IsobaricQuantifiedPeptide) {
							ratio2 = ((IsobaricQuantifiedPeptide) pep2).getCountRatio(cond1, cond2);
						}
						if (ratio2 == Double.NEGATIVE_INFINITY) {
							Output2.append("NEG_INF" + "\t");
						} else if (ratio2 == Double.POSITIVE_INFINITY) {
							Output2.append("POS_INF" + "\t");
						} else {
							Output2.append(ratio2 + "\t");
						}

						Output2.append(pepResults.getSequenceIdentity() + "\t");
						Output2.append("1" + "\n");
						Output2.flush();

						/*
						 * // prints p1 and pep2
						 * System.out.print(quantifiedProtein.getAccession() +
						 * "\t"); System.out.print(pep2.getSequence() + "\t");
						 * System.out.print("1" + "\t"); System.out.print("0" +
						 * "\t"); log.info("1");
						 */

						// }
					}

				}
			}
			return gAM;
		} finally

		{
			if (Output2 != null) {
				Output2.close();
			}
		}

	}

	public static Map<String, Set<NWResult>> putResultInMap(Map<String, Set<NWResult>> gAM, NWResult result,
			Set<NWResult> gAS, QuantifiedPeptideInterface pep) {
		String seq = pep.getSequence();

		// if gAM has the sequence, use that sequence to get the result
		if ((gAM.containsKey(seq))) {
			gAM.get(seq).add(result);
		}
		// if gAM does not have the sequence, make a new set of results, add the
		// results to this set, and then put the sequence and the set in gAM
		else {
			Set<NWResult> set = new HashSet<NWResult>();
			set.add(result);
			gAM.put(seq, set);
		}
		// return to alignPeptides
		return gAM;
	}

	// merges the clusters if there is a similar peptide pair between the two
	// clusters
	public static ProteinCluster mergeClusters(ProteinCluster cluster, ProteinCluster cluster2) {
		for (QuantifiedPeptideInterface peptide : cluster2.getPeptideSet()) {
			cluster.addPeptides(peptide);
		}

		for (QuantifiedProteinInterface protein : cluster2.getProteinSet()) {
			cluster.addProteins(protein);
		}
		return cluster;
	}

	public static CensusChroParser getCensusChroParser(File fastaFile, File inputFileFolder, String[] fileNames,
			Map<QuantCondition, QuantificationLabel> labelsByConditions, char[] enzymeArray, int missedCleavages,
			File uniprotReleasesFolder, String decoyRegexp) throws FileNotFoundException {
		List<Map<QuantCondition, QuantificationLabel>> list = new ArrayList<Map<QuantCondition, QuantificationLabel>>();
		for (String fileName : fileNames) {
			list.add(labelsByConditions);
		}
		return getCensusChroParser(fastaFile, inputFileFolder, fileNames, list, enzymeArray, missedCleavages,
				uniprotReleasesFolder, decoyRegexp);
	}

	public static CensusChroParser getCensusChroParser(File fastaFile, File inputFilefolder, String[] fileNames,
			List<Map<QuantCondition, QuantificationLabel>> labelsByConditions, char[] enzymeArray, int missedCleavages,
			File uniprotReleasesFolder, String decoyRegexp) throws FileNotFoundException {
		// Set parser (6 files) to peptides
		List<RemoteSSHFileReference> xmlFiles = new ArrayList<RemoteSSHFileReference>();

		for (String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.add(new RemoteSSHFileReference(inputXmlFile));
		}

		CensusChroParser parser = new CensusChroParser(xmlFiles, labelsByConditions);
		parser.addIonExclusion(IonSerieType.B, 1);
		parser.addIonExclusion(IonSerieType.Y, 1);
		// gets rid of decoys
		parser.setDecoyPattern(decoyRegexp);
		parser.setDistinguishModifiedPeptides(false);

		if (fastaFile != null) {
			DBIndexSearchParams defaultDBIndexParams = DBIndexInterface.getDefaultDBIndexParams(fastaFile);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeArr(enzymeArray);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeOffset(0);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeNocutResidues("");
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setMaxMissedCleavages(missedCleavages);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setH2OPlusProtonAdded(true);
			DBIndexInterface dbIndex = new DBIndexInterface(defaultDBIndexParams);
			parser.setDbIndex(dbIndex);
		}
		// change primary accession to the latest in Uniprot
		UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);

		Set<String> accessions = new HashSet<String>();
		accessions.addAll(parser.getProteinMap().keySet());
		log.info("Getting annotations from Uniprot...");
		Map<String, Entry> annotatedProteins = uplr.getAnnotatedProteins(null, accessions);

		for (String accession : accessions) {
			QuantifiedProteinInterface quantifiedProtein = parser.getProteinMap().get(accession);
			Entry entry = annotatedProteins.get(accession);
			if (entry != null && entry.getAccession() != null && !entry.getAccession().isEmpty()) {
				String primaryAccession = entry.getAccession().get(0);
				if (!accession.equals(primaryAccession) && !accession.contains(primaryAccession)) {
					log.info("Replacing Uniprot accession " + quantifiedProtein.getAccession() + " by "
							+ primaryAccession);
					quantifiedProtein.setAccession(primaryAccession);
					if (parser.getProteinMap().containsKey(primaryAccession)) {
						// there was already a protein with that
						// primaryAccession
						QuantifiedProteinInterface quantifiedProtein2 = parser.getProteinMap().get(primaryAccession);
						// merge quantifiedPRotein and quantifiedPRotein2
						mergeProteins(quantifiedProtein, quantifiedProtein2);
					}
					parser.getProteinMap().put(primaryAccession, quantifiedProtein);
				}
			} else {
				// // remove the protein because is obsolete
				// log.info(quantifiedProtein.getAccession());
				// parser.getProteinMap().remove(accession);
			}
		}

		return parser;
	}

	private static void mergeProteins(QuantifiedProteinInterface proteinReceiver,
			QuantifiedProteinInterface proteinDonor) {
		// PSMS
		for (QuantifiedPSMInterface psm : proteinDonor.getQuantifiedPSMs()) {
			proteinReceiver.addPSM(psm);
			psm.getQuantifiedProteins().remove(proteinDonor);
			psm.addQuantifiedProtein(proteinReceiver);
		}
		// Peptides
		for (QuantifiedPeptideInterface peptide : proteinDonor.getQuantifiedPeptides()) {
			proteinReceiver.addPeptide(peptide);
			peptide.getQuantifiedProteins().remove(proteinDonor);

		}
	}

	/**
	 * This will create a new protein which will contains all peptides and psms
	 * of the passed proteins. To be consistent, all psms and peptides of the
	 * passed proteins will be 'disconnected' from the original protein and will
	 * be connected to the merged one.
	 *
	 * @param proteins
	 * @return
	 */
	public static QuantifiedProteinInterface mergeProteins(Set<QuantifiedProteinInterface> proteins) {
		String consensusACC = Utils.getAccessionString(proteins);
		QuantifiedProteinInterface mergedProtein = new QuantifiedProteinFromCensusOut(consensusACC);
		// this is necessary because when passed in the constructor, it is being
		// parsed to only one:
		mergedProtein.setAccession(consensusACC);
		mergedProtein.setDescription(Utils.getDescriptionString(proteins));
		List<String> taxonomies = new ArrayList<String>();
		for (QuantifiedProteinInterface quantifiedProtein : proteins) {
			if (!taxonomies.contains(quantifiedProtein.getTaxonomy())) {
				taxonomies.add(quantifiedProtein.getTaxonomy());
			}
			// add peptides to new protein and
			// remove old protein from peptides
			final Set<QuantifiedPeptideInterface> quantifiedPeptides = quantifiedProtein.getQuantifiedPeptides();
			for (QuantifiedPeptideInterface quantifiedPeptide : quantifiedPeptides) {
				mergedProtein.addPeptide(quantifiedPeptide);
				quantifiedPeptide.getQuantifiedProteins().remove(quantifiedProtein);
			}

			// add psms to new protein and
			// remove old protein from psms
			final Set<QuantifiedPSMInterface> quantifiedPSMs = quantifiedProtein.getQuantifiedPSMs();
			for (QuantifiedPSMInterface quantifiedPSM : quantifiedPSMs) {
				mergedProtein.addPSM(quantifiedPSM);
				quantifiedPSM.getQuantifiedProteins().remove(quantifiedProtein);
			}
		}

		Collections.sort(taxonomies);
		StringBuilder sb = new StringBuilder();
		for (String taxonomy : taxonomies) {
			if (!"".equals(sb.toString())) {
				sb.append("-");
			}
			sb.append(taxonomy);
		}
		mergedProtein.setTaxonomy(sb.toString());

		return mergedProtein;
	}

	public static boolean shareAtLeastOnePeptide(QuantifiedProteinInterface protein1,
			QuantifiedProteinInterface protein2) {

		Set<QuantifiedPeptideInterface> peptides1 = protein1.getQuantifiedPeptides();
		for (QuantifiedPeptideInterface peptide : peptides1) {
			if (protein2.getQuantifiedPeptides().contains(peptide)) {
				return true;
			}
		}
		return false;
	}

	public static boolean shareAtLeastOnePeptideBySimilarity(QuantifiedProteinInterface protein1,
			QuantifiedProteinInterface protein2, Map<String, Set<NWResult>> gAM) {
		Set<QuantifiedPeptideInterface> peptides1 = protein1.getQuantifiedPeptides();

		for (QuantifiedPeptideInterface peptide1 : peptides1) {
			String peptideSeq1 = peptide1.getSequence();
			if (gAM.containsKey(peptideSeq1)) {
				Set<NWResult> alignments = gAM.get(peptideSeq1);

				for (NWResult nwResult : alignments) {
					String peptideSeq2 = nwResult.getSeq1();
					if (peptideSeq1.equals(peptideSeq2)) {
						peptideSeq2 = nwResult.getSeq2();
					}

					if (proteinContainsPeptideSequence(protein2, peptideSeq2)) {
						return true;
					}
				}
			}
		}
		return false;
	}

	private static Boolean proteinContainsPeptideSequence(QuantifiedProteinInterface protein, String peptideSeq) {
		Set<QuantifiedPeptideInterface> peptides = protein.getQuantifiedPeptides();
		for (QuantifiedPeptideInterface peptide : peptides) {
			if (peptide.getSequence().equals(peptideSeq)) {
				return true;
			}
		}
		return false;
	}

	public static Double findLargestNumber(List<Double> ratios) {

		Double largest = -Double.MAX_VALUE;
		for (Double ratio : ratios) {
			if (ratio > largest) {
				largest = ratio;
			}
		}
		return largest;
	}

	public static Double findSmallestNumber(List<Double> ratios) {

		Double smallest = Double.MAX_VALUE;
		for (Double ratio : ratios) {
			if (ratio < smallest) {
				smallest = ratio;
			}
		}
		return smallest;
	}

	// qTest threshold comparer
	public static boolean qTestThreshold(Double qResult1, Double qResult2, List<Double> ratios) {
		if (ratios.size() < 3) {
			return false;
		}

		if (ratios.size() == 3) {
			if (qResult1 > 0.97 || qResult2 > 0.97) {
				return true;
			}
		}
		return false;
	}

	// gets all the count ratios for protein
	public static List<Double> getCountRatio(QuantifiedProteinInterface protein1, QuantCondition cond1,
			QuantCondition cond2) {
		List<Double> ratioList = new ArrayList<Double>();
		Set<QuantifiedPeptideInterface> pepSet1 = protein1.getQuantifiedPeptides();
		for (QuantifiedPeptideInterface quantifiedPeptide : pepSet1) {
			if (quantifiedPeptide instanceof IsobaricQuantifiedPeptide) {
				IsobaricQuantifiedPeptide isoPeptide = (IsobaricQuantifiedPeptide) quantifiedPeptide;
				double countRatio = isoPeptide.getCountRatio(cond1, cond2);
				if (Double.isInfinite(countRatio) || Double.isNaN(countRatio)) {
					continue;
				}
				ratioList.add(countRatio);
			}
		}
		return ratioList;
	}

	// gets all the count ratios for protein
	public static List<Double> getAllCountRatio(QuantifiedProteinInterface protein1, QuantCondition cond1,
			QuantCondition cond2) {
		List<Double> ratioList = new ArrayList<Double>();
		Set<QuantifiedPeptideInterface> pepSet1 = protein1.getQuantifiedPeptides();
		for (QuantifiedPeptideInterface quantifiedPeptide : pepSet1) {
			if (quantifiedPeptide instanceof IsobaricQuantifiedPeptide) {
				IsobaricQuantifiedPeptide isoPeptide = (IsobaricQuantifiedPeptide) quantifiedPeptide;
				double countRatio = isoPeptide.getCountRatio(cond1, cond2);
				if (Double.isNaN(countRatio)) {
					continue;
				}
				ratioList.add(countRatio);
			}
		}
		return ratioList;
	}

	public static Ratio getSharedPepRatio(QuantifiedProteinInterface protein1, QuantifiedProteinInterface protein2,
			QuantCondition cond1, QuantCondition cond2) {
		final Set<QuantifiedPeptideInterface> sharedPeptides = getSharedPeptides(protein1, protein2, false);
		return getPepRatio(sharedPeptides, cond1, cond2, null);
	}

	public static Set<QuantifiedPeptideInterface> getSharedPeptides(QuantifiedProteinInterface protein1,
			QuantifiedProteinInterface protein2, boolean onlysharedByThisToProteins) {
		Set<QuantifiedPeptideInterface> peptides1 = protein1.getQuantifiedPeptides();
		Set<QuantifiedPeptideInterface> peptides2 = protein2.getQuantifiedPeptides();
		Set<QuantifiedPeptideInterface> ret = new HashSet<QuantifiedPeptideInterface>();
		for (QuantifiedPeptideInterface peptide1 : peptides1) {
			if (peptides2.contains(peptide1)) {
				boolean valid = true;
				if (!onlysharedByThisToProteins) {
					ret.add(peptide1);
				} else {
					// only is valid if the peptide only contains 2 proteins
					if (peptide1.getQuantifiedProteins().size() == 2) {
						ret.add(peptide1);
					}
				}
			}
		}
		return ret;
	}

	public static Ratio getPepRatio(Collection<QuantifiedPeptideInterface> peptides, QuantCondition cond1,
			QuantCondition cond2, String replicateName) {

		Ratio ratio = new Ratio();

		for (QuantifiedPeptideInterface peptide : peptides) {
			Set<QuantifiedPSMInterface> replicatePsms = new HashSet<QuantifiedPSMInterface>();
			for (QuantifiedPSMInterface quantifiedPSM : peptide.getQuantifiedPSMs()) {
				final String fileName = quantifiedPSM.getFileName();
				if (replicateName != null) {
					if (fileName.contains(replicateName)) {
						replicatePsms.add(quantifiedPSM);
					}
				} else {
					replicatePsms.add(quantifiedPSM);
				}
			}
			int numPSMs = replicatePsms.size();
			if (numPSMs == 0) {
				continue;
			}
			// get number of ions in one condition, and normalize by the
			// number of PSMs
			int peakCount1 = 0;
			for (QuantifiedPSMInterface quantifiedPSM : replicatePsms) {
				if (quantifiedPSM instanceof IsobaricQuantifiedPSM) {
					IsobaricQuantifiedPSM isoPSM = (IsobaricQuantifiedPSM) quantifiedPSM;
					if (isoPSM.getIonsByCondition().containsKey(cond1)) {
						peakCount1 += isoPSM.getIonsByCondition().get(cond1).size();
					}
				}
			}
			double normalizedPeakCount1 = peakCount1 * 1.0 / numPSMs;

			int peakCount2 = 0;
			for (QuantifiedPSMInterface quantifiedPSM : replicatePsms) {
				if (quantifiedPSM instanceof IsobaricQuantifiedPSM) {
					IsobaricQuantifiedPSM isoPeptide = (IsobaricQuantifiedPSM) quantifiedPSM;
					if (isoPeptide.getIonsByCondition().containsKey(cond2)) {
						peakCount2 += isoPeptide.getIonsByCondition().get(cond2).size();
					}
				}
			}
			double normalizedPeakCount2 = peakCount2 * 1.0 / numPSMs;

			ratio.addIonCount(cond1, normalizedPeakCount1);
			ratio.addIonCount(cond2, normalizedPeakCount2);
		}

		return ratio;
	}

	public static List<Double> getSharedPepRatioValues(QuantifiedProteinInterface protein1,
			QuantifiedProteinInterface protein2, QuantCondition cond1, QuantCondition cond2) {
		Set<QuantifiedPeptideInterface> peptides1 = protein1.getQuantifiedPeptides();
		Set<QuantifiedPeptideInterface> peptides2 = protein2.getQuantifiedPeptides();
		List<Double> ratios = new ArrayList<Double>();

		for (QuantifiedPeptideInterface peptide1 : peptides1) {
			if (peptides2.contains(peptide1)) {
				if (peptide1 instanceof QuantifiedPeptideInterface) {
					IsobaricQuantifiedPeptide isoPeptide = (IsobaricQuantifiedPeptide) peptide1;
					ratios.add(isoPeptide.getCountRatio(cond1, cond2));
				}
			}
		}
		return ratios;
	}

	/**
	 * Get the shared peptides between protein1 and protein2. Each peptide will
	 * be mapped to the proteins that are mapped.
	 *
	 * @param protein1
	 * @param protein2
	 * @param onlySharedByTheseTwoProteins
	 * @return
	 */
	public static Map<String, Set<QuantifiedPeptideInterface>> getSharedPeptidesByProtein(
			QuantifiedProteinInterface protein1, QuantifiedProteinInterface protein2,
			boolean onlySharedByTheseTwoProteins) {

		Set<QuantifiedPeptideInterface> peptides1 = protein1.getQuantifiedPeptides();
		Set<QuantifiedPeptideInterface> peptides2 = protein2.getQuantifiedPeptides();
		Set<QuantifiedPeptideInterface> totalPeptides = new HashSet<QuantifiedPeptideInterface>();
		totalPeptides.addAll(peptides1);
		totalPeptides.addAll(peptides2);
		Map<String, Set<QuantifiedPeptideInterface>> map = new HashMap<String, Set<QuantifiedPeptideInterface>>();
		for (QuantifiedPeptideInterface peptide : totalPeptides) {
			final Set<QuantifiedProteinInterface> quantifiedProteins = peptide.getQuantifiedProteins();
			if (quantifiedProteins.contains(protein1) && quantifiedProteins.contains(protein2)) {
				boolean include = false;
				if (onlySharedByTheseTwoProteins) {
					if (quantifiedProteins.size() == 2) {
						include = true;
					}
				} else {
					include = true;
				}
				if (include) {
					// peptide shared by protein1 and protein2
					for (QuantifiedProteinInterface quantifiedProtein : quantifiedProteins) {
						String proteinAccKey = quantifiedProtein.getAccession();
						if (map.containsKey(proteinAccKey)) {
							map.get(proteinAccKey).add(peptide);
						} else {
							Set<QuantifiedPeptideInterface> set = new HashSet<QuantifiedPeptideInterface>();
							set.add(peptide);
							map.put(proteinAccKey, set);
						}
					}

				}
			}
		}
		return map;
	}

	/**
	 * Get the shared peptides between protein1 and protein2. It is returned in
	 * a Map, in which each entry is a set of peptides that share the same
	 * proteins.<br>
	 * For example, one set will be the peptides shared ONLY by protein1 and
	 * protein2. Other set would be the peptides shared by protein1, protein2
	 * and a third protein3. Consequently, another set would be the peptides
	 * shared by protein1, protein2, protein3 and a fourth protein4, for
	 * example. And so on.<br>
	 * For each peptide set the key of the map will be the protein accessions
	 * where they belong sorted alphabetically.
	 *
	 * @param protein1
	 * @param protein2
	 * @param onlySharedByTheseTwoProteins
	 * @return
	 */
	public static Map<String, Set<QuantifiedPeptideInterface>> getSharedPeptidesMap(QuantifiedProteinInterface protein1,
			QuantifiedProteinInterface protein2, boolean onlySharedByTheseTwoProteins) {
		Map<String, Set<QuantifiedPeptideInterface>> map = new HashMap<String, Set<QuantifiedPeptideInterface>>();
		if (protein1 == null || protein2 == null) {
			return map;
		}
		Set<QuantifiedPeptideInterface> peptides1 = protein1.getQuantifiedPeptides();
		Set<QuantifiedPeptideInterface> peptides2 = protein2.getQuantifiedPeptides();
		Set<QuantifiedPeptideInterface> totalPeptides = new HashSet<QuantifiedPeptideInterface>();
		totalPeptides.addAll(peptides1);
		totalPeptides.addAll(peptides2);

		for (QuantifiedPeptideInterface peptide : totalPeptides) {
			final Set<QuantifiedProteinInterface> quantifiedProteins = peptide.getQuantifiedProteins();
			if (quantifiedProteins.contains(protein1) && quantifiedProteins.contains(protein2)) {
				boolean include = false;
				if (onlySharedByTheseTwoProteins) {
					if (quantifiedProteins.size() == 2) {
						include = true;
					}
				} else {
					include = true;
				}
				if (include) {
					// peptide shared by protein1 and protein2

					String proteinAccKey = Utils.getAccessionString(quantifiedProteins);
					if (map.containsKey(proteinAccKey)) {
						map.get(proteinAccKey).add(peptide);
					} else {
						Set<QuantifiedPeptideInterface> set = new HashSet<QuantifiedPeptideInterface>();
						set.add(peptide);
						map.put(proteinAccKey, set);
					}

				}
			}
		}
		return map;
	}

	/**
	 * explainwhatis doing
	 *
	 * @param protein1
	 * @param protein2
	 * @param cond1
	 * @param cond2
	 * @param uniquePepOnly
	 * @return
	 */
	public static Ratio getUniquePepRatio(QuantifiedProteinInterface protein1, QuantifiedProteinInterface protein2,
			QuantCondition cond1, QuantCondition cond2, boolean uniquePepOnly) {
		Set<QuantifiedPeptideInterface> peptides1 = protein1.getQuantifiedPeptides();
		Ratio ratio = new Ratio();

		for (QuantifiedPeptideInterface peptide1 : peptides1) {
			if (uniquePepOnly) {
				if (peptide1.getQuantifiedProteins().size() != 1) {
					continue;
				}

			} else {
				if (peptide1.getQuantifiedProteins().contains(protein2)) {
					continue;
				}
			}

			if (peptide1 instanceof IsobaricQuantifiedPeptide) {
				IsobaricQuantifiedPeptide isoPeptide = (IsobaricQuantifiedPeptide) peptide1;
				final int numPSMs = peptide1.getQuantifiedPSMs().size();
				// get number of ions in one condition, and normalize by the
				// number of PSMs
				int peakCount1 = 0;
				if (isoPeptide.getIonsByCondition().containsKey(cond1)) {
					peakCount1 = isoPeptide.getIonsByCondition().get(cond1).size();
				}
				double normalizedPeakCount1 = peakCount1 * 1.0 / numPSMs;
				int peakCount2 = 0;
				if (isoPeptide.getIonsByCondition().containsKey(cond2)) {
					peakCount2 = isoPeptide.getIonsByCondition().get(cond2).size();
				}
				double normalizedPeakCount2 = peakCount2 * 1.0 / numPSMs;
				ratio.addIonCount(cond1, normalizedPeakCount1);
				ratio.addIonCount(cond2, normalizedPeakCount2);
			}
		}
		return ratio;
	}

	/**
	 * explainwhatis doing
	 *
	 * @param protein1
	 * @param protein2
	 * @param cond1
	 * @param cond2
	 * @return
	 */
	public static List<Double> getUniquePepRatioValues(QuantifiedProteinInterface protein1,
			QuantifiedProteinInterface protein2, QuantCondition cond1, QuantCondition cond2, boolean uniquePepOnly) {
		Set<QuantifiedPeptideInterface> peptides1 = protein1.getQuantifiedPeptides();
		List<Double> ratios = new ArrayList<Double>();

		for (QuantifiedPeptideInterface peptide1 : peptides1) {
			if (uniquePepOnly) {
				if (peptide1.getQuantifiedProteins().size() == 1) {
					if (peptide1 instanceof IsobaricQuantifiedPeptide) {
						IsobaricQuantifiedPeptide isoPeptide = (IsobaricQuantifiedPeptide) peptide1;
						final double countRatio = isoPeptide.getCountRatio(cond1, cond2);
						ratios.add(countRatio);
					}
				} else {
					continue;
				}
			} else {
				if (peptide1.getQuantifiedProteins().contains(protein2)) {
					continue;
				} else {
					if (peptide1 instanceof IsobaricQuantifiedPeptide) {
						IsobaricQuantifiedPeptide isoPeptide = (IsobaricQuantifiedPeptide) peptide1;
						final double countRatio = isoPeptide.getCountRatio(cond1, cond2);
						ratios.add(countRatio);
					}
				}
			}
		}
		return ratios;
	}

	/**
	 * explainwhatis doing
	 *
	 * @param protein1
	 * @param protein2
	 * @param cond1
	 * @param cond2
	 * @return
	 */
	public static List<QuantifiedPeptideInterface> getUniquePeptides(QuantifiedProteinInterface protein1,
			QuantifiedProteinInterface protein2, boolean uniquePepOnly) {
		Set<QuantifiedPeptideInterface> peptides1 = protein1.getQuantifiedPeptides();
		List<QuantifiedPeptideInterface> ret = new ArrayList<QuantifiedPeptideInterface>();

		for (QuantifiedPeptideInterface peptide1 : peptides1) {
			if (uniquePepOnly) {
				if (peptide1.getQuantifiedProteins().size() == 1) {
					ret.add(peptide1);
				} else {
					continue;
				}
			} else {
				if (peptide1.getQuantifiedProteins().contains(protein2)) {
					continue;
				} else {
					ret.add(peptide1);
				}
			}
		}
		return ret;
	}

	/**
	 * Checks to see if the list has both ratios and infinities in it. Returns
	 * true if mixed. Returns false if only one kind.
	 *
	 * @param ratios
	 * @return
	 */
	public static boolean hasRatiosAndINF(List<Double> ratios) {
		int counterINF = 0;
		int counterRatio = 0;
		for (Double ratio : ratios) {
			if (ratio != null && Double.compare(ratio, Double.NaN) != 0) {
				if (Double.isInfinite(ratio)) {
					counterINF++;
				} else {
					counterRatio++;
				}
			} else {
				// means that one of the threeCompare is null!
			}
		}
		if (counterINF > 0 && counterRatio > 0) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Gets a list of all the infinities from the original list of ratios. If
	 * the ratio is null, it is ignored.
	 *
	 * @param ratios
	 * @return
	 */
	public static List<Double> getINFRatiosValues(List<Double> ratios) {
		List<Double> INFList = new ArrayList<Double>();
		for (Double ratio : ratios) {
			if (ratio != null && Double.compare(ratio, Double.NaN) != 0) {
				if (Double.isInfinite(ratio)) {
					INFList.add(ratio);
				}
			}
		}
		return INFList;
	}

	/**
	 * Gets a list of all the infinities from the original list of ratios. If
	 * the ratio is null, it is ignored.
	 *
	 * @param ratios
	 * @return
	 */
	public static List<Ratio> getINFRatios(List<Ratio> ratios, QuantCondition cond1, QuantCondition cond2) {
		List<Ratio> INFList = new ArrayList<Ratio>();
		for (Ratio ratio : ratios) {
			final double log2CountRatio = ratio.getLog2CountRatio(cond1, cond2);
			if (ratio != null && !Double.isNaN(log2CountRatio)) {
				if (Double.isInfinite(ratio.getCountRatio(cond1, cond2))) {
					INFList.add(ratio);
				}
			}
		}
		return INFList;
	}

	/**
	 * Gets a list of all the ratios (excludes INF) from the original list of
	 * ratios. If the ratio is null, it is ignored.
	 *
	 * @param ratios
	 * @return
	 */
	public static List<Double> getNonINFRatiosValues(List<Double> ratios) {
		List<Double> RatioList = new ArrayList<Double>();
		for (Double ratio : ratios) {
			if (ratio != null && Double.compare(ratio, Double.NaN) != 0) {
				if (!Double.isInfinite(ratio)) {
					RatioList.add(ratio);
				}
			}
		}
		return RatioList;
	}

	/**
	 * Gets a list of all the ratios (excludes INF) from the original list of
	 * ratios. If the ratio is null, it is ignored.
	 *
	 * @param ratios
	 * @return
	 */
	public static List<Ratio> getNonINFRatios(List<Ratio> ratios, QuantCondition cond1, QuantCondition cond2) {
		List<Ratio> RatioList = new ArrayList<Ratio>();
		for (Ratio ratio : ratios) {
			final double log2CountRatio = ratio.getLog2CountRatio(cond1, cond2);
			if (ratio != null && !Double.isNaN(log2CountRatio)) {
				if (!Double.isInfinite(log2CountRatio)) {
					RatioList.add(ratio);
				}
			}
		}
		return RatioList;
	}

	/**
	 * Checks to see if all the INF passed through are the same
	 *
	 * @param ratios
	 * @return Double, returns INF if all are INF, -INF if all are -INF, and
	 *         null if not the same.
	 */
	public static Double areAllINFValuesSame(List<Double> ratios) {
		int posCount = 0;
		int negCount = 0;
		for (Double ratio : ratios) {
			if (Double.compare(Double.POSITIVE_INFINITY, ratio) == 0) {
				posCount++;
			}
			if (Double.compare(Double.NEGATIVE_INFINITY, ratio) == 0) {
				negCount++;
			}
		}
		if (posCount > 0 && negCount > 0) {
			return null;
		}
		if (posCount > 0 && negCount == 0) {
			return Double.POSITIVE_INFINITY;
		}

		if (posCount == 0 && negCount > 0) {
			return Double.NEGATIVE_INFINITY;
		}
		return null;
	}

	/**
	 * Checks to see if all the INF passed through are the same
	 *
	 * @param ratios
	 * @return Double, returns INF if all are INF, -INF if all are -INF, and
	 *         null if not the same.
	 */
	public static Double areAllINFSame(List<Ratio> ratios, QuantCondition cond1, QuantCondition cond2) {
		int posCount = 0;
		int negCount = 0;
		for (Ratio ratio : ratios) {
			if (Double.compare(Double.POSITIVE_INFINITY, ratio.getLog2CountRatio(cond1, cond2)) == 0) {
				posCount++;
			}
			if (Double.compare(Double.NEGATIVE_INFINITY, ratio.getLog2CountRatio(cond1, cond2)) == 0) {
				negCount++;
			}
		}
		if (posCount > 0 && negCount > 0) {
			return null;
		}
		if (posCount > 0 && negCount == 0) {
			return Double.POSITIVE_INFINITY;
		}

		if (posCount == 0 && negCount > 0) {
			return Double.NEGATIVE_INFINITY;
		}
		return null;
	}

	/**
	 * Gets the number of ions (light or heavy) in the peptide
	 *
	 * @param peptide
	 * @return
	 */
	public static int getIonCount(QuantifiedPeptideInterface peptide, QuantificationLabel label) {
		int total = 0;
		if (peptide instanceof IsobaricQuantifiedPeptide) {
			Set<IsobaricQuantifiedPSM> psms = ((IsobaricQuantifiedPeptide) peptide).getIsobaricQuantifiedPSMs();

			for (IsobaricQuantifiedPSM quantifiedPSM : psms) {
				Set<Ion> ions = quantifiedPSM.getIonsByLabel(label);

				if (ions != null) {
					total = total + ions.size();
				}

			}
		}
		return total;
	}

	/**
	 * Gets the number of ions (light or heavy) in the peptide
	 *
	 * @param peptide
	 * @return
	 */
	public static int getIonCount(QuantifiedPeptideInterface peptide, QuantCondition condition) {
		int total = 0;
		if (peptide instanceof IsobaricQuantifiedPeptide) {
			Set<IsobaricQuantifiedPSM> psms = ((IsobaricQuantifiedPeptide) peptide).getIsobaricQuantifiedPSMs();

			for (IsobaricQuantifiedPSM quantifiedPSM : psms) {
				Set<Ion> ions = quantifiedPSM.getIonsByCondition().get(condition);
				if (ions != null) {
					total = total + ions.size();
				}

			}
		}
		return total;
	}

	/**
	 * Gets the number of ions (light or heavy) in the peptide
	 *
	 * @param peptide
	 * @return
	 */
	public static int getIonCount(QuantifiedPeptideInterface peptide) {
		Set<QuantifiedPeptideInterface> set = new HashSet<QuantifiedPeptideInterface>();
		set.add(peptide);
		return getIonCount(set);
	}

	public static List<QuantifiedPeptideInterface> getSortedPeptidesBySequence(
			Collection<QuantifiedPeptideInterface> peptides) {
		List<QuantifiedPeptideInterface> ret = new ArrayList<QuantifiedPeptideInterface>();
		ret.addAll(peptides);
		Collections.sort(ret, new Comparator<QuantifiedPeptideInterface>() {

			public int compare(QuantifiedPeptideInterface o1, QuantifiedPeptideInterface o2) {
				return o1.getSequence().compareTo(o2.getSequence());
			}
		});
		return ret;
	}

	public static List<QuantifiedProteinInterface> getSortedProteinsByAcc(
			Collection<QuantifiedProteinInterface> proteinsToSort) {
		List<QuantifiedProteinInterface> list = new ArrayList<QuantifiedProteinInterface>();
		list.addAll(proteinsToSort);
		Collections.sort(list, new Comparator<QuantifiedProteinInterface>() {
			public int compare(QuantifiedProteinInterface o1, QuantifiedProteinInterface o2) {
				return o1.getAccession().compareTo(o2.getAccession());
			}
		});
		return list;
	}

	/**
	 * Get a CVS list of protein accessions after sorting them alphabetically
	 *
	 * @param proteins
	 * @return
	 */
	public static String getAccessionString(Collection<QuantifiedProteinInterface> proteins) {

		StringBuilder sb = new StringBuilder();
		for (QuantifiedProteinInterface protein : Utils.getSortedProteinsByAcc(proteins)) {
			if (!"".equals(sb.toString()))
				sb.append(PROTEIN_ACC_SEPARATOR);
			sb.append(protein.getAccession());
		}
		return sb.toString();
	}

	/**
	 * Get a CVS list of peptide sequences after sorting them alphabetically by
	 * the sequence
	 *
	 * @param peptides
	 * @return
	 */
	public static String getPeptidesSequenceString(Collection<QuantifiedPeptideInterface> peptides) {

		StringBuilder sb = new StringBuilder();
		for (QuantifiedPeptideInterface peptide : Utils.getSortedPeptidesBySequence(peptides)) {
			if (!"".equals(sb.toString()))
				sb.append(",");
			sb.append(peptide.getSequence());
		}
		return sb.toString();
	}

	/**
	 * Gets a CSV list of protein descriptions, after sorting them by accession
	 * alphabetically
	 *
	 * @param proteins
	 * @return
	 */
	public static String getDescriptionString(Collection<QuantifiedProteinInterface> proteins) {
		StringBuilder sb = new StringBuilder();
		for (QuantifiedProteinInterface protein : Utils.getSortedProteinsByAcc(proteins)) {
			if (!"".equals(sb.toString())) {
				sb.append(PROTEIN_DESCRIPTION_SEPARATOR);
			}
			sb.append(protein.getDescription());
		}
		return sb.toString();
	}

	public static boolean individualProt(QuantifiedProteinInterface quantifiedProtein) {
		Set<QuantifiedPeptideInterface> peptides = quantifiedProtein.getQuantifiedPeptides();
		for (QuantifiedPeptideInterface quantifiedPeptide : peptides) {
			Set<QuantifiedProteinInterface> proteins = quantifiedPeptide.getQuantifiedProteins();
			if (proteins.size() != 1) {
				return false;
			}
		}
		return true;

	}

	public static String format(Double num) {
		if (num != null && !Double.isInfinite(num) && !Double.isNaN(num)) {
			return df.format(num);
		} else if (num != null && Double.isInfinite(num)) {
			if (Double.compare(num, Double.POSITIVE_INFINITY) == 0) {
				return "POS_INF";
			} else if (Double.compare(num, Double.NEGATIVE_INFINITY) == 0) {
				return "NEG_INF";
			}
		} else if (Double.isNaN(num)) {
			return "DNQ"; // means detected, not quantified
		}
		return String.valueOf(num);
	}

	public static String getPeptidesSequenceString(Map<String, Set<QuantifiedPeptideInterface>> peptideMap) {

		Set<QuantifiedPeptideInterface> set = new HashSet<QuantifiedPeptideInterface>();
		for (Set<QuantifiedPeptideInterface> quantifiedPeptideSet : peptideMap.values()) {
			set.addAll(quantifiedPeptideSet);
		}

		return getPeptidesSequenceString(set);
	}

	public static Set<QuantifiedPeptideInterface> getUniquePeptides(QuantifiedProteinInterface protein) {
		Set<QuantifiedPeptideInterface> ret = new HashSet<QuantifiedPeptideInterface>();
		for (QuantifiedPeptideInterface quantifiedPeptide : protein.getQuantifiedPeptides()) {
			if (quantifiedPeptide.getQuantifiedProteins().size() == 1) {
				ret.add(quantifiedPeptide);
			}
		}
		return ret;
	}

	// public static Double getAverageRatio(
	// Set<QuantifiedPeptideInterface> quantifiedPeptides, QuantCondition cond1,
	// QuantCondition cond2) {
	// List<Double> pepRatioList = new ArrayList<Double>();
	// for (QuantifiedPeptideInterface peps : quantifiedPeptides) {
	// Double ratio = peps.getCountRatio(cond1, cond2);
	// pepRatioList.add(ratio);
	// }
	// return ProteinPair.getConsensusRatio(pepRatioList);
	// }
	public static String getGeneNameString(Map<String, Entry> annotatedProteins, ProteinCluster cluster,
			Set<String> validTaxonomies, boolean onlyFirst) {
		final Set<QuantifiedProteinInterface> proteinSet = cluster.getProteinSet();

		return getGeneNameString(annotatedProteins, getProteinMap(proteinSet), validTaxonomies, onlyFirst);
	}

	public static String getGeneNameString(Map<String, Entry> annotatedProteins, QuantifiedProteinInterface protein,
			Set<String> validTaxonomies, boolean onlyFirst) {
		return getGeneNameString(annotatedProteins, getProteinMap(protein), validTaxonomies, onlyFirst);
	}

	public static Map<String, QuantifiedProteinInterface> getProteinMap(QuantifiedProteinInterface protein) {
		Set<QuantifiedProteinInterface> set = new HashSet<QuantifiedProteinInterface>();
		set.add(protein);
		return getProteinMap(set);
	}

	public static Map<String, QuantifiedProteinInterface> getProteinMap(
			Collection<QuantifiedProteinInterface> proteins) {
		Map<String, QuantifiedProteinInterface> map = new HashMap<String, QuantifiedProteinInterface>();
		for (QuantifiedProteinInterface quantifiedProtein : proteins) {
			final String accession = quantifiedProtein.getAccession();
			if (accession.contains(" ")) {
				final String[] split = accession.split(" ");
				for (String string : split) {
					if (string.length() == 1)
						continue;
					map.put(string, quantifiedProtein);
				}
			} else if (accession.contains("-")) {
				final String[] split = accession.split("-");
				for (String string : split) {
					if (string.length() == 1)
						continue;
					map.put(string, quantifiedProtein);
				}
			} else {
				map.put(accession, quantifiedProtein);
			}
		}
		return map;
	}

	public static String getGeneNameString(Map<String, Entry> annotatedProteins,
			Map<String, QuantifiedProteinInterface> proteinSet, Set<String> validTaxonomies, boolean onlyFirst) {

		Set<String> set = new HashSet<String>();
		for (String rawAcc : proteinSet.keySet()) {
			QuantifiedProteinInterface quantifiedProtein = proteinSet.get(rawAcc);

			List<String> accs = new ArrayList<String>();
			if (rawAcc.contains(" ")) {
				final String[] split = rawAcc.split(" ");
				for (String string : split) {
					accs.add(string);
				}
			} else {
				accs.add(rawAcc);
			}

			int index = 0;
			for (String acc : accs) {
				if (annotatedProteins.containsKey(acc)) {
					String taxon = getTaxonomy(acc, annotatedProteins.get(acc).getOrganism());
					if (taxon != null) {
						boolean valid = false;
						if (validTaxonomies != null && !validTaxonomies.isEmpty()) {
							for (String skipTaxonomy : validTaxonomies) {
								if (taxon.contains(skipTaxonomy)) {
									valid = true;
								}
							}
						} else {
							valid = true;
						}
						if (!valid) {
							continue;
						}
					}

					final String geneName = getGeneName(annotatedProteins.get(acc).getGene());
					set.add(geneName);
				} else {
					log.warn(acc + " not annotated");
				}
			}
		}
		List<String> list = new ArrayList<String>();
		list.addAll(set);
		Collections.sort(list);
		if (onlyFirst && !list.isEmpty()) {
			return list.get(0);
		}
		StringBuilder sb = new StringBuilder();
		if (!ProteinClusterQuantParameters.getInstance().isPrintOnlyFirstGene()) {
			for (String string : list) {
				if (!"".equals(sb.toString()))
					sb.append(",");
				sb.append(string);
			}
		} else {
			if (!list.isEmpty()) {
				sb.append(list.get(0));
			}
		}

		return sb.toString();
	}

	private static String getTaxonomy(String acc, OrganismType organism) {
		if (organism != null && acc != null) {
			if (organism.getName() != null) {
				for (OrganismNameType organismType : organism.getName()) {
					if (organismType.getType().equals("scientific")) {
						return organismType.getValue();
					}
				}
				if (!organism.getName().isEmpty()) {
					return organism.getName().get(0).getValue();
				}
			}
		}
		return null;
	}

	private static String getGeneName(List<GeneType> gene) {
		for (GeneType geneType : gene) {
			for (GeneNameType geneName : geneType.getName()) {
				if (geneName.getType().equals("primary"))
					return geneName.getValue();
			}
		}
		if (gene.isEmpty() || gene.get(0).getName().isEmpty()) {
			return "";
		}
		return gene.get(0).getName().get(0).getValue();
	}

	public static int getIonCount(Collection<QuantifiedPeptideInterface> peptides, QuantificationLabel label) {
		int ionCount = 0;
		if (peptides != null) {
			for (QuantifiedPeptideInterface quantifiedPeptide : peptides) {
				if (quantifiedPeptide instanceof IsobaricQuantifiedPeptide) {
					ionCount += Utils.getIonCount(quantifiedPeptide, label);
				}
			}
		}
		return ionCount;
	}

	public static int getIonCount(Collection<QuantifiedPeptideInterface> peptides, QuantCondition condition) {
		int ionCount = 0;
		if (peptides != null) {
			for (QuantifiedPeptideInterface quantifiedPeptide : peptides) {
				if (quantifiedPeptide instanceof IsobaricQuantifiedPeptide) {
					ionCount += Utils.getIonCount(quantifiedPeptide, condition);
				}
			}
		}
		return ionCount;
	}

	public static int getIonCount(Collection<QuantifiedPeptideInterface> peptides) {
		int ionCount = 0;
		if (peptides != null) {
			for (QuantifiedPeptideInterface quantifiedPeptide : peptides) {
				if (quantifiedPeptide instanceof IsobaricQuantifiedPeptide) {
					for (QuantificationLabel label : QuantificationLabel.values()) {
						ionCount += Utils.getIonCount(quantifiedPeptide, label);
					}
				}
			}
		}
		return ionCount;
	}

	public static Map<String, QuantifiedPeptideInterface> getPeptideMapFromClusters(Set<ProteinCluster> clusterSet) {
		Map<String, QuantifiedPeptideInterface> map = new HashMap<String, QuantifiedPeptideInterface>();
		for (ProteinCluster proteinCluster : clusterSet) {
			final Set<QuantifiedPeptideInterface> peptideSet = proteinCluster.getPeptideSet();
			for (QuantifiedPeptideInterface quantifiedPeptide : peptideSet) {
				final String peptideKey = quantifiedPeptide.getKey();
				if (map.containsKey(peptideKey)) {
					throw new IllegalArgumentException(
							"peptide key cannot be repeated. There is an incosistence b etween how Census parser is reading PSMs and creating peptides and how the peptideKey is created here");
				}
				map.put(peptideKey, quantifiedPeptide);
			}
		}
		return map;
	}

	public static boolean proteinsShareAllPeptides(QuantifiedProteinInterface protein1,
			QuantifiedProteinInterface protein2) {
		if (protein1.getQuantifiedPeptides().size() == protein2.getQuantifiedPeptides().size()) {
			for (QuantifiedPeptideInterface peptide : protein1.getQuantifiedPeptides()) {
				if (!protein2.getQuantifiedPeptides().contains(peptide)) {
					return false;
				}
			}
			return true;
		}
		return false;
	}

	public static Set<QuantifiedPSMInterface> getPSMs(Collection<QuantifiedPeptideInterface> peptides) {
		Set<QuantifiedPSMInterface> ret = new HashSet<QuantifiedPSMInterface>();
		for (QuantifiedPeptideInterface peptide : peptides) {
			ret.addAll(peptide.getQuantifiedPSMs());
		}
		return ret;
	}
}