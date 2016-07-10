package edu.scripps.yates.proteinclusters;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.proteinclusters.util.ProteinClusterQuantParameters;
import edu.scripps.yates.proteinclusters.util.Utils;
import edu.scripps.yates.utilities.alignment.nwalign.NWAlign;
import edu.scripps.yates.utilities.alignment.nwalign.NWResult;

public class ProteinCluster {

	// Two sets for proteins and peptides
	private Set<QuantifiedProteinInterface> proteinSet = new HashSet<QuantifiedProteinInterface>();
	private final Set<QuantifiedPeptideInterface> peptideSet = new HashSet<QuantifiedPeptideInterface>();
	private final Set<ProteinPair> proteinPairs = new HashSet<ProteinPair>();
	private final Map<String, Set<NWResult>> alignmentsByPeptide = new HashMap<String, Set<NWResult>>();
	private Map<String, QuantifiedPeptideInterface> peptideMap;
	public static int discardedPeptides = 0;

	public ProteinCluster() {
	}

	// Adds proteins to the protein sets
	public void addProteins(QuantifiedProteinInterface protein) {
		getProteinSet().add(protein);
	}

	// Adds peptides to the peptide sets
	public void addPeptides(QuantifiedPeptideInterface peptide) {
		peptideSet.add(peptide);
	}

	public NWResult[][] align() {
		// get peptides from clusters in form of set

		if (peptideSet.size() == 1) {
			return null;
		}

		NWResult[][] results = new NWResult[peptideSet.size()][peptideSet.size()];

		List<QuantifiedPeptideInterface> pepList = new ArrayList<QuantifiedPeptideInterface>();
		pepList.addAll(peptideSet);

		// two for loops: compare two different peptides
		for (int i = 0; i < pepList.size(); i++) {
			QuantifiedPeptideInterface pep1 = pepList.get(i);

			for (int j = i + 1; j < pepList.size(); j++) {
				if (i == j) {
					// same peptide, move on
					continue;
				}

				// NWAlign (give back score)
				QuantifiedPeptideInterface pep2 = pepList.get(j);
				NWResult result = NWAlign.needlemanWunsch(pep1.getSequence(), pep2.getSequence(), -11, -1);
				results[i][j] = result;
			}

		}
		return results;
	}

	public Set<QuantifiedProteinInterface> getProteinSet() {
		return proteinSet;
	}

	public void setProteinSet(Set<QuantifiedProteinInterface> proteinSet) {
		this.proteinSet = proteinSet;
	}

	public Set<QuantifiedPeptideInterface> getPeptideSet() {
		return peptideSet;
	}

	public void createPairs(Map<String, Set<NWResult>> gAM) {
		proteinPairs.clear();
		List<QuantifiedProteinInterface> proteinList = new ArrayList<QuantifiedProteinInterface>();
		proteinList.addAll(proteinSet);

		for (int i = 0; i < proteinList.size(); i++) {
			for (int j = i + 1; j < proteinList.size(); j++) {
				QuantifiedProteinInterface protein1 = proteinList.get(i);
				QuantifiedProteinInterface protein2 = proteinList.get(j);
				if (Utils.shareAtLeastOnePeptide(protein1, protein2)
						|| Utils.shareAtLeastOnePeptideBySimilarity(protein1, protein2, gAM)) {
					ProteinPair pair = new ProteinPair(protein1, protein2, this);
					proteinPairs.add(pair);
				}
			}
		}

	}

	public void createPairsCollapsingIndistinguisibleProteins(Map<String, Set<NWResult>> gAM) {
		proteinPairs.clear();
		boolean keepLooping = true;
		while (keepLooping) {
			keepLooping = false;
			List<QuantifiedProteinInterface> proteinList = new ArrayList<QuantifiedProteinInterface>();
			proteinList.addAll(proteinSet);

			// new proteins merged
			Set<QuantifiedProteinInterface> newProteinsMerged = new HashSet<QuantifiedProteinInterface>();
			// indistinguisableProteinsToDelete
			Set<QuantifiedProteinInterface> indistinguisableProteinsToDelete = new HashSet<QuantifiedProteinInterface>();
			// look for indistinguisible proteins

			for (int i = 0; i < proteinList.size(); i++) {
				QuantifiedProteinInterface protein1 = proteinList.get(i);
				// if (protein1.getAccession().equals("B4MH50")) {
				// System.out.println(this);
				// }
				for (int j = i + 1; j < proteinList.size(); j++) {
					QuantifiedProteinInterface protein2 = proteinList.get(j);
					if (Utils.proteinsShareAllPeptides(protein1, protein2)) {
						keepLooping = true;
						indistinguisableProteinsToDelete.add(protein1);
						indistinguisableProteinsToDelete.add(protein2);
						// merge all proteins in one
						Set<QuantifiedProteinInterface> indistinguisableProteins2 = new HashSet<QuantifiedProteinInterface>();
						indistinguisableProteins2.add(protein1);
						indistinguisableProteins2.add(protein2);
						QuantifiedProteinInterface proteinMerged = Utils.mergeProteins(indistinguisableProteins2);
						newProteinsMerged.add(proteinMerged);
						break;
					}
				}
				if (keepLooping) {
					break;
				}
			}

			// remove from proteinSet the indistinguisable proteins
			for (QuantifiedProteinInterface quantifiedProtein : indistinguisableProteinsToDelete) {
				boolean removed = proteinSet.remove(quantifiedProtein);
				// System.out.println(quantifiedProtein + " removed=" +
				// removed);
			}
			// add the new proteins merged to the proteinSet
			for (QuantifiedProteinInterface quantifiedProtein : newProteinsMerged) {
				proteinSet.add(quantifiedProtein);
			}
		}

		List<QuantifiedProteinInterface> proteinList1 = new ArrayList<QuantifiedProteinInterface>();
		List<QuantifiedProteinInterface> proteinList2 = new ArrayList<QuantifiedProteinInterface>();
		proteinList1.addAll(proteinSet);
		proteinList2.addAll(proteinSet);
		for (int i = 0; i < proteinList1.size(); i++) {
			for (int j = i + 1; j < proteinList2.size(); j++) {
				QuantifiedProteinInterface protein1 = proteinList1.get(i);
				QuantifiedProteinInterface protein2 = proteinList2.get(j);
				if (Utils.shareAtLeastOnePeptide(protein1, protein2)
						|| Utils.shareAtLeastOnePeptideBySimilarity(protein1, protein2, gAM)) {

					ProteinPair pair = new ProteinPair(protein1, protein2, this);
					proteinPairs.add(pair);
				}
			}
		}
	}

	/**
	 * Returns true if at least one of the {@link ProteinPair} is inconsistent
	 * according to a QTest
	 *
	 * @param cond1
	 * @param cond2
	 * @return
	 */
	public boolean isInconsistenceWithQTest(QuantCondition cond1, QuantCondition cond2) {
		for (ProteinPair proteinPair : proteinPairs) {
			if (proteinPair.isInconsistenceWithQTest(cond1, cond2)) {
				return true;
			}
		}
		return false;
	}

	public Set<ProteinPair> getProteinPairs() {
		return proteinPairs;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("[");
		for (QuantifiedProteinInterface protein : Utils.getSortedProteinsByAcc(proteinSet)) {
			if (!"[".equals(sb.toString())) {
				sb.append(",");
			}
			sb.append(protein.getAccession());
		}
		sb.append("]");
		StringBuilder sb2 = new StringBuilder();
		sb2.append("[");
		for (QuantifiedPeptideInterface peptide : Utils.getSortedPeptidesBySequence(peptideSet)) {
			if (!"[".equals(sb2.toString())) {
				sb2.append(",");
			}
			sb2.append(peptide.getSequence());
		}
		sb2.append("]");
		return sb.toString() + "\t" + sb2.toString();
	}

	public void addAlignment(NWResult alignResult) {
		if (alignResult != null) {
			final String seq1 = alignResult.getSeq1();
			final String seq2 = alignResult.getSeq2();
			addAlignmentByPeptide(seq1, alignResult);
			addAlignmentByPeptide(seq2, alignResult);
		}

	}

	private void addAlignmentByPeptide(String seq1, NWResult alignResult) {
		if (alignResult != null && seq1 != null) {
			if (alignmentsByPeptide.containsKey(seq1)) {
				alignmentsByPeptide.get(seq1).add(alignResult);
			} else {
				Set<NWResult> set = new HashSet<NWResult>();
				set.add(alignResult);
				alignmentsByPeptide.put(seq1, set);
			}
		}
	}

	/**
	 * If align was performed, this will return the corresponding peptides
	 * aligned to a given one in the cluster.
	 *
	 * @param peptide
	 * @return
	 */
	public Set<QuantifiedPeptideInterface> getAlignedPeptides(QuantifiedPeptideInterface peptide) {
		Map<String, QuantifiedPeptideInterface> peptideMap = getPeptideMap();
		Set<QuantifiedPeptideInterface> ret = new HashSet<QuantifiedPeptideInterface>();
		if (alignmentsByPeptide.containsKey(peptide.getSequence())) {
			final Set<NWResult> alignments = alignmentsByPeptide.get(peptide.getSequence());
			for (NWResult nwResult : alignments) {
				final QuantifiedPeptideInterface pep1 = peptideMap.get(nwResult.getSeq1());
				final QuantifiedPeptideInterface pep2 = peptideMap.get(nwResult.getSeq2());
				if (pep1.equals(peptide) && !pep2.equals(peptide)) {
					ret.add(pep2);
				} else if (!pep1.equals(peptide) && pep2.equals(peptide)) {
					ret.add(pep1);
				}
			}
		}
		return ret;
	}

	private Map<String, QuantifiedPeptideInterface> getPeptideMap() {
		if (peptideMap == null) {
			peptideMap = new HashMap<String, QuantifiedPeptideInterface>();
			for (QuantifiedPeptideInterface quantifiedPeptide : peptideSet) {
				peptideMap.put(quantifiedPeptide.getSequence(), quantifiedPeptide);
			}
		}
		return peptideMap;
	}

	/**
	 * Get the {@link NWResult} from two specific peptides if available.
	 *
	 * @param quantifiedPeptide
	 * @param quantifiedPeptide2
	 * @return
	 */
	public NWResult getAlignmentResult(QuantifiedPeptideInterface quantifiedPeptide,
			QuantifiedPeptideInterface quantifiedPeptide2) {
		if (alignmentsByPeptide.containsKey(quantifiedPeptide.getSequence())) {
			final Set<NWResult> alignments = alignmentsByPeptide.get(quantifiedPeptide.getSequence());
			for (NWResult nwResult : alignments) {
				if (quantifiedPeptide.getSequence().equals(nwResult.getSeq1())) {
					if (quantifiedPeptide2.getSequence().equals(nwResult.getSeq2())) {
						return nwResult;
					}
				} else if (quantifiedPeptide.getSequence().equals(nwResult.getSeq2())) {
					if (quantifiedPeptide2.getSequence().equals(nwResult.getSeq1())) {
						return nwResult;
					}
				}
			}
		}
		return null;
	}

	public void applyIonsThreshold(int ionsPerPeptideThreshold) {

		Set<QuantifiedPeptideInterface> peptidesToDelete = new HashSet<QuantifiedPeptideInterface>();
		// get all nodes and apply threshold
		if (proteinPairs.isEmpty()) {
			// if there is no protein pairs, there is only one protein but maybe
			// more than one peptide
			if (ProteinClusterQuantParameters.getInstance().isCollapseIndistinguishablePeptides()) {
				if (Utils.getIonCount(peptideSet) < ionsPerPeptideThreshold) {
					peptidesToDelete.addAll(peptideSet);
				}
			} else {
				for (QuantifiedPeptideInterface quantifiedPeptide : peptideSet) {
					if (Utils.getIonCount(quantifiedPeptide) < ionsPerPeptideThreshold) {
						peptidesToDelete.add(quantifiedPeptide);
					}
				}
			}
		} else {
			for (ProteinPair proteinPair : proteinPairs) {
				QuantifiedProteinInterface protein1 = proteinPair.getProtein1();
				QuantifiedProteinInterface protein2 = proteinPair.getProtein2();
				// peptides U1
				List<List<QuantifiedPeptideInterface>> pep1DoubleList = new ArrayList<List<QuantifiedPeptideInterface>>();
				final List<QuantifiedPeptideInterface> peptidesU1 = Utils.getUniquePeptides(protein1, protein2, true);
				if (ProteinClusterQuantParameters.getInstance().isCollapseIndistinguishablePeptides()) {
					pep1DoubleList.add(peptidesU1);
				} else {
					for (QuantifiedPeptideInterface quantifiedPeptide : peptidesU1) {
						List<QuantifiedPeptideInterface> list = new ArrayList<QuantifiedPeptideInterface>();
						list.add(quantifiedPeptide);
						pep1DoubleList.add(list);
					}
				}

				// peptides for deletion
				for (List<QuantifiedPeptideInterface> uniquePeptides_U1 : pep1DoubleList) {
					int ionCount = Utils.getIonCount(uniquePeptides_U1);
					if (ionCount < ionsPerPeptideThreshold) {
						peptidesToDelete.addAll(uniquePeptides_U1);
					}
				}

				// peptides U2
				List<List<QuantifiedPeptideInterface>> pep2DoubleList = new ArrayList<List<QuantifiedPeptideInterface>>();
				final List<QuantifiedPeptideInterface> peptidesU2 = Utils.getUniquePeptides(protein2, protein1, true);
				if (ProteinClusterQuantParameters.getInstance().isCollapseIndistinguishablePeptides()) {
					pep2DoubleList.add(peptidesU2);
				} else {
					for (QuantifiedPeptideInterface quantifiedPeptide : peptidesU2) {
						List<QuantifiedPeptideInterface> list = new ArrayList<QuantifiedPeptideInterface>();
						list.add(quantifiedPeptide);
						pep2DoubleList.add(list);
					}
				}
				for (List<QuantifiedPeptideInterface> uniquePeptides_U2 : pep2DoubleList) {
					int ionCount = Utils.getIonCount(uniquePeptides_U2);
					if (ionCount < ionsPerPeptideThreshold) {
						peptidesToDelete.addAll(uniquePeptides_U2);
					}
				}

				// S12 peptides shared by protein 1 and protein 2
				final Map<String, Set<QuantifiedPeptideInterface>> sharedPeptidesMap_S12 = Utils
						.getSharedPeptidesMap(protein1, protein2, false);
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
						int ionCount = Utils.getIonCount(sharedPeptides_S12);
						if (ionCount < ionsPerPeptideThreshold) {
							peptidesToDelete.addAll(sharedPeptides_S12);
						}
					}
				}

			}
		}
		// discard all peptides in peptidesToDelete
		discardedPeptides += peptidesToDelete.size();
		for (QuantifiedPeptideInterface quantifiedPeptide : peptidesToDelete) {
			// unlink from proteins and PSMs
			discardPeptide(quantifiedPeptide);
			// remove from peptide set
			peptideSet.remove(quantifiedPeptide);
		}
		// remove proteins in proteinSet not having peptides
		final Iterator<QuantifiedProteinInterface> proteinsIterator = proteinSet.iterator();
		while (proteinsIterator.hasNext()) {
			QuantifiedProteinInterface protein = proteinsIterator.next();
			if (protein.getQuantifiedPeptides().isEmpty()) {
				proteinsIterator.remove();
			}
		}
	}

	public Set<String> getPeptideNodeKeys() {
		Set<String> ret = new HashSet<String>();

		if (proteinPairs.isEmpty()) {
			if (ProteinClusterQuantParameters.getInstance().isCollapseIndistinguishablePeptides()) {
				ret.add(Utils.getPeptidesSequenceString(peptideSet));
			} else {
				for (QuantifiedPeptideInterface peptide : peptideSet) {
					ret.add(peptide.getSequence());
				}
			}
		} else {
			for (ProteinPair proteinPair : proteinPairs) {
				QuantifiedProteinInterface protein1 = proteinPair.getProtein1();
				QuantifiedProteinInterface protein2 = proteinPair.getProtein2();
				// peptides U1
				final List<QuantifiedPeptideInterface> peptidesU1 = Utils.getUniquePeptides(protein1, protein2, true);
				if (!peptidesU1.isEmpty()) {
					if (ProteinClusterQuantParameters.getInstance().isCollapseIndistinguishablePeptides()) {
						ret.add(Utils.getPeptidesSequenceString(peptidesU1));
					} else {
						for (QuantifiedPeptideInterface quantifiedPeptide : peptidesU1) {
							ret.add(quantifiedPeptide.getSequence());
						}
					}
				}
				// peptides U2
				final List<QuantifiedPeptideInterface> peptidesU2 = Utils.getUniquePeptides(protein2, protein1, true);
				if (!peptidesU2.isEmpty()) {
					if (ProteinClusterQuantParameters.getInstance().isCollapseIndistinguishablePeptides()) {
						ret.add(Utils.getPeptidesSequenceString(peptidesU2));
					} else {
						for (QuantifiedPeptideInterface quantifiedPeptide : peptidesU2) {
							ret.add(quantifiedPeptide.getSequence());
						}
					}
				}
				// S12 peptides shared by protein 1 and protein 2
				final Map<String, Set<QuantifiedPeptideInterface>> sharedPeptidesMap_S12 = Utils
						.getSharedPeptidesMap(protein1, protein2, false);
				if (!sharedPeptidesMap_S12.isEmpty()) {
					for (Set<QuantifiedPeptideInterface> peptidesS12 : sharedPeptidesMap_S12.values()) {
						if (ProteinClusterQuantParameters.getInstance().isCollapseIndistinguishablePeptides()) {
							ret.add(Utils.getPeptidesSequenceString(peptidesS12));
						} else {
							for (QuantifiedPeptideInterface peptide : peptidesS12) {
								ret.add(peptide.getSequence());
							}
						}
					}
				}
			}
		}
		return ret;
	}

	private void discardPeptide(QuantifiedPeptideInterface quantifiedPeptide) {

		// remove this peptide from its proteins
		final Set<QuantifiedProteinInterface> quantifiedProteins = quantifiedPeptide.getQuantifiedProteins();
		for (QuantifiedProteinInterface quantifiedProtein : quantifiedProteins) {
			quantifiedProtein.getQuantifiedPeptides().remove(quantifiedPeptide);
		}
		// get all psms of that peptide
		final Set<QuantifiedPSMInterface> quantifiedPSMs = quantifiedPeptide.getQuantifiedPSMs();
		for (QuantifiedPSMInterface quantifiedPSM : quantifiedPSMs) {
			// remove this psms from its proteins
			final Iterator<QuantifiedProteinInterface> quantifiedProteins2 = quantifiedPSM.getQuantifiedProteins()
					.iterator();
			while (quantifiedProteins2.hasNext()) {
				final QuantifiedProteinInterface quantifiedProtein = quantifiedProteins2.next();
				quantifiedProtein.getQuantifiedPSMs().remove(quantifiedPSM);
			}
			// remove this psm from teh parser
			// parser.getPSMMap().remove(quantifiedPSM.getPSMIdentifier());

		}

	}

	/**
	 * Get the number of individual proteins in the cluster. In case of
	 * collapsing the proteins that share all peptides, it detects them and
	 * count them even if they are represented in a merged
	 * QuantifiedProteinInterface object
	 *
	 * @return
	 */
	public int getNumIndividualProteins() {
		int num = 0;
		for (QuantifiedProteinInterface protein : proteinSet) {
			final String accession = protein.getAccession();
			// if accession contains an space, it means that it represents to
			// more than one protein
			if (accession.contains(" ")) {
				num += accession.split(" ").length;
			} else {
				num++;
			}
		}
		return num;
	}

	/**
	 * Get the set of {@link QuantifiedPeptideInterface} in the cluster.
	 *
	 * @param excludeUniquePeptides
	 *            if true, the unique peptides (belonging to only one protein)
	 *            will not be returned.
	 * @return
	 */
	public Collection<QuantifiedPeptideInterface> getPeptideSet(boolean excludeUniquePeptides) {
		final Set<QuantifiedPeptideInterface> wholePeptideSet = getPeptideSet();
		if (!excludeUniquePeptides) {
			return wholePeptideSet;
		}
		Set<QuantifiedPeptideInterface> ret = new HashSet<QuantifiedPeptideInterface>();
		for (QuantifiedPeptideInterface quantifiedPeptide : wholePeptideSet) {
			if (quantifiedPeptide.getQuantifiedProteins().size() > 1) {
				ret.add(quantifiedPeptide);
			}
		}
		return ret;
	}
}
