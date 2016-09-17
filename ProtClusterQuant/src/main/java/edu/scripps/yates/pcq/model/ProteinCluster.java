package edu.scripps.yates.pcq.model;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.pcq.ProteinClusterQuantParameters;
import edu.scripps.yates.pcq.util.PCQUtils;
import edu.scripps.yates.utilities.alignment.nwalign.NWAlign;
import edu.scripps.yates.utilities.alignment.nwalign.NWResult;

public class ProteinCluster {
	private final static Logger log = Logger.getLogger(ProteinCluster.class);
	private static int proteinClusterCounter = 0;
	// Two sets for proteins and peptides
	private final Set<QuantifiedProteinInterface> individualQuantifiedProteinSet = new HashSet<QuantifiedProteinInterface>();
	private final Set<QuantifiedPeptideInterface> individualQuantifiedPeptideSet = new HashSet<QuantifiedPeptideInterface>();

	private final Set<PCQProteinNode> proteinNodes = new HashSet<PCQProteinNode>();
	private final Set<PCQPeptideNode> peptideNodes = new HashSet<PCQPeptideNode>();

	private final Set<ProteinPair> proteinPairs = new HashSet<ProteinPair>();
	private final Map<String, Set<NWResult>> alignmentsByPeptide = new HashMap<String, Set<NWResult>>();
	private Map<String, QuantifiedPeptideInterface> peptideMap;

	private final HashMap<String, PCQProteinNode> proteinNodesByProteinAcc = new HashMap<String, PCQProteinNode>();
	private final HashMap<String, PCQPeptideNode> peptideNodesByPeptideSequence = new HashMap<String, PCQPeptideNode>();

	private final int clusterID;

	public ProteinCluster() {
		clusterID = ++proteinClusterCounter;
	}

	private ProteinClusterQuantParameters getParams() {
		return ProteinClusterQuantParameters.getInstance();
	}

	/**
	 * By calling this function, the cluster will reorganize by collapsing the
	 * nodes appropiatelly and calculating the ratios accordingly
	 */
	public void reorganize() {
		// reset nodes
		resetNodes();

		// create nodes with the individual proteins and peptides

		// first create protein nodes
		createProteinNodes();

		// create nodes peptide nodes
		createPeptideNodes();

		// connect protein and peptide nodes
		connectProteinAndPeptideNodes();

	}

	private void connectProteinAndPeptideNodes() {
		for (PCQPeptideNode peptideNode : getPeptideNodes()) {
			final Set<QuantifiedPeptideInterface> quantifiedPeptides = peptideNode.getQuantifiedPeptides();
			for (QuantifiedPeptideInterface quantifiedPeptideInterface : quantifiedPeptides) {
				for (QuantifiedProteinInterface protein : quantifiedPeptideInterface.getQuantifiedProteins()) {

					PCQProteinNode proteinNode = proteinNodesByProteinAcc.get(protein.getAccession());
					// connect both
					proteinNode.addPeptideNode(peptideNode);
					peptideNode.addProteinNode(proteinNode);
				}
			}
		}
		// log.info(this);
	}

	private void createPeptideNodes() {
		// log.debug("Creating peptide nodes in cluster");
		List<QuantifiedPeptideInterface> peptides = new ArrayList<QuantifiedPeptideInterface>();
		peptides.addAll(individualQuantifiedPeptideSet);
		Set<String> set = new HashSet<String>();
		for (QuantifiedPeptideInterface pep : peptides) {
			if (set.contains(pep.getSequence())) {
				log.info("error");
			} else {
				set.add(pep.getSequence());
			}
		}
		if (peptides.size() > 1) {
			for (int i = 0; i < peptides.size(); i++) {
				QuantifiedPeptideInterface peptide1 = peptides.get(i);
				for (int j = i + 1; j < peptides.size(); j++) {
					QuantifiedPeptideInterface peptide2 = peptides.get(j);
					if (getParams().isCollapseIndistinguishablePeptides()
							&& PCQUtils.peptidesShareAllProteins(peptide1, peptide2)) {
						PCQPeptideNode peptideNode = null;
						PCQPeptideNode peptideNode2 = null;
						// if the protein was already associated with some
						// protein
						// node
						if (peptideNodesByPeptideSequence.containsKey(peptide1.getSequence())) {
							peptideNode = peptideNodesByPeptideSequence.get(peptide1.getSequence());
						}
						if (peptideNodesByPeptideSequence.containsKey(peptide2.getSequence())) {
							peptideNode2 = peptideNodesByPeptideSequence.get(peptide2.getSequence());
						}
						if (peptideNode == null && peptideNode2 == null) {
							peptideNode = new PCQPeptideNode(this, peptide1, peptide2);
						} else {
							if (peptideNode != null && peptideNode2 != null) {
								// remove both nodes from set
								peptideNodes.remove(peptideNode);
								peptideNodes.remove(peptideNode2);
								// merge nodes into peptideNode
								peptideNode = PCQUtils.mergePeptideNodes(peptideNode, peptideNode2);
							} else if (peptideNode == null && peptideNode2 != null) {
								peptideNode = peptideNode2;
							}
						}
						// add the two peptides to the peptide node
						peptideNode.addQuantifiedPeptide(peptide1);
						peptideNode.addQuantifiedPeptide(peptide2);
						// add to the set of nodes
						peptideNodes.add(peptideNode);
						// add to the map
						peptideNodesByPeptideSequence.put(peptide1.getSequence(), peptideNode);
						peptideNodesByPeptideSequence.put(peptide2.getSequence(), peptideNode);

					} else {
						// create a peptide node for each peptide separated
						// peptide node for psm1
						PCQPeptideNode peptideNode = null;
						if (peptideNodesByPeptideSequence.containsKey(peptide1.getSequence())) {
							peptideNode = peptideNodesByPeptideSequence.get(peptide1.getSequence());
						} else {
							peptideNode = new PCQPeptideNode(this, peptide1);
						}
						peptideNodes.add(peptideNode);
						peptideNodesByPeptideSequence.put(peptide1.getSequence(), peptideNode);

						// protein node for protein2
						PCQPeptideNode peptideNode2 = null;
						if (peptideNodesByPeptideSequence.containsKey(peptide2.getSequence())) {
							peptideNode2 = peptideNodesByPeptideSequence.get(peptide2.getSequence());
						} else {
							peptideNode2 = new PCQPeptideNode(this, peptide2);
						}
						peptideNodes.add(peptideNode2);
						peptideNodesByPeptideSequence.put(peptide2.getSequence(), peptideNode2);
					}
				}
			}
		} else {
			// only one peptide
			// create a peptide node for the peptide
			QuantifiedPeptideInterface peptide = peptides.iterator().next();
			PCQPeptideNode peptideNode = new PCQPeptideNode(this, peptide);
			peptideNodes.add(peptideNode);
			peptideNodesByPeptideSequence.put(peptide.getSequence(), peptideNode);
		}

	}

	private void createProteinNodes() {
		Map<String, Set<QuantifiedProteinInterface>> proteinMap = new HashMap<String, Set<QuantifiedProteinInterface>>();
		for (QuantifiedProteinInterface protein : individualQuantifiedProteinSet) {
			if (proteinMap.containsKey(protein.getAccession())) {
				proteinMap.get(protein.getAccession()).add(protein);
			} else {
				Set<QuantifiedProteinInterface> set = new HashSet<QuantifiedProteinInterface>();
				set.add(protein);
				proteinMap.put(protein.getAccession(), set);
			}
		}
		List<String> accList = new ArrayList<String>();
		accList.addAll(proteinMap.keySet());
		if (proteinMap.size() > 1) {
			for (int i = 0; i < accList.size(); i++) {
				final String acc1 = accList.get(i);
				Set<QuantifiedProteinInterface> proteins1 = proteinMap.get(acc1);
				for (int j = i + 1; j < accList.size(); j++) {
					final String acc2 = accList.get(j);
					Set<QuantifiedProteinInterface> proteins2 = proteinMap.get(acc2);

					if (getParams().isCollapseIndistinguishableProteins()
							&& PCQUtils.proteinsShareAllPeptides(proteins1, proteins2)) {
						PCQProteinNode proteinNode = null;
						PCQProteinNode proteinNode2 = null;
						// if the protein was already associated with some
						// protein
						// node
						if (proteinNodesByProteinAcc.containsKey(acc1)) {
							proteinNode = proteinNodesByProteinAcc.get(acc1);
						}
						if (proteinNodesByProteinAcc.containsKey(acc2)) {
							proteinNode2 = proteinNodesByProteinAcc.get(acc2);
						}
						if (proteinNode == null && proteinNode2 == null) {
							// if non of the proteins are in any protein node
							// yet, create a new one
							proteinNode = new PCQProteinNode(this, proteins1, proteins2);
						} else {
							// some of the proteins is assigned to some protein
							// node already
							if (proteinNode != null && proteinNode2 != null) {
								// both proteins are assigned to a node
								if (!proteinNode.equals(proteinNode2)) {
									// if each protein is assigned to a
									// different protein node
									// remove the second node
									proteinNodes.remove(proteinNode2);

									// merge protein nodes in the first protein
									// node
									PCQUtils.mergeProteinNodes(proteinNode, proteinNode2);

									// look for any other protein pointing to
									// proteinNode2 and assign it to
									// proteinNode1
									for (String acc : proteinNodesByProteinAcc.keySet()) {
										if (proteinNodesByProteinAcc.get(acc).equals(proteinNode2)) {
											proteinNodesByProteinAcc.put(acc, proteinNode);
										}
									}
								}
							} else if (proteinNode == null && proteinNode2 != null) {
								// protein 2 is assigned to a protein node, but
								// protein 1 is not. Use that node.
								proteinNode = proteinNode2;
							}
						}
						// add the two proteins to the protein node
						proteinNode.addProteins(proteins1);
						proteinNode.addProteins(proteins2);

						// add to the set of nodes
						proteinNodes.add(proteinNode);
						// add to the map
						proteinNodesByProteinAcc.put(acc1, proteinNode);
						proteinNodesByProteinAcc.put(acc2, proteinNode);

					} else {
						// create a protein node for each protein separated
						// protein node for protein1
						PCQProteinNode proteinNode = null;
						if (proteinNodesByProteinAcc.containsKey(acc1)) {
							proteinNode = proteinNodesByProteinAcc.get(acc1);
						} else {
							proteinNode = new PCQProteinNode(this, proteins1);
						}
						proteinNodes.add(proteinNode);
						proteinNodesByProteinAcc.put(acc1, proteinNode);

						// protein node for protein2
						PCQProteinNode proteinNode2 = null;
						if (proteinNodesByProteinAcc.containsKey(acc2)) {
							proteinNode2 = proteinNodesByProteinAcc.get(acc2);
						} else {
							proteinNode2 = new PCQProteinNode(this, proteins2);
						}
						proteinNodes.add(proteinNode2);
						proteinNodesByProteinAcc.put(acc2, proteinNode2);

					}
				}
			}
		} else {
			// only one protein
			// create a protein node for the protein
			// protein node for protein1
			final Collection<QuantifiedProteinInterface> proteins = proteinMap.values().iterator().next();
			PCQProteinNode proteinNode = new PCQProteinNode(this, proteins);
			proteinNodes.add(proteinNode);
			proteinNodesByProteinAcc.put(proteinMap.keySet().iterator().next(), proteinNode);

		}
		// log.debug(proteinNodes.size() + " protein nodes created in cluster");
		// for (PCQProteinNode proteinNode : proteinNodes) {
		// log.info(proteinNode);
		// }

	}

	private void resetNodes() {
		resetPeptideNodes();
		resetProteinNodes();
	}

	private void resetProteinNodes() {
		proteinNodes.clear();
		proteinNodesByProteinAcc.clear();
	}

	private void resetPeptideNodes() {
		peptideNodes.clear();
		peptideNodesByPeptideSequence.clear();
	}

	/**
	 * @return the proteinNodes
	 */
	public Set<PCQProteinNode> getProteinNodes() {
		return proteinNodes;
	}

	/**
	 * @return the peptideNodes
	 */
	public Set<PCQPeptideNode> getPeptideNodes() {
		return peptideNodes;
	}

	// Adds proteins to the protein sets
	public boolean addIndividualQuantifiedProtein(QuantifiedProteinInterface protein) {
		return individualQuantifiedProteinSet.add(protein);
	}

	// Adds peptides to the peptide sets
	public boolean addIndividualQuantifiedPeptide(QuantifiedPeptideInterface peptide) {
		return individualQuantifiedPeptideSet.add(peptide);
	}

	public NWResult[][] align() {
		// get peptides from clusters in form of set

		if (individualQuantifiedPeptideSet.size() == 1) {
			return null;
		}

		NWResult[][] results = new NWResult[individualQuantifiedPeptideSet.size()][individualQuantifiedPeptideSet
				.size()];

		List<QuantifiedPeptideInterface> pepList = new ArrayList<QuantifiedPeptideInterface>();
		pepList.addAll(individualQuantifiedPeptideSet);

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
		return individualQuantifiedProteinSet;
	}

	public Set<QuantifiedPeptideInterface> getPeptideSet() {
		return individualQuantifiedPeptideSet;
	}

	public void createPairs(Map<String, Set<NWResult>> gAM) {
		proteinPairs.clear();
		// if
		// (ProteinClusterQuantParameters.getInstance().isCollapseIndistinguishableProteins())
		// {
		// createPairsCollapsingIndistinguisibleProteins(gAM);
		// } else {
		// createPairsWithoutCollapsing(gAM);
		// }
		createProteinPairs(gAM);
	}

	private void createProteinPairs(Map<String, Set<NWResult>> gAM) {

		List<PCQProteinNode> proteinNodeList = new ArrayList<PCQProteinNode>();
		proteinNodeList.addAll(proteinNodes);

		for (int i = 0; i < proteinNodeList.size(); i++) {
			for (int j = i + 1; j < proteinNodeList.size(); j++) {
				PCQProteinNode proteinNode1 = proteinNodeList.get(i);

				PCQProteinNode proteinNode2 = proteinNodeList.get(j);

				// DO NOT CREATE PROTEIN PAIRS WITH DISCARDED PROTEIN NODES
				if (proteinNode1.isDiscarded() || proteinNode2.isDiscarded()) {
					continue;
				}
				if (PCQUtils.shareAtLeastOnePeptideNode(proteinNode1, proteinNode2) || (gAM != null
						&& PCQUtils.shareAtLeastOnePeptideBySimilarity(proteinNode1, proteinNode2, gAM))) {
					ProteinPair pair = new ProteinPair(proteinNode1, proteinNode2);
					// final String string = pair.toString();
					proteinPairs.add(pair);
				}
			}
		}
	}

	public Set<ProteinPair> getProteinPairs() {
		return proteinPairs;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("PROTN[");
		for (PCQProteinNode protein : PCQUtils.getSortedProteinNodesByAcc(proteinNodes)) {
			if (!"PROTN[".equals(sb.toString())) {
				sb.append(",");
			}
			sb.append(protein.getAccession());
		}
		sb.append("]NTORP");
		StringBuilder sb2 = new StringBuilder();
		sb2.append("PEPN[");
		for (PCQPeptideNode peptide : PCQUtils.getSortedPeptideNodesBySequence(peptideNodes)) {
			if (!"PEPN[".equals(sb2.toString())) {
				sb2.append(",");
			}
			sb2.append(peptide.getSequence());
		}
		sb2.append("]NPEP");
		final String string = sb.append("\t").append(sb2).toString();
		return string;
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
	public Set<PCQPeptideNode> getAlignedPeptides(PCQPeptideNode peptideNode) {

		if (peptideNode.getQuantifiedPeptides().size() > 1) {
			// If peptide nodes contain more than one different peptide,
			// then the aligments cannot be done, so there is something
			// wrong
			return Collections.emptySet();
		}
		QuantifiedPeptideInterface peptide = peptideNode.getQuantifiedPeptides().iterator().next();

		Map<String, QuantifiedPeptideInterface> peptideMap = getPeptideMap();
		Set<PCQPeptideNode> ret = new HashSet<PCQPeptideNode>();
		if (alignmentsByPeptide.containsKey(peptide.getSequence())) {
			final Set<NWResult> alignments = alignmentsByPeptide.get(peptide.getSequence());
			for (NWResult nwResult : alignments) {
				final QuantifiedPeptideInterface pep1 = peptideMap.get(nwResult.getSeq1());
				final QuantifiedPeptideInterface pep2 = peptideMap.get(nwResult.getSeq2());
				if (pep1.equals(peptide) && !pep2.equals(peptide)) {
					ret.add(peptideNodesByPeptideSequence.get(pep2.getSequence()));
				} else if (!pep1.equals(peptide) && pep2.equals(peptide)) {
					ret.add(peptideNodesByPeptideSequence.get(pep1.getSequence()));
				}
			}
		}
		return ret;
	}

	private Map<String, QuantifiedPeptideInterface> getPeptideMap() {
		if (peptideMap == null) {
			peptideMap = new HashMap<String, QuantifiedPeptideInterface>();
			for (QuantifiedPeptideInterface quantifiedPeptide : individualQuantifiedPeptideSet) {
				peptideMap.put(quantifiedPeptide.getSequence(), quantifiedPeptide);
			}
		}
		return peptideMap;
	}

	/**
	 * Get the {@link NWResult} from two specific peptides if available.
	 *
	 * @param peptideNode
	 * @param peptideNode2
	 * @return
	 */
	public NWResult getAlignmentResult(PCQPeptideNode peptideNode, PCQPeptideNode peptideNode2) {
		if (alignmentsByPeptide.containsKey(peptideNode.getSequence())) {
			final Set<NWResult> alignments = alignmentsByPeptide.get(peptideNode.getSequence());
			for (NWResult nwResult : alignments) {
				if (peptideNode.getSequence().equals(nwResult.getSeq1())) {
					if (peptideNode2.getSequence().equals(nwResult.getSeq2())) {
						return nwResult;
					}
				} else if (peptideNode.getSequence().equals(nwResult.getSeq2())) {
					if (peptideNode2.getSequence().equals(nwResult.getSeq1())) {
						return nwResult;
					}
				}
			}
		}
		return null;
	}

	// public void applyIonsThreshold(int ionsPerPeptideThreshold) {
	//
	// Set<QuantifiedPeptideInterface> peptidesToDelete = new
	// HashSet<QuantifiedPeptideInterface>();
	// // get all nodes and apply threshold
	// if (proteinPairs.isEmpty()) {
	// // if there is no protein pairs, there is only one protein but maybe
	// // more than one peptide
	// if
	// (ProteinClusterQuantParameters.getInstance().isCollapseIndistinguishablePeptides())
	// {
	// if (Utils.getIonCount(peptideSet) < ionsPerPeptideThreshold) {
	// peptidesToDelete.addAll(peptideSet);
	// }
	// } else {
	// for (QuantifiedPeptideInterface quantifiedPeptide : peptideSet) {
	// if (Utils.getIonCount(quantifiedPeptide) < ionsPerPeptideThreshold) {
	// peptidesToDelete.add(quantifiedPeptide);
	// }
	// }
	// }
	// } else {
	// for (ProteinPair proteinPair : proteinPairs) {
	// QuantifiedProteinInterface protein1 = proteinPair.getProtein1();
	// QuantifiedProteinInterface protein2 = proteinPair.getProtein2();
	// // peptides U1
	// List<List<QuantifiedPeptideInterface>> pep1DoubleList = new
	// ArrayList<List<QuantifiedPeptideInterface>>();
	// final List<QuantifiedPeptideInterface> peptidesU1 =
	// Utils.getUniquePeptides(protein1, protein2, true);
	// if
	// (ProteinClusterQuantParameters.getInstance().isCollapseIndistinguishablePeptides())
	// {
	// pep1DoubleList.add(peptidesU1);
	// } else {
	// for (QuantifiedPeptideInterface quantifiedPeptide : peptidesU1) {
	// List<QuantifiedPeptideInterface> list = new
	// ArrayList<QuantifiedPeptideInterface>();
	// list.add(quantifiedPeptide);
	// pep1DoubleList.add(list);
	// }
	// }
	//
	// // peptides for deletion
	// for (List<QuantifiedPeptideInterface> uniquePeptides_U1 : pep1DoubleList)
	// {
	// int ionCount = Utils.getIonCount(uniquePeptides_U1);
	// if (ionCount < ionsPerPeptideThreshold) {
	// peptidesToDelete.addAll(uniquePeptides_U1);
	// }
	// }
	//
	// // peptides U2
	// List<List<QuantifiedPeptideInterface>> pep2DoubleList = new
	// ArrayList<List<QuantifiedPeptideInterface>>();
	// final List<QuantifiedPeptideInterface> peptidesU2 =
	// Utils.getUniquePeptides(protein2, protein1, true);
	// if
	// (ProteinClusterQuantParameters.getInstance().isCollapseIndistinguishablePeptides())
	// {
	// pep2DoubleList.add(peptidesU2);
	// } else {
	// for (QuantifiedPeptideInterface quantifiedPeptide : peptidesU2) {
	// List<QuantifiedPeptideInterface> list = new
	// ArrayList<QuantifiedPeptideInterface>();
	// list.add(quantifiedPeptide);
	// pep2DoubleList.add(list);
	// }
	// }
	// for (List<QuantifiedPeptideInterface> uniquePeptides_U2 : pep2DoubleList)
	// {
	// int ionCount = Utils.getIonCount(uniquePeptides_U2);
	// if (ionCount < ionsPerPeptideThreshold) {
	// peptidesToDelete.addAll(uniquePeptides_U2);
	// }
	// }
	//
	// // S12 peptides shared by protein 1 and protein 2
	// final Map<String, Set<QuantifiedPeptideInterface>> sharedPeptidesMap_S12
	// = Utils
	// .getSharedPeptidesMap(protein1, protein2, false);
	// for (Set<QuantifiedPeptideInterface> peptidesS12 :
	// sharedPeptidesMap_S12.values()) {
	// List<Set<QuantifiedPeptideInterface>> pep12DoubleList = new
	// ArrayList<Set<QuantifiedPeptideInterface>>();
	// if
	// (ProteinClusterQuantParameters.getInstance().isCollapseIndistinguishablePeptides())
	// {
	// pep12DoubleList.add(peptidesS12);
	// } else {
	// for (QuantifiedPeptideInterface peptide : peptidesS12) {
	// Set<QuantifiedPeptideInterface> set = new
	// HashSet<QuantifiedPeptideInterface>();
	// set.add(peptide);
	// pep12DoubleList.add(set);
	// }
	// }
	// for (Set<QuantifiedPeptideInterface> sharedPeptides_S12 :
	// pep12DoubleList) {
	// int ionCount = Utils.getIonCount(sharedPeptides_S12);
	// if (ionCount < ionsPerPeptideThreshold) {
	// peptidesToDelete.addAll(sharedPeptides_S12);
	// }
	// }
	// }
	// }
	// }
	// // discard all peptides in peptidesToDelete
	// discardedPeptides += peptidesToDelete.size();
	// for (QuantifiedPeptideInterface quantifiedPeptide : peptidesToDelete) {
	// // unlink from proteins and PSMs
	// discardPeptide(quantifiedPeptide);
	// // remove from peptide set
	// peptideSet.remove(quantifiedPeptide);
	// }
	// // remove proteins in proteinSet not having peptides
	// final Iterator<QuantifiedProteinInterface> proteinsIterator =
	// proteinSet.iterator();
	// while (proteinsIterator.hasNext()) {
	// QuantifiedProteinInterface protein = proteinsIterator.next();
	// if (protein.getQuantifiedPeptides().isEmpty()) {
	// proteinsIterator.remove();
	// }
	// }
	// }

	public Set<String> getPeptideNodeKeys() {
		Set<String> ret = new HashSet<String>();
		for (PCQPeptideNode peptideNode : peptideNodes) {
			ret.add(peptideNode.getKey());
		}

		return ret;
	}

	/**
	 * Get the number of different individual peptides in the cluster
	 *
	 * @return
	 */
	public int getNumDifferentNonDiscardedIndividualPeptides() {
		Set<String> set = new HashSet<String>();
		for (PCQPeptideNode peptideNode : getNonDiscardedPeptideNodes()) {
			for (QuantifiedPeptideInterface peptide : peptideNode.getItemsInNode()) {
				set.add(peptide.getSequence());
			}
		}
		return set.size();
	}

	/**
	 * Get the number of different individual proteins in the cluster
	 *
	 * @return
	 */
	public int getNumDifferentNonDiscardedIndividualProteins() {
		Set<String> set = new HashSet<String>();
		for (PCQProteinNode proteinNode : getNonDiscardedProteinNodes()) {
			for (QuantifiedProteinInterface protein : proteinNode.getItemsInNode()) {
				set.add(protein.getAccession());
			}
		}
		return set.size();
	}

	/**
	 * @return the clusterID
	 */
	public int getClusterID() {
		return clusterID;
	}

	public PCQPeptideNode getPeptideNodeByKey(String peptideNodeID) {

		for (PCQPeptideNode peptideNode : peptideNodes) {
			if (peptideNode.getKey().equals(peptideNodeID)) {
				return peptideNode;
			}
		}

		return null;
	}

	public Set<PCQProteinNode> getNonDiscardedProteinNodes() {
		Set<PCQProteinNode> ret = new HashSet<PCQProteinNode>();
		for (PCQProteinNode proteinNode : getProteinNodes()) {
			if (!proteinNode.isDiscarded()) {
				ret.add(proteinNode);
			}
		}
		return ret;
	}

	public Set<PCQPeptideNode> getNonDiscardedPeptideNodes() {
		Set<PCQPeptideNode> ret = new HashSet<PCQPeptideNode>();
		for (PCQPeptideNode peptideNode : getPeptideNodes()) {
			if (!peptideNode.isDiscarded()) {
				ret.add(peptideNode);
			}
		}
		return ret;
	}

	public Set<QuantifiedProteinInterface> getNonDiscardedProteinSet() {
		Set<QuantifiedProteinInterface> ret = new HashSet<QuantifiedProteinInterface>();
		final Set<QuantifiedProteinInterface> proteinSet = getProteinSet();
		for (QuantifiedProteinInterface protein : proteinSet) {
			if (!protein.isDiscarded()) {
				ret.add(protein);
			}
		}
		return ret;
	}
}
