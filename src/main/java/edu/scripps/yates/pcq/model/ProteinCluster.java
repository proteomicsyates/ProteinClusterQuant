package edu.scripps.yates.pcq.model;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.census.read.model.IsobaricQuantifiedPeptide;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.pcq.filter.PCQFilter;
import edu.scripps.yates.pcq.params.ProteinClusterQuantParameters;
import edu.scripps.yates.pcq.util.PCQUtils;
import edu.scripps.yates.pcq.xgmml.util.AlignedPeptides;
import edu.scripps.yates.pcq.xgmml.util.AlignmentSet;
import edu.scripps.yates.utilities.alignment.nwalign.NWAlign;
import edu.scripps.yates.utilities.alignment.nwalign.NWResult;
import edu.scripps.yates.utilities.sequence.PositionInPeptide;
import edu.scripps.yates.utilities.sequence.PositionInProtein;
import edu.scripps.yates.utilities.strings.StringUtils;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class ProteinCluster {
	private final static Logger log = Logger.getLogger(ProteinCluster.class);
	private static int proteinClusterCounter = 0;
	// Two sets for proteins and peptides
	private final Set<QuantifiedProteinInterface> individualQuantifiedProteinSet = new THashSet<QuantifiedProteinInterface>();
	private final Set<QuantifiedPeptideInterface> individualQuantifiedPeptideSet = new THashSet<QuantifiedPeptideInterface>();

	private final Set<PCQProteinNode> proteinNodes = new THashSet<PCQProteinNode>();
	private final Set<PCQPeptideNode> peptideNodes = new THashSet<PCQPeptideNode>();

	private final Set<ProteinPair> proteinPairs = new THashSet<ProteinPair>();
	// proteinPairs with discarded proteins
	private final Set<ProteinPair> discardedProteinPairs = new THashSet<ProteinPair>();
	private final AlignmentSet alignmentSet = new AlignmentSet();

	private final Map<String, PCQProteinNode> proteinNodesByProteinKey = new THashMap<String, PCQProteinNode>();
	private final Map<String, PCQPeptideNode> peptideNodesByPeptideNodeKey = new THashMap<String, PCQPeptideNode>();

	private final int clusterID;

	public ProteinCluster() {
		clusterID = ++proteinClusterCounter;
	}

	private ProteinClusterQuantParameters getParams() {
		return ProteinClusterQuantParameters.getInstance();
	}

	/**
	 * By calling this function, the cluster will reorganize by collapsing the
	 * nodes properly and calculating the ratios accordingly
	 */
	public void createNodes() {
		// reset nodes
		// removes all nodes in the cluster and clear all maps
		resetNodes();

		// create nodes with the individual proteins and peptides

		// first create protein nodes
		createProteinNodes();

		// create nodes peptide nodes
		if (getParams().isCollapseBySites()) {
			createPeptideNodesByQuantifiedSites(getParams().getAaQuantified());
		} else {
			createPeptideNodes();
		}

		// connect protein and peptide nodes
		connectProteinAndPeptideNodes();

	}

	private void connectProteinAndPeptideNodes() {

		final Set<PCQPeptideNode> peptideNodes2 = getPeptideNodes();
		for (PCQPeptideNode peptideNode : peptideNodes2) {

			final Set<QuantifiedPeptideInterface> quantifiedPeptides = peptideNode.getQuantifiedPeptides();
			for (QuantifiedPeptideInterface peptide : quantifiedPeptides) {
				for (QuantifiedProteinInterface protein : peptide.getQuantifiedProteins()) {
					PCQProteinNode proteinNode = proteinNodesByProteinKey.get(protein.getKey());
					if (proteinNode == null) {
						log.info(this);
						boolean found = false;
						for (QuantifiedProteinInterface protein2 : getProteinSet()) {
							if (protein2.getKey().equals(protein.getKey())) {
								found = true;
							}
						}
						if (!found) {
							log.info(protein + " found " + found);
						}
					}
					// if (peptide.getPtms().isEmpty()) {
					// connect both
					proteinNode.addPeptideNode(peptideNode);
					peptideNode.addProteinNode(proteinNode);

				}

			}
		}

	}

	private void createPeptideNodes() {
		// log.debug("Creating peptide nodes in cluster");
		List<QuantifiedPeptideInterface> peptides = new ArrayList<QuantifiedPeptideInterface>();
		peptides.addAll(individualQuantifiedPeptideSet);
		Set<String> set = new THashSet<String>();
		for (QuantifiedPeptideInterface pep : peptides) {
			if (set.contains(pep.getKey())) {
				log.info("Inconsistency error. 2 peptides with the same key cannot exist!");
			} else {
				set.add(pep.getKey());
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
						if (peptideNodesByPeptideNodeKey.containsKey(peptide1.getKey())) {
							peptideNode = peptideNodesByPeptideNodeKey.get(peptide1.getKey());
						}
						if (peptideNodesByPeptideNodeKey.containsKey(peptide2.getKey())) {
							peptideNode2 = peptideNodesByPeptideNodeKey.get(peptide2.getKey());
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
						peptideNodesByPeptideNodeKey.put(peptide1.getKey(), peptideNode);
						peptideNodesByPeptideNodeKey.put(peptide2.getKey(), peptideNode);

					} else {
						// create a peptide node for each peptide separated
						// peptide node for psm1
						PCQPeptideNode peptideNode = null;
						if (peptideNodesByPeptideNodeKey.containsKey(peptide1.getKey())) {
							peptideNode = peptideNodesByPeptideNodeKey.get(peptide1.getKey());
						} else {
							peptideNode = new PCQPeptideNode(this, peptide1);
						}
						peptideNode.addQuantifiedPeptide(peptide1);
						peptideNodes.add(peptideNode);
						peptideNodesByPeptideNodeKey.put(peptide1.getKey(), peptideNode);

						// protein node for protein2
						PCQPeptideNode peptideNode2 = null;
						if (peptideNodesByPeptideNodeKey.containsKey(peptide2.getKey())) {
							peptideNode2 = peptideNodesByPeptideNodeKey.get(peptide2.getKey());
						} else {
							peptideNode2 = new PCQPeptideNode(this, peptide2);
						}
						peptideNode2.addQuantifiedPeptide(peptide2);
						peptideNodes.add(peptideNode2);
						peptideNodesByPeptideNodeKey.put(peptide2.getKey(), peptideNode2);
					}
				}
			}
		} else {
			// only one peptide
			// create a peptide node for the peptide
			QuantifiedPeptideInterface peptide = peptides.iterator().next();
			PCQPeptideNode peptideNode = new PCQPeptideNode(this, peptide);
			peptideNodes.add(peptideNode);
			peptideNodesByPeptideNodeKey.put(peptide.getKey(), peptideNode);
		}

	}

	/**
	 * Iterates all peptides in the cluster over all others and creates a node
	 * for all peptides that have the same key, which, in this case is the
	 * position of the quantifiedAA in the protein.
	 * 
	 * @param quantifiedAAs
	 */
	private void createPeptideNodesByQuantifiedSites(char[] quantifiedAAs) {
		UniprotProteinLocalRetriever uplr = PCQUtils
				.getUniprotProteinLocalRetrieverByFolder(getParams().getUniprotReleasesFolder());
		List<QuantifiedPeptideInterface> peptides = new ArrayList<QuantifiedPeptideInterface>();
		peptides.addAll(individualQuantifiedPeptideSet);
		Set<String> set = new THashSet<String>();
		for (QuantifiedPeptideInterface pep : peptides) {
			if (set.contains(pep.getKey())) {
				log.warn("Inconsistency error. 2 peptides with the same key cannot exist!");
			} else {
				set.add(pep.getKey());
			}
		}
		Set<String> peptideSequencesDiscarded = new HashSet<String>();
		if (peptides.size() > 1) {
			for (int i = 0; i < peptides.size(); i++) {
				QuantifiedPeptideInterface peptide1 = peptides.get(i);
				final String sequence1 = peptide1.getSequence();
				// discard if it doesn't contain any quantified aa
				if (!PCQUtils.containsAny(sequence1, quantifiedAAs)) {
					if (!peptideSequencesDiscarded.contains(sequence1)) {
						peptideSequencesDiscarded.add(sequence1);
						log.info(sequence1 + " discarded for not having one quantitation sites from '"
								+ StringUtils.getSeparatedValueStringFromChars(quantifiedAAs, ",") + "'");
					}
					individualQuantifiedPeptideSet.remove(peptide1);
					continue;
				}
				// if it is not isobaric peptide and it has more than one
				// quantified site, it cannot be distinguish, so we discard it.
				if (!(peptide1 instanceof IsobaricQuantifiedPeptide)
						&& PCQUtils.howManyContains(sequence1, quantifiedAAs) > 1) {
					if (!peptideSequencesDiscarded.contains(sequence1)) {
						peptideSequencesDiscarded.add(sequence1);
						log.info(sequence1 + " discarded for having ambiguous quantitation sites from "
								+ StringUtils.getSeparatedValueStringFromChars(quantifiedAAs, ",") + "'");
					}
					individualQuantifiedPeptideSet.remove(peptide1);
					continue;
				}

				// get the keys from the peptide.
				// not that the peptide could have more than one key because 2
				// reasons:
				// - it could be shared by more than one protein
				// - it could have more than one quantified aminoacid in its
				// sequence
				Map<PositionInPeptide, List<PositionInProtein>> proteinKeysByPeptide1Keys = peptide1
						.getProteinKeysByPeptideKeysForQuantifiedAAs(quantifiedAAs, uplr);
				if (proteinKeysByPeptide1Keys.isEmpty()) {
					// peptides without a site mapped to a protein are discarded
					if (!peptideSequencesDiscarded.contains(sequence1)) {
						peptideSequencesDiscarded.add(sequence1);
						log.warn("Peptide '" + sequence1 + "' cannot be mapped to any protein sequence");
					}
					individualQuantifiedPeptideSet.remove(peptide1);
					continue;
				}
				for (int j = i + 1; j < peptides.size(); j++) {
					QuantifiedPeptideInterface peptide2 = peptides.get(j);
					final String sequence2 = peptide2.getSequence();
					// discard if it doesn't contain any quantified aa
					if (!PCQUtils.containsAny(sequence2, quantifiedAAs)) {
						if (!peptideSequencesDiscarded.contains(sequence2)) {
							peptideSequencesDiscarded.add(sequence2);
							log.info(sequence2 + " discarded for not having one quantitation sites from '"
									+ StringUtils.getSeparatedValueStringFromChars(quantifiedAAs, ",") + "'");
						}
						individualQuantifiedPeptideSet.remove(peptide2);
						continue;
					}
					// if it is not isobaric peptide and it has more than one
					// quantified site, it cannot be distinguish, so we discard
					// it.
					if (!(peptide2 instanceof IsobaricQuantifiedPeptide)
							&& PCQUtils.howManyContains(sequence2, quantifiedAAs) > 1) {
						if (!peptideSequencesDiscarded.contains(sequence2)) {
							peptideSequencesDiscarded.add(sequence2);
							log.info(sequence2 + " discarded for having ambiguous quantitation sites from '"
									+ StringUtils.getSeparatedValueStringFromChars(quantifiedAAs, ",") + "'");
						}
						individualQuantifiedPeptideSet.remove(peptide2);
						continue;
					}
					// get the keys from the peptide.
					// not that the peptide could have more than one key because
					// 2
					// reasons:
					// - it could be shared by more than one protein
					// - it could have more than one quantified aminoacid in its
					// sequence
					Map<PositionInPeptide, List<PositionInProtein>> proteinKeysByPeptide2Keys = peptide2
							.getProteinKeysByPeptideKeysForQuantifiedAAs(quantifiedAAs, uplr);
					if (proteinKeysByPeptide2Keys.isEmpty()) {
						if (!peptideSequencesDiscarded.contains(sequence2)) {
							peptideSequencesDiscarded.add(sequence2);
							// peptides without a site mapped to a protein are
							// discarded
							log.warn("Peptide '" + sequence2 + "' cannot be mapped to any protein sequence");
						}
						individualQuantifiedPeptideSet.remove(peptide2);
						continue;
					}

					// iterate over all the keys
					for (PositionInPeptide positionInPeptide1 : proteinKeysByPeptide1Keys.keySet()) {
						List<PositionInProtein> proteinKeysFromPeptide1 = proteinKeysByPeptide1Keys
								.get(positionInPeptide1);
						for (PositionInPeptide positionInPeptide2 : proteinKeysByPeptide2Keys.keySet()) {
							List<PositionInProtein> proteinKeysFromPeptide2 = proteinKeysByPeptide2Keys
									.get(positionInPeptide2);

							// now I compare the two list of keys
							// if they are equal, that means, if they share the
							// same sites of the same proteins
							if (PCQUtils.areEquals(proteinKeysFromPeptide1, proteinKeysFromPeptide2)) {

								String key = PCQUtils.getPositionsInProteinsKey(proteinKeysFromPeptide1);

								// collapsed in the same peptide node
								PCQPeptideNode peptideNode = null;
								// if the protein was already associated with
								// some protein node
								if (peptideNodesByPeptideNodeKey.containsKey(key)) {
									peptideNode = peptideNodesByPeptideNodeKey.get(key);
								} else {
									peptideNode = new PCQPeptideNode(this, key,
											new Pair<QuantifiedPeptideInterface, PositionInPeptide>(peptide1,
													positionInPeptide1),
											new Pair<QuantifiedPeptideInterface, PositionInPeptide>(peptide2,
													positionInPeptide2));
								}
								// add the two peptides to the peptide node
								peptideNode.addQuantifiedPeptide(peptide1, positionInPeptide1);
								peptideNode.addQuantifiedPeptide(peptide2, positionInPeptide2);
								// add to the set of nodes
								peptideNodes.add(peptideNode);
								// add to the map
								peptideNodesByPeptideNodeKey.put(key, peptideNode);

							} else if (PCQUtils.shareAny(proteinKeysFromPeptide1, proteinKeysFromPeptide2)) {
								// they have keys in common, but not all of them
								// now, create a node with all this keys unique
								// to the first peptide
								List<PositionInProtein> uniqueKeysForPeptide1 = PCQUtils
										.getUniqueToFirst(proteinKeysFromPeptide1, proteinKeysFromPeptide2);
								if (!uniqueKeysForPeptide1.isEmpty()) {

									String key1 = PCQUtils.getPositionsInProteinsKey(uniqueKeysForPeptide1);
									PCQPeptideNode peptideNode = null;
									if (peptideNodesByPeptideNodeKey.containsKey(key1)) {
										peptideNode = peptideNodesByPeptideNodeKey.get(key1);
									} else {
										peptideNode = new PCQPeptideNode(this, key1,
												new Pair<QuantifiedPeptideInterface, PositionInPeptide>(peptide1,
														positionInPeptide1));
									}
									// add the two peptides to the peptide node
									peptideNode.addQuantifiedPeptide(peptide1, positionInPeptide1);
									// add to the set of nodes
									peptideNodes.add(peptideNode);
									// add to the map
									peptideNodesByPeptideNodeKey.put(key1, peptideNode);
								}
								// now, create a node with all this keys unique
								// to the second peptide
								List<PositionInProtein> uniqueKeysForPeptide2 = PCQUtils
										.getUniqueToFirst(proteinKeysFromPeptide2, proteinKeysFromPeptide1);
								if (!uniqueKeysForPeptide2.isEmpty()) {
									String key2 = PCQUtils.getPositionsInProteinsKey(uniqueKeysForPeptide2);
									PCQPeptideNode peptideNode = null;
									if (peptideNodesByPeptideNodeKey.containsKey(key2)) {
										peptideNode = peptideNodesByPeptideNodeKey.get(key2);
									} else {
										peptideNode = new PCQPeptideNode(this, key2,
												new Pair<QuantifiedPeptideInterface, PositionInPeptide>(peptide2,
														positionInPeptide2));
									}
									// add the two peptides to the peptide node
									peptideNode.addQuantifiedPeptide(peptide2, positionInPeptide2);
									// add to the set of nodes
									peptideNodes.add(peptideNode);
									// add to the map
									peptideNodesByPeptideNodeKey.put(key2, peptideNode);
								}
								// now, create a node with all this keys in
								// common between peptide1 and peptide2
								List<PositionInProtein> sharedKeys = PCQUtils.getInCommon(proteinKeysFromPeptide2,
										proteinKeysFromPeptide1);
								if (!sharedKeys.isEmpty()) {
									String sharedKey = PCQUtils.getPositionsInProteinsKey(sharedKeys);
									PCQPeptideNode peptideNode = null;
									if (peptideNodesByPeptideNodeKey.containsKey(sharedKey)) {
										peptideNode = peptideNodesByPeptideNodeKey.get(sharedKey);
									} else {
										peptideNode = new PCQPeptideNode(this, sharedKey,
												new Pair<QuantifiedPeptideInterface, PositionInPeptide>(peptide1,
														positionInPeptide1),
												new Pair<QuantifiedPeptideInterface, PositionInPeptide>(peptide2,
														positionInPeptide2));
									}
									// add the two peptides to the peptide node
									peptideNode.addQuantifiedPeptide(peptide1, positionInPeptide1);
									peptideNode.addQuantifiedPeptide(peptide2, positionInPeptide2);
									// add to the set of nodes
									peptideNodes.add(peptideNode);
									// add to the map
									peptideNodesByPeptideNodeKey.put(sharedKey, peptideNode);
								}
							} else {
								// they dont share any key create a peptide node
								// for each peptide separated peptide node for
								// psm1
								String key1 = PCQUtils.getPositionsInProteinsKey(proteinKeysFromPeptide1);
								PCQPeptideNode peptideNode = null;
								if (peptideNodesByPeptideNodeKey.containsKey(key1)) {
									peptideNode = peptideNodesByPeptideNodeKey.get(key1);
								} else {
									peptideNode = new PCQPeptideNode(this, key1,
											new Pair<QuantifiedPeptideInterface, PositionInPeptide>(peptide1,
													positionInPeptide1));
								}
								// add the two peptides to the peptide node
								peptideNode.addQuantifiedPeptide(peptide1, positionInPeptide1);
								peptideNodes.add(peptideNode);
								peptideNodesByPeptideNodeKey.put(key1, peptideNode);

								// protein node for protein2
								String key2 = PCQUtils.getPositionsInProteinsKey(proteinKeysFromPeptide2);
								PCQPeptideNode peptideNode2 = null;
								if (peptideNodesByPeptideNodeKey.containsKey(key2)) {
									peptideNode2 = peptideNodesByPeptideNodeKey.get(key2);
								} else {
									peptideNode2 = new PCQPeptideNode(this, key2,
											new Pair<QuantifiedPeptideInterface, PositionInPeptide>(peptide2,
													positionInPeptide2));
								}
								peptideNode2.addQuantifiedPeptide(peptide2, positionInPeptide2);
								peptideNodes.add(peptideNode2);
								peptideNodesByPeptideNodeKey.put(key2, peptideNode2);
							}
						}
					}
				}
			}
		} else {
			// only one peptide
			// create a peptide node for each position in the peptide
			QuantifiedPeptideInterface peptide = peptides.iterator().next();
			Map<PositionInPeptide, List<PositionInProtein>> proteinKeysByPeptideKeys = peptide
					.getProteinKeysByPeptideKeysForQuantifiedAAs(quantifiedAAs, uplr);

			for (PositionInPeptide positionInPeptide : proteinKeysByPeptideKeys.keySet()) {
				List<PositionInProtein> positionsInProtein = proteinKeysByPeptideKeys.get(positionInPeptide);

				String key = PCQUtils.getPositionsInProteinsKey(positionsInProtein);
				PCQPeptideNode peptideNode = null;
				if (peptideNodesByPeptideNodeKey.containsKey(key)) {
					peptideNode = peptideNodesByPeptideNodeKey.get(key);
				} else {
					peptideNode = new PCQPeptideNode(this, key,
							new Pair<QuantifiedPeptideInterface, PositionInPeptide>(peptide, positionInPeptide));
				}
				peptideNode.addQuantifiedPeptide(peptide, positionInPeptide);
				peptideNodes.add(peptideNode);
				peptideNodesByPeptideNodeKey.put(key, peptideNode);
			}
		}

	}

	private void createProteinNodes() {
		// create a map to store proteins by accession
		Map<String, Set<QuantifiedProteinInterface>> proteinMap = new THashMap<String, Set<QuantifiedProteinInterface>>();
		for (QuantifiedProteinInterface protein : individualQuantifiedProteinSet) {
			if (proteinMap.containsKey(protein.getKey())) {
				proteinMap.get(protein.getKey()).add(protein);
			} else {
				Set<QuantifiedProteinInterface> set = new THashSet<QuantifiedProteinInterface>();
				set.add(protein);
				proteinMap.put(protein.getKey(), set);
			}
		}
		List<String> keyList = new ArrayList<String>();
		keyList.addAll(proteinMap.keySet());
		if (proteinMap.size() > 1) {
			for (int i = 0; i < keyList.size(); i++) {
				final String key1 = keyList.get(i);
				Set<QuantifiedProteinInterface> proteins1 = proteinMap.get(key1);
				for (int j = i + 1; j < keyList.size(); j++) {
					final String key2 = keyList.get(j);
					Set<QuantifiedProteinInterface> proteins2 = proteinMap.get(key2);

					if (getParams().isCollapseIndistinguishableProteins()
							&& PCQUtils.proteinsShareAllPeptides(proteins1, proteins2)) {
						PCQProteinNode proteinNode = null;
						PCQProteinNode proteinNode2 = null;
						// if the protein was already associated with some
						// protein node
						if (proteinNodesByProteinKey.containsKey(key1)) {
							proteinNode = proteinNodesByProteinKey.get(key1);
						}
						if (proteinNodesByProteinKey.containsKey(key2)) {
							proteinNode2 = proteinNodesByProteinKey.get(key2);
						}
						if (proteinNode == null && proteinNode2 == null) {
							// if non of the proteins are in any protein node
							// yet, create a new one
							proteinNode = new PCQProteinNode(this, proteins1, proteins2);
						} else {
							// some of the proteins are assigned to some protein
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
									for (String key : proteinNodesByProteinKey.keySet()) {
										if (proteinNodesByProteinKey.get(key).equals(proteinNode2)) {
											proteinNodesByProteinKey.put(key, proteinNode);
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
						proteinNodesByProteinKey.put(key1, proteinNode);
						proteinNodesByProteinKey.put(key2, proteinNode);

					} else {
						// create a protein node for each protein separated
						// protein node for protein1
						PCQProteinNode proteinNode = null;
						if (proteinNodesByProteinKey.containsKey(key1)) {
							proteinNode = proteinNodesByProteinKey.get(key1);
						} else {
							proteinNode = new PCQProteinNode(this, proteins1);
						}
						proteinNodes.add(proteinNode);
						proteinNodesByProteinKey.put(key1, proteinNode);

						// protein node for protein2
						PCQProteinNode proteinNode2 = null;
						if (proteinNodesByProteinKey.containsKey(key2)) {
							proteinNode2 = proteinNodesByProteinKey.get(key2);
						} else {
							proteinNode2 = new PCQProteinNode(this, proteins2);
						}
						proteinNodes.add(proteinNode2);
						proteinNodesByProteinKey.put(key2, proteinNode2);

					}
				}
			}
		} else {
			// only one protein
			// create a protein node for the protein
			// protein node for protein1
			if (proteinMap.isEmpty()) {
				log.info(this);
			}
			final Collection<QuantifiedProteinInterface> proteins = proteinMap.values().iterator().next();
			PCQProteinNode proteinNode = new PCQProteinNode(this, proteins);
			proteinNodes.add(proteinNode);
			proteinNodesByProteinKey.put(proteinMap.keySet().iterator().next(), proteinNode);

		}
		// log.debug(proteinNodes.size() + " protein nodes created in cluster");
		// for (PCQProteinNode proteinNode : proteinNodes) {
		// log.info(proteinNode);
		// }

	}

	/**
	 * removes any node stored in the object
	 */
	private void resetNodes() {
		resetPeptideNodes();
		resetProteinNodes();
	}

	private void resetProteinNodes() {
		proteinNodes.clear();
		proteinNodesByProteinKey.clear();
	}

	private void resetPeptideNodes() {
		peptideNodes.clear();
		peptideNodesByPeptideNodeKey.clear();
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
				NWResult result = NWAlign.needlemanWunsch(pep1.getSequence(), pep2.getSequence());
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

	public void createPairs(AlignmentSet peptideAlignments) {
		proteinPairs.clear();
		discardedProteinPairs.clear();
		// if
		// (ProteinClusterQuantParameters.getInstance().isCollapseIndistinguishableProteins())
		// {
		// createPairsCollapsingIndistinguisibleProteins(gAM);
		// } else {
		// createPairsWithoutCollapsing(gAM);
		// }
		createProteinPairs(peptideAlignments);
	}

	private void createProteinPairs(AlignmentSet peptideAlignments) {

		List<PCQProteinNode> proteinNodeList = new ArrayList<PCQProteinNode>();
		proteinNodeList.addAll(proteinNodes);

		for (int i = 0; i < proteinNodeList.size(); i++) {
			for (int j = i + 1; j < proteinNodeList.size(); j++) {
				PCQProteinNode proteinNode1 = proteinNodeList.get(i);

				PCQProteinNode proteinNode2 = proteinNodeList.get(j);
				if (PCQUtils.shareAtLeastOnePeptideNode(proteinNode1, proteinNode2, false)
						|| (peptideAlignments != null && PCQUtils.shareAtLeastOnePeptideBySimilarity(proteinNode1,
								proteinNode2, peptideAlignments, false))) {
					ProteinPair pair = new ProteinPair(proteinNode1, proteinNode2);
					if (pair.isContainsDiscardedProteinNode()) {
						discardedProteinPairs.add(pair);
					} else {
						proteinPairs.add(pair);
					}
				}
			}
		}
	}

	public Set<ProteinPair> getProteinPairs() {
		return proteinPairs;
	}

	/**
	 * Get {@link ProteinPair}s containing at least a discarded
	 * {@link PCQProteinNode}
	 *
	 * @return
	 */
	public Set<ProteinPair> getDiscardedProteinPairs() {
		return discardedProteinPairs;
	}

	/**
	 * Get all protein pairs, the regular ones and the discarded ones
	 * (containing a discarded protein node)
	 *
	 * @return
	 */
	public Set<ProteinPair> getAllProteinPairs() {
		Set<ProteinPair> ret = new THashSet<ProteinPair>();
		ret.addAll(getProteinPairs());
		ret.addAll(getDiscardedProteinPairs());
		return ret;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("PROTN[");
		for (PCQProteinNode protein : PCQUtils.getSortedProteinNodesByAcc(proteinNodes)) {
			if (!"PROTN[".equals(sb.toString())) {
				sb.append(",");
			}
			sb.append(protein.getKey());
		}
		sb.append("]NTORP");
		StringBuilder sb2 = new StringBuilder();
		sb2.append("PEPN[");
		for (PCQPeptideNode peptide : PCQUtils.getSortedPeptideNodesBySequence(peptideNodes)) {
			if (!"PEPN[".equals(sb2.toString())) {
				sb2.append(",");
			}
			sb2.append(peptide.getKey());
		}
		sb2.append("]NPEP");
		final String string = sb.append("\t").append(sb2).toString();
		StringBuilder sb3 = new StringBuilder();
		List<QuantifiedProteinInterface> proteinList = new ArrayList<QuantifiedProteinInterface>();
		proteinList.addAll(getProteinSet());
		Collections.sort(proteinList, new Comparator<QuantifiedProteinInterface>() {
			@Override
			public int compare(QuantifiedProteinInterface o1, QuantifiedProteinInterface o2) {
				return o1.getKey().compareTo(o2.getKey());
			}
		});
		for (QuantifiedProteinInterface protein : proteinList) {
			if (!"".equals(sb3.toString())) {
				sb.append(",");
			}
			sb3.append(protein + "\n");
		}
		return string + "\n" + sb3.toString();
	}

	public void addAlignment(AlignedPeptides alignment) {
		alignmentSet.addAlignment(alignment);
	}

	/**
	 * If align was performed, this will return the corresponding peptide nodes
	 * aligned to a given one in the cluster.
	 *
	 * @param peptide
	 * @return
	 */
	public Set<PCQPeptideNode> getAlignedPeptideNodes(PCQPeptideNode peptideNode) {

		if (peptideNode.getQuantifiedPeptides().size() > 1) {
			// If peptide nodes contain more than one different peptide,
			// then the aligments cannot be done, so there is something
			// wrong
			return Collections.emptySet();
		}
		QuantifiedPeptideInterface peptide = peptideNode.getQuantifiedPeptides().iterator().next();

		Set<PCQPeptideNode> ret = new THashSet<PCQPeptideNode>();

		for (AlignedPeptides alignment : alignmentSet.getAlignmentsForPeptide(peptide)) {
			ret.add(peptideNodesByPeptideNodeKey.get(alignment.getPeptide1().getKey()));
			ret.add(peptideNodesByPeptideNodeKey.get(alignment.getPeptide2().getKey()));
		}
		// not include the same as peptideNode
		ret.remove(peptideNode);

		return ret;
	}

	/**
	 * Get the {@link NWResult} from two specific peptides if available.
	 *
	 * @param peptideNode
	 * @param peptideNode2
	 * @return
	 */
	public NWResult getAlignmentResult(PCQPeptideNode peptideNode, PCQPeptideNode peptideNode2) {
		final QuantifiedPeptideInterface peptide1 = peptideNode.getQuantifiedPeptides().iterator().next();
		final QuantifiedPeptideInterface peptide2 = peptideNode2.getQuantifiedPeptides().iterator().next();

		final Set<AlignedPeptides> alignments = alignmentSet.getAlignmentsForPeptide(peptide1);
		for (AlignedPeptides alignment : alignments) {
			if (alignment.getPeptideAligned(peptide1).equals(peptide2)) {
				return alignment.getAlignmentResult();
			}
		}

		return null;
	}

	public Set<String> getPeptideNodeKeys() {
		Set<String> ret = new THashSet<String>();
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
		Set<String> set = new THashSet<String>();
		for (PCQPeptideNode peptideNode : getNonDiscardedPeptideNodes()) {
			for (QuantifiedPeptideInterface peptide : peptideNode.getItemsInNode()) {
				set.add(peptide.getSequence());
			}
		}
		return set.size();
	}

	/**
	 * Get the number of different individual peptides in the cluster
	 *
	 * @return
	 */
	public int getNumDifferentDiscardedIndividualPeptides() {
		Set<String> set = new THashSet<String>();
		for (PCQPeptideNode peptideNode : getPeptideNodes()) {
			if (peptideNode.isDiscarded()) {
				for (QuantifiedPeptideInterface peptide : peptideNode.getItemsInNode()) {
					set.add(peptide.getSequence());
				}
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
		Set<String> set = new THashSet<String>();
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
		Set<PCQProteinNode> ret = new THashSet<PCQProteinNode>();
		for (PCQProteinNode proteinNode : getProteinNodes()) {
			if (!proteinNode.isDiscarded()) {
				ret.add(proteinNode);
			}
		}
		return ret;
	}

	public Set<PCQPeptideNode> getNonDiscardedPeptideNodes() {
		Set<PCQPeptideNode> ret = new THashSet<PCQPeptideNode>();
		for (PCQPeptideNode peptideNode : getPeptideNodes()) {
			if (!peptideNode.isDiscarded()) {
				ret.add(peptideNode);
			} else {
				// log.info("asdf");
			}
		}
		return ret;
	}

	public Set<QuantifiedProteinInterface> getNonDiscardedProteinSet() {
		Set<QuantifiedProteinInterface> ret = new THashSet<QuantifiedProteinInterface>();
		final Set<QuantifiedProteinInterface> proteinSet = getProteinSet();
		for (QuantifiedProteinInterface protein : proteinSet) {
			if (!protein.isDiscarded()) {
				ret.add(protein);
			}
		}
		return ret;
	}

	public void removeIndividualProteinsWithNoPSMs() {
		final Iterator<PCQProteinNode> proteinNodeIterator = proteinNodes.iterator();
		while (proteinNodeIterator.hasNext()) {
			final PCQProteinNode proteinNode = proteinNodeIterator.next();
			final Iterator<QuantifiedProteinInterface> proteinIterator = proteinNode.getItemsInNode().iterator();
			while (proteinIterator.hasNext()) {
				final QuantifiedProteinInterface protein = proteinIterator.next();
				if (protein.getQuantifiedPSMs().isEmpty()) {
					proteinIterator.remove();
				}
			}
			if (proteinNode.getQuantifiedProteins().isEmpty()) {
				proteinNodeIterator.remove();
				PCQFilter.getDiscardedProteinNodes().add(proteinNode);
			}
		}
		// check if some psms have proteins with no psms linked, and if so,
		// unlink them from the psms
		final Iterator<PCQPeptideNode> peptideNodeIterator = peptideNodes.iterator();
		while (peptideNodeIterator.hasNext()) {
			final PCQPeptideNode peptideNode = peptideNodeIterator.next();
			final Iterator<QuantifiedPeptideInterface> peptideIterator = peptideNode.getItemsInNode().iterator();
			while (peptideIterator.hasNext()) {
				final QuantifiedPeptideInterface peptide = peptideIterator.next();
				final Set<QuantifiedPSMInterface> quantifiedPSMs = peptide.getQuantifiedPSMs();
				if (quantifiedPSMs.isEmpty()) {
					peptideIterator.remove();
				}
				for (QuantifiedPSMInterface psm : quantifiedPSMs) {
					final Iterator<QuantifiedProteinInterface> proteinIterator = psm.getQuantifiedProteins().iterator();
					while (proteinIterator.hasNext()) {
						QuantifiedProteinInterface protein = proteinIterator.next();

						if (protein.getQuantifiedPSMs().isEmpty()) {
							// unlink psm to that protein
							proteinIterator.remove();
						}
					}

				}
			}
			if (peptideNode.getQuantifiedPSMs().isEmpty()) {
				peptideNodeIterator.remove();
				PCQFilter.getDiscardedPeptideNodes().add(peptideNode);
			}
		}
	}

	public boolean removeÍndividualProtein(QuantifiedProteinInterface protein) {
		if (protein.getKey().equals("A0A075B7C1")) {
			log.info(protein);
		}
		proteinNodesByProteinKey.remove(protein.getKey());
		return getProteinSet().remove(protein);
	}

	public Map<String, PCQProteinNode> getProteinNodesByProteinKey() {
		return proteinNodesByProteinKey;

	}

	/**
	 * @return the peptideNodesByPeptideSequence
	 */
	public Map<String, PCQPeptideNode> getPeptideNodesByPeptideSequence() {
		return peptideNodesByPeptideNodeKey;
	}

}
