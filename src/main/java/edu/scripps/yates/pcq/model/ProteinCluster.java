package edu.scripps.yates.pcq.model;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.census.read.model.IsobaricQuantifiedPeptide;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.pcq.filter.PCQFilter;
import edu.scripps.yates.pcq.params.ProteinClusterQuantParameters;
import edu.scripps.yates.pcq.util.DiscardedPeptidesSet;
import edu.scripps.yates.pcq.util.DiscardedPeptidesSet.DISCARD_REASON;
import edu.scripps.yates.pcq.util.PCQUtils;
import edu.scripps.yates.pcq.xgmml.util.AlignedPeptides;
import edu.scripps.yates.pcq.xgmml.util.AlignmentSet;
import edu.scripps.yates.pcq.xgmml.util.ProteinNodeLabel;
import edu.scripps.yates.utilities.alignment.nwalign.NWAlign;
import edu.scripps.yates.utilities.alignment.nwalign.NWResult;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.sequence.PTMInProtein;
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
	 * By calling this function, the cluster will reorganize by collapsing the nodes
	 * properly and calculating the ratios accordingly
	 * 
	 * @param annotatedProteins
	 * @param proteoFormFastaReader
	 * 
	 * @throws IOException
	 */
	public void createNodes(Map<String, Entry> annotatedProteins) throws IOException {
		// reset nodes
		// removes all nodes in the cluster and clear all maps
		resetNodes();

		// create nodes with the individual proteins and peptides

		// first create protein nodes
		createProteinNodes();

		// create nodes peptide nodes
		if (getParams().isCollapseBySites()) {
			createPeptideNodesByQuantifiedSites(getParams().getAaQuantified());
		} else if (getParams().isCollapseByPTMs()) {
			createPeptideNodesByPTMs();
		} else {
			createPeptideNodes();
		}

		// connect protein and peptide nodes
		connectProteinAndPeptideNodes();

		// set description and taxonomy if they dont have it
		for (final PCQProteinNode proteinNode : getProteinNodes()) {
			proteinNode.annotateProteins(annotatedProteins);
		}

	}

	private void connectProteinAndPeptideNodes() {

		final Set<PCQPeptideNode> peptideNodes2 = getPeptideNodes();
		for (final PCQPeptideNode peptideNode : peptideNodes2) {

			final Set<QuantifiedPeptideInterface> quantifiedPeptides = peptideNode.getQuantifiedPeptides();
			for (final QuantifiedPeptideInterface peptide : quantifiedPeptides) {
				for (final QuantifiedProteinInterface protein : peptide.getQuantifiedProteins()) {
					final PCQProteinNode proteinNode = proteinNodesByProteinKey.get(protein.getKey());
					if (proteinNode == null) {
						log.info(this);
						boolean found = false;
						for (final QuantifiedProteinInterface protein2 : getProteinSet()) {
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
		final List<QuantifiedPeptideInterface> peptides = new ArrayList<QuantifiedPeptideInterface>();
		peptides.addAll(individualQuantifiedPeptideSet);
		final Set<String> set = new THashSet<String>();
		for (final QuantifiedPeptideInterface pep : peptides) {
			if (pep.getSequence().equals("FEELCSDLFR")) {
				log.info(pep);
			}
			if (set.contains(pep.getKey())) {
				log.info("Inconsistency error. 2 peptides with the same key cannot exist!");
			} else {
				set.add(pep.getKey());
			}
		}
		if (peptides.size() > 1) {
			for (int i = 0; i < peptides.size(); i++) {
				final QuantifiedPeptideInterface peptide1 = peptides.get(i);
				if (peptide1.getSequence().equals("FEELCSDLFR")) {
					log.info(peptide1);
				}
				for (int j = i + 1; j < peptides.size(); j++) {
					final QuantifiedPeptideInterface peptide2 = peptides.get(j);
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
			final QuantifiedPeptideInterface peptide = peptides.iterator().next();
			final PCQPeptideNode peptideNode = new PCQPeptideNode(this, peptide);
			peptideNodes.add(peptideNode);
			peptideNodesByPeptideNodeKey.put(peptide.getKey(), peptideNode);
		}

	}

	/**
	 * Iterates all peptides in the cluster over all others and creates a node for
	 * all peptides that have the same key, which, in this case is the position of
	 * the quantifiedAA in the protein.
	 * 
	 * @param quantifiedAAs
	 */
	private void createPeptideNodesByQuantifiedSites(char[] quantifiedAAs) {
		final StringBuilder quantifiedAAsString = new StringBuilder();
		for (final char c : quantifiedAAs) {
			quantifiedAAsString.append(c);
		}
		final UniprotProteinLocalRetriever uplr = PCQUtils
				.getUniprotProteinLocalRetrieverByFolder(getParams().getUniprotReleasesFolder());
		final boolean useProteinGeneName = getParams().getProteinLabel() == ProteinNodeLabel.GENE;
		final boolean useProteinID = getParams().getProteinLabel() == ProteinNodeLabel.ID;

		final List<QuantifiedPeptideInterface> peptides = new ArrayList<QuantifiedPeptideInterface>();
		peptides.addAll(individualQuantifiedPeptideSet);
		final Set<String> set = new THashSet<String>();
		for (final QuantifiedPeptideInterface pep : peptides) {
			if (set.contains(pep.getKey())) {
				log.warn("Inconsistency error. 2 peptides with the same key cannot exist!");
			} else {
				set.add(pep.getKey());
			}
		}
		// to not include the same peptide in different peptide nodes
		final Map<QuantifiedPeptideInterface, PCQPeptideNode> peptideNodesByPeptides = new THashMap<QuantifiedPeptideInterface, PCQPeptideNode>();
		final Set<String> peptideSequencesDiscarded = new THashSet<String>();
		if (peptides.size() > 1) {
			for (int i = 0; i < peptides.size(); i++) {
				final QuantifiedPeptideInterface peptide1 = peptides.get(i);
				final String sequence1 = peptide1.getSequence();
				// discard if it doesn't contain any quantified aa
				if (getParams().isRemoveFilteredNodes() && !PCQUtils.containsAny(sequence1, quantifiedAAs)) {
					if (!peptideSequencesDiscarded.contains(sequence1)) {
						DiscardedPeptidesSet.getInstance().add(peptide1, DISCARD_REASON.PEPTIDE_WITH_NO_QUANT_SITE,
								"quant site(s): " + quantifiedAAsString);
						peptideSequencesDiscarded.add(sequence1);
						log.warn(sequence1 + " discarded for not having one quantitation sites from '"
								+ StringUtils.getSeparatedValueStringFromChars(quantifiedAAs, ",") + "'");
					}
					individualQuantifiedPeptideSet.remove(peptide1);
					continue;
				}
				// if it is not isobaric peptide and it has more than one
				// quantified site, it cannot be distinguish, so we discard it.
				if (getParams().isRemoveFilteredNodes() && !(peptide1 instanceof IsobaricQuantifiedPeptide)
						&& PCQUtils.howManyContains(sequence1, quantifiedAAs) > 1) {
					if (!peptideSequencesDiscarded.contains(sequence1)) {
						peptideSequencesDiscarded.add(sequence1);
						DiscardedPeptidesSet.getInstance().add(peptide1, DISCARD_REASON.AMBIGOUS_QUANT,
								"quant site(s): " + quantifiedAAsString);
						log.warn(sequence1 + " discarded for having ambiguous quantitation sites from "
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
				final Map<PositionInPeptide, List<PositionInProtein>> proteinKeysByPeptide1Keys = peptide1
						.getProteinKeysByPeptideKeysForQuantifiedAAs(quantifiedAAs, uplr, PCQUtils.proteinSequences);
				if (proteinKeysByPeptide1Keys.isEmpty()) {

					if (getParams().isRemoveFilteredNodes() || getParams().isCollapseBySites()) {
						// peptides without a site mapped to a protein are discarded
						if (!peptideSequencesDiscarded.contains(sequence1)) {
							peptideSequencesDiscarded.add(sequence1);
							DiscardedPeptidesSet.getInstance().add(peptide1, DISCARD_REASON.PEPTIDE_WITH_NO_PROTEIN);
							log.warn("Peptide '" + sequence1
									+ "' cannot be mapped to any protein sequence and it will be ignored");
						}
						individualQuantifiedPeptideSet.remove(peptide1);
						continue;
					} else {
						log.warn("Peptide '" + peptide1.getSequence() + "' cannot be mapped to any protein sequence");
					}
				}
				for (int j = i + 1; j < peptides.size(); j++) {
					final QuantifiedPeptideInterface peptide2 = peptides.get(j);
					final String sequence2 = peptide2.getSequence();

					// discard if it doesn't contain any quantified aa
					if (getParams().isRemoveFilteredNodes() && !PCQUtils.containsAny(sequence2, quantifiedAAs)) {
						if (!peptideSequencesDiscarded.contains(sequence2)) {
							DiscardedPeptidesSet.getInstance().add(peptide2, DISCARD_REASON.PEPTIDE_WITH_NO_QUANT_SITE,
									"quant site(s): " + quantifiedAAsString);

							peptideSequencesDiscarded.add(sequence2);
							log.warn(sequence2 + " discarded for not having one quantitation sites from '"
									+ StringUtils.getSeparatedValueStringFromChars(quantifiedAAs, ",") + "'");
						}
						individualQuantifiedPeptideSet.remove(peptide2);
						continue;
					}
					// if it is not isobaric peptide and it has more than one
					// quantified site, it cannot be distinguish, so we discard
					// it.
					if (getParams().isRemoveFilteredNodes() && !(peptide2 instanceof IsobaricQuantifiedPeptide)
							&& PCQUtils.howManyContains(sequence2, quantifiedAAs) > 1) {
						if (!peptideSequencesDiscarded.contains(sequence2)) {
							DiscardedPeptidesSet.getInstance().add(peptide2, DISCARD_REASON.AMBIGOUS_QUANT,
									"quant site(s): " + quantifiedAAsString);
							peptideSequencesDiscarded.add(sequence2);
							log.warn(sequence2 + " discarded for having ambiguous quantitation sites from '"
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

					final Map<PositionInPeptide, List<PositionInProtein>> proteinKeysByPeptide2Keys = peptide2
							.getProteinKeysByPeptideKeysForQuantifiedAAs(quantifiedAAs, uplr,
									PCQUtils.proteinSequences);
					if (proteinKeysByPeptide2Keys.isEmpty()) {
						if (getParams().isRemoveFilteredNodes() || getParams().isCollapseBySites()) {
							if (!peptideSequencesDiscarded.contains(sequence2)) {
								DiscardedPeptidesSet.getInstance().add(peptide2,
										DISCARD_REASON.PEPTIDE_WITH_NO_PROTEIN);
								peptideSequencesDiscarded.add(sequence2);
								// peptides without a site mapped to a protein are
								// discarded
								log.warn("Peptide '" + sequence2
										+ "' cannot be mapped to any protein sequence and it will be ignored");
							}
							individualQuantifiedPeptideSet.remove(peptide2);
							continue;
						} else {
							log.warn("Peptide '" + peptide2.getSequence()
									+ "' cannot be mapped to any protein sequence");

						}
					}
					// if we have isobaric isotopologues, we could have several
					// ratios in a single peptide with multiple quantified
					// sites<br>Otherwise, we have to consider the
					// positionsInPeptides all together per peptide
					if (peptide1 instanceof IsobaricQuantifiedPeptide
							&& peptide2 instanceof IsobaricQuantifiedPeptide) {
						// iterate over all the keys
						for (final PositionInPeptide positionInPeptide1 : proteinKeysByPeptide1Keys.keySet()) {
							final List<PositionInProtein> proteinKeysFromPeptide1 = proteinKeysByPeptide1Keys
									.get(positionInPeptide1);
							for (final PositionInPeptide positionInPeptide2 : proteinKeysByPeptide2Keys.keySet()) {
								final List<PositionInProtein> proteinKeysFromPeptide2 = proteinKeysByPeptide2Keys
										.get(positionInPeptide2);

								final String key1 = QuantUtils.getPositionsInProteinsKey(proteinKeysFromPeptide1,
										useProteinGeneName, useProteinID, uplr, getParams().getUniprotVersion());
								final String key2 = QuantUtils.getPositionsInProteinsKey(proteinKeysFromPeptide2,
										useProteinGeneName, useProteinID, uplr, getParams().getUniprotVersion());
								// now I compare the two list of keys
								// if they are equal, that means, if they share
								// the
								// same sites of the same proteins
								if (PCQUtils.areEquals(proteinKeysFromPeptide1, proteinKeysFromPeptide2)) {

									// key1 and key2 should be the same

									// collapsed in the same peptide node
									PCQPeptideNode peptideNode = null;
									// if the protein was already associated
									// with
									// some protein node
									if (peptideNodesByPeptideNodeKey.containsKey(key1)) {
										peptideNode = peptideNodesByPeptideNodeKey.get(key1);
									} else {
										peptideNode = new PCQPeptideNode(this, key1,
												new Pair<QuantifiedPeptideInterface, PositionInPeptide>(peptide1,
														positionInPeptide1),
												new Pair<QuantifiedPeptideInterface, PositionInPeptide>(peptide2,
														positionInPeptide2));
									}
									peptideNodesByPeptides.put(peptide1, peptideNode);
									peptideNodesByPeptides.put(peptide2, peptideNode);
									// add the two peptides to the peptide node
									peptideNode.addQuantifiedPeptide(peptide1, positionInPeptide1);
									peptideNode.addQuantifiedPeptide(peptide2, positionInPeptide2);
									if (!peptideNode.getItemsInNode().isEmpty()) {
										// add to the set of nodes
										peptideNodes.add(peptideNode);
										// add to the map
										peptideNodesByPeptideNodeKey.put(key1, peptideNode);
									}

								} else {
									// they dont share any key create a peptide
									// node
									// for each peptide separated peptide node
									// for
									// psm1
									PCQPeptideNode peptideNode = null;
									if (peptideNodesByPeptideNodeKey.containsKey(key1)) {
										peptideNode = peptideNodesByPeptideNodeKey.get(key1);
									} else {
										peptideNode = new PCQPeptideNode(this, key1,
												new Pair<QuantifiedPeptideInterface, PositionInPeptide>(peptide1,
														positionInPeptide1));
									}
									peptideNodesByPeptides.put(peptide1, peptideNode);
									// add the two peptides to the peptide node
									peptideNode.addQuantifiedPeptide(peptide1, positionInPeptide1);
									if (!peptideNode.getItemsInNode().isEmpty()) {

										peptideNodes.add(peptideNode);
										peptideNodesByPeptideNodeKey.put(key1, peptideNode);
									}
									// protein node for protein2
									PCQPeptideNode peptideNode2 = null;
									if (peptideNodesByPeptideNodeKey.containsKey(key2)) {
										peptideNode2 = peptideNodesByPeptideNodeKey.get(key2);
									} else {
										peptideNode2 = new PCQPeptideNode(this, key2,
												new Pair<QuantifiedPeptideInterface, PositionInPeptide>(peptide2,
														positionInPeptide2));
									}
									peptideNodesByPeptides.put(peptide2, peptideNode2);

									peptideNode2.addQuantifiedPeptide(peptide2, positionInPeptide2);
									if (!peptideNode2.getItemsInNode().isEmpty()) {

										peptideNodes.add(peptideNode2);
										peptideNodesByPeptideNodeKey.put(key2, peptideNode2);
									}
								}
							}
						}
					} else {
						// in this case, it is not isobaric isotopologues, so we
						// have to take all the positionsInPeptide together
						// because each peptide will only have one ratio and if
						// that peptide has more than one position quantified,
						// it will be an ambiguous one
						// iterate over all the keys
						final List<PositionInProtein> proteinKeysFromPeptide1 = new ArrayList<PositionInProtein>();
						final List<PositionInProtein> proteinKeysFromPeptide2 = new ArrayList<PositionInProtein>();
						for (final PositionInPeptide positionInPeptide1 : proteinKeysByPeptide1Keys.keySet()) {
							proteinKeysFromPeptide1.addAll(proteinKeysByPeptide1Keys.get(positionInPeptide1));
						}
						for (final PositionInPeptide positionInPeptide2 : proteinKeysByPeptide2Keys.keySet()) {
							proteinKeysFromPeptide2.addAll(proteinKeysByPeptide2Keys.get(positionInPeptide2));
						}
						final String key1 = QuantUtils.getPositionsInProteinsKey(proteinKeysFromPeptide1,
								useProteinGeneName, useProteinID, uplr, getParams().getUniprotVersion());
						final String key2 = QuantUtils.getPositionsInProteinsKey(proteinKeysFromPeptide2,
								useProteinGeneName, useProteinID, uplr, getParams().getUniprotVersion());

						// now I compare the two list of keys
						// if they are equal, that means, if they share
						// the
						// same sites of the same proteins
						if (PCQUtils.areEquals(proteinKeysFromPeptide1, proteinKeysFromPeptide2)) {

							// key1 and key2 should be the same

							// collapsed in the same peptide node
							PCQPeptideNode peptideNode = null;
							// if the protein was already associated
							// with
							// some protein node
							if (peptideNodesByPeptideNodeKey.containsKey(key1)) {
								peptideNode = peptideNodesByPeptideNodeKey.get(key1);
							} else {
								peptideNode = new PCQPeptideNode(this, key1);
								if (proteinKeysByPeptide1Keys.isEmpty()) {
									log.info("asdf");
								}
								for (final PositionInPeptide positionInPeptide : proteinKeysByPeptide1Keys.keySet()) {
									peptideNode.addQuantifiedPeptide(peptide1, positionInPeptide);
								}
								if (proteinKeysByPeptide2Keys.isEmpty()) {
									log.info("asdf");
								}
								for (final PositionInPeptide positionInPeptide : proteinKeysByPeptide2Keys.keySet()) {
									peptideNode.addQuantifiedPeptide(peptide2, positionInPeptide);
								}
							}
							peptideNodesByPeptides.put(peptide1, peptideNode);
							peptideNodesByPeptides.put(peptide2, peptideNode);
							// add the two peptides to the peptide node
							for (final PositionInPeptide positionInPeptide : proteinKeysByPeptide1Keys.keySet()) {
								peptideNode.addQuantifiedPeptide(peptide1, positionInPeptide);
							}
							for (final PositionInPeptide positionInPeptide : proteinKeysByPeptide2Keys.keySet()) {
								peptideNode.addQuantifiedPeptide(peptide2, positionInPeptide);
							}
							if (!peptideNode.getItemsInNode().isEmpty()) {

								// add to the set of nodes
								peptideNodes.add(peptideNode);
								// add to the map
								peptideNodesByPeptideNodeKey.put(key1, peptideNode);
							}

						} else {
							// they dont share any key create a peptide node for
							// each peptide separated peptide node for psm1
							PCQPeptideNode peptideNode = null;
							if (peptideNodesByPeptideNodeKey.containsKey(key1)) {
								peptideNode = peptideNodesByPeptideNodeKey.get(key1);
							} else {
								peptideNode = new PCQPeptideNode(this, key1);
								for (final PositionInPeptide positionInPeptide : proteinKeysByPeptide1Keys.keySet()) {
									peptideNode.addQuantifiedPeptide(peptide1, positionInPeptide);
								}
							}
							peptideNodesByPeptides.put(peptide1, peptideNode);
							// add the two peptides to the peptide node
							for (final PositionInPeptide positionInPeptide : proteinKeysByPeptide1Keys.keySet()) {
								peptideNode.addQuantifiedPeptide(peptide1, positionInPeptide);
							}
							if (!peptideNode.getItemsInNode().isEmpty()) {

								peptideNodes.add(peptideNode);
								peptideNodesByPeptideNodeKey.put(key1, peptideNode);
							}
							// protein node for protein2
							PCQPeptideNode peptideNode2 = null;
							if (peptideNodesByPeptideNodeKey.containsKey(key2)) {
								peptideNode2 = peptideNodesByPeptideNodeKey.get(key2);
							} else {
								peptideNode2 = new PCQPeptideNode(this, key2);
								for (final PositionInPeptide positionInPeptide : proteinKeysByPeptide2Keys.keySet()) {
									peptideNode2.addQuantifiedPeptide(peptide2, positionInPeptide);
								}
							}
							peptideNodesByPeptides.put(peptide2, peptideNode2);

							for (final PositionInPeptide positionInPeptide : proteinKeysByPeptide2Keys.keySet()) {
								peptideNode2.addQuantifiedPeptide(peptide2, positionInPeptide);
							}
							if (!peptideNode2.getItemsInNode().isEmpty()) {

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
			final QuantifiedPeptideInterface peptide = peptides.iterator().next();
			final Map<PositionInPeptide, List<PositionInProtein>> proteinKeysByPeptideKeys = peptide
					.getProteinKeysByPeptideKeysForQuantifiedAAs(quantifiedAAs, uplr, PCQUtils.proteinSequences);

			for (final PositionInPeptide positionInPeptide : proteinKeysByPeptideKeys.keySet()) {
				final List<PositionInProtein> positionsInProtein = proteinKeysByPeptideKeys.get(positionInPeptide);

				final String key = QuantUtils.getPositionsInProteinsKey(positionsInProtein, useProteinGeneName,
						useProteinID, uplr, getParams().getUniprotVersion());
				PCQPeptideNode peptideNode = null;
				if (peptideNodesByPeptideNodeKey.containsKey(key)) {
					peptideNode = peptideNodesByPeptideNodeKey.get(key);
				} else {
					peptideNode = new PCQPeptideNode(this, key,
							new Pair<QuantifiedPeptideInterface, PositionInPeptide>(peptide, positionInPeptide));
				}
				peptideNode.addQuantifiedPeptide(peptide, positionInPeptide);
				if (!peptideNode.getItemsInNode().isEmpty()) {

					peptideNodes.add(peptideNode);
					peptideNodesByPeptideNodeKey.put(key, peptideNode);
				}
			}
		}

	}

	/**
	 * Iterates all peptides in the cluster over all others and creates a node for
	 * all peptides that have the same key, which, in this case is the position of
	 * the quantifiedAA in the protein.
	 * 
	 * @param quantifiedAAs
	 */
	private void createPeptideNodesByPTMs() {

		final UniprotProteinLocalRetriever uplr = PCQUtils
				.getUniprotProteinLocalRetrieverByFolder(getParams().getUniprotReleasesFolder());

		final boolean useProteinGeneName = getParams().getProteinLabel() == ProteinNodeLabel.GENE;
		final boolean useProteinID = getParams().getProteinLabel() == ProteinNodeLabel.ID;
		final List<QuantifiedPeptideInterface> peptides = new ArrayList<QuantifiedPeptideInterface>();
		peptides.addAll(individualQuantifiedPeptideSet);
		final Set<String> set = new THashSet<String>();
		for (final QuantifiedPeptideInterface pep : peptides) {
			if (set.contains(pep.getKey())) {
				log.warn("Inconsistency error. 2 peptides with the same key cannot exist!");
			} else {
				set.add(pep.getKey());
			}
		}
		// to not include the same peptide in different peptide nodes
		final Map<QuantifiedPeptideInterface, PCQPeptideNode> peptideNodesByPeptides = new THashMap<QuantifiedPeptideInterface, PCQPeptideNode>();
		if (peptides.size() > 1) {
			for (int i = 0; i < peptides.size(); i++) {
				final QuantifiedPeptideInterface peptide1 = peptides.get(i);

				// get the keys from the peptide.
				// not that the peptide could have more than one key because 2
				// reasons:
				// - it could be shared by more than one protein
				// - it could have more than one quantified aminoacid in its
				// sequence
				final List<PTMInProtein> ptmsInProteinFromPeptide1 = filterPTMsOfInterest(
						peptide1.getPTMsInProtein(uplr, PCQUtils.proteinSequences));
				for (int j = i + 1; j < peptides.size(); j++) {
					final QuantifiedPeptideInterface peptide2 = peptides.get(j);
					// get the keys from the peptide.
					// not that the peptide could have more than one key because
					// 2
					// reasons:
					// - it could be shared by more than one protein
					// - it could have more than one quantified aminoacid in its
					// sequence

					final List<PTMInProtein> ptmsInProteinFromPeptide2 = filterPTMsOfInterest(
							peptide2.getPTMsInProtein(uplr, PCQUtils.proteinSequences));

					// iterate over all the keys

					final String key1 = QuantUtils.getPositionsInProteinsKey(
							QuantUtils.getAsPositionInProtein(ptmsInProteinFromPeptide1),
							peptide1.getQuantifiedProteins().stream().map(p -> p.getAccession())
									.collect(Collectors.toSet()),
							useProteinGeneName, useProteinID, uplr, getParams().getUniprotVersion());

					final String key2 = QuantUtils.getPositionsInProteinsKey(
							QuantUtils.getAsPositionInProtein(ptmsInProteinFromPeptide2),
							peptide2.getQuantifiedProteins().stream().map(p -> p.getAccession())
									.collect(Collectors.toSet()),
							useProteinGeneName, useProteinID, uplr, getParams().getUniprotVersion());

					// now I compare the two list of keys
					// if they are equal, that means, if they share the
					// same sites of the same proteins
					if (PCQUtils.areEquals(QuantUtils.getAsPositionInProtein(ptmsInProteinFromPeptide1),
							QuantUtils.getAsPositionInProtein(ptmsInProteinFromPeptide2))) {

						// key1 and key2 should be the same

						// collapsed in the same peptide node
						PCQPeptideNode peptideNode = null;
						// if the protein was already associated with
						// some protein node
						if (peptideNodesByPeptideNodeKey.containsKey(key1)) {
							peptideNode = peptideNodesByPeptideNodeKey.get(key1);
						} else {
							peptideNode = new PCQPeptideNode(this, key1);
						}
						peptideNodesByPeptides.put(peptide1, peptideNode);
						peptideNodesByPeptides.put(peptide2, peptideNode);
						// add the two peptides to the peptide node
						peptideNode.addQuantifiedPeptide(peptide1);
						peptideNode.addQuantifiedPeptide(peptide2);

						// add to the set of nodes
						peptideNodes.add(peptideNode);
						// add to the map
						peptideNodesByPeptideNodeKey.put(key1, peptideNode);

					} else {
						// they dont share any key create a peptide node
						// for each peptide separated peptide node for
						// psm1
						PCQPeptideNode peptideNode = null;
						if (peptideNodesByPeptideNodeKey.containsKey(key1)) {
							peptideNode = peptideNodesByPeptideNodeKey.get(key1);
						} else {
							peptideNode = new PCQPeptideNode(this, key1);
						}
						peptideNodesByPeptides.put(peptide1, peptideNode);
						// add the two peptides to the peptide node
						peptideNode.addQuantifiedPeptide(peptide1);

						peptideNodes.add(peptideNode);
						peptideNodesByPeptideNodeKey.put(key1, peptideNode);

						// protein node for protein2
						PCQPeptideNode peptideNode2 = null;
						if (peptideNodesByPeptideNodeKey.containsKey(key2)) {
							peptideNode2 = peptideNodesByPeptideNodeKey.get(key2);
						} else {
							peptideNode2 = new PCQPeptideNode(this, key2);
						}
						peptideNodesByPeptides.put(peptide2, peptideNode2);
						peptideNode2.addQuantifiedPeptide(peptide2);

						peptideNodes.add(peptideNode2);
						peptideNodesByPeptideNodeKey.put(key2, peptideNode2);
					}
				}
			}

		} else {
			// only one peptide
			// create a peptide node for each position in the peptide
			final QuantifiedPeptideInterface peptide = peptides.iterator().next();
			final List<PTMInProtein> positionsInProtein = filterPTMsOfInterest(
					peptide.getPTMsInProtein(uplr, PCQUtils.proteinSequences));

			final String key = QuantUtils.getPositionsInProteinsKey(
					QuantUtils.getAsPositionInProtein(positionsInProtein), useProteinGeneName, useProteinID, uplr,
					getParams().getUniprotVersion());
			PCQPeptideNode peptideNode = null;
			if (peptideNodesByPeptideNodeKey.containsKey(key)) {
				peptideNode = peptideNodesByPeptideNodeKey.get(key);
			} else {
				peptideNode = new PCQPeptideNode(this, key);
			}
			peptideNode.addQuantifiedPeptide(peptide);

			peptideNodes.add(peptideNode);
			peptideNodesByPeptideNodeKey.put(key, peptideNode);

		}

	}

	/**
	 * Selects the {@link PTMInProtein} objects that are among the ones stated as
	 * quantified in the input parameter file
	 * 
	 * @param ptmsInProteins
	 * @return
	 */
	private List<PTMInProtein> filterPTMsOfInterest(List<PTMInProtein> ptmsInProteins) {
		if (getParams().isCollapseByPTMs()) {
			if (getParams().getPTMsQuantified() != null) {
				final List<PTM> quantifiedPTMs = getParams().getPTMsQuantified();

				final List<PTMInProtein> ret = new ArrayList<PTMInProtein>();
				for (final PTMInProtein ptmInProtein : ptmsInProteins) {
					boolean valid = false;
					for (final PTM ptm : quantifiedPTMs) {
						if (ptm.isEquivalent(ptmInProtein)) {
							valid = true;
							break;
						}
					}
					if (valid) {
						ret.add(ptmInProtein);
					}
				}
				return ret;
			}
		}
		return Collections.emptyList();
	}

	private void createProteinNodes() throws IOException {
		// create a map to store proteins by accession
		final Map<String, Set<QuantifiedProteinInterface>> proteinMapByKey = new THashMap<String, Set<QuantifiedProteinInterface>>();
		final Map<String, Set<QuantifiedProteinInterface>> proteinMapByAcc = new THashMap<String, Set<QuantifiedProteinInterface>>();
		for (final QuantifiedProteinInterface protein : individualQuantifiedProteinSet) {
			if (proteinMapByKey.containsKey(protein.getKey())) {
				proteinMapByKey.get(protein.getKey()).add(protein);
			} else {
				final Set<QuantifiedProteinInterface> set = new THashSet<QuantifiedProteinInterface>();
				set.add(protein);
				proteinMapByKey.put(protein.getKey(), set);
			}
			if (proteinMapByAcc.containsKey(protein.getAccession())) {
				proteinMapByAcc.get(protein.getAccession()).add(protein);
			} else {
				final Set<QuantifiedProteinInterface> set = new THashSet<QuantifiedProteinInterface>();
				set.add(protein);
				proteinMapByAcc.put(protein.getAccession(), set);
			}
		}

		final List<String> keyList = new ArrayList<String>();
		keyList.addAll(proteinMapByKey.keySet());
		if (proteinMapByKey.size() > 1) {
			for (int i = 0; i < keyList.size(); i++) {
				final String key1 = keyList.get(i);
				final Set<QuantifiedProteinInterface> proteins1 = proteinMapByKey.get(key1);
				for (int j = i + 1; j < keyList.size(); j++) {
					final String key2 = keyList.get(j);
					final Set<QuantifiedProteinInterface> proteins2 = proteinMapByKey.get(key2);

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
									for (final String key : proteinNodesByProteinKey.keySet()) {
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
			if (proteinMapByKey.isEmpty()) {
				log.info(this);
			}
			final Collection<QuantifiedProteinInterface> proteins = proteinMapByKey.values().iterator().next();
			final PCQProteinNode proteinNode = new PCQProteinNode(this, proteins);
			proteinNodes.add(proteinNode);
			proteinNodesByProteinKey.put(proteinMapByKey.keySet().iterator().next(), proteinNode);

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

		final NWResult[][] results = new NWResult[individualQuantifiedPeptideSet.size()][individualQuantifiedPeptideSet
				.size()];

		final List<QuantifiedPeptideInterface> pepList = new ArrayList<QuantifiedPeptideInterface>();
		pepList.addAll(individualQuantifiedPeptideSet);

		// two for loops: compare two different peptides
		for (int i = 0; i < pepList.size(); i++) {
			final QuantifiedPeptideInterface pep1 = pepList.get(i);
			for (int j = i + 1; j < pepList.size(); j++) {
				if (i == j) {
					// same peptide, move on
					continue;
				}
				// NWAlign (give back score)
				final QuantifiedPeptideInterface pep2 = pepList.get(j);
				final NWResult result = NWAlign.needlemanWunsch(pep1.getSequence(), pep2.getSequence());
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
		// (getParams().isCollapseIndistinguishableProteins())
		// {
		// createPairsCollapsingIndistinguisibleProteins(gAM);
		// } else {
		// createPairsWithoutCollapsing(gAM);
		// }
		createProteinPairs(peptideAlignments);
	}

	private void createProteinPairs(AlignmentSet peptideAlignments) {

		final List<PCQProteinNode> proteinNodeList = new ArrayList<PCQProteinNode>();
		proteinNodeList.addAll(proteinNodes);

		for (int i = 0; i < proteinNodeList.size(); i++) {
			for (int j = i + 1; j < proteinNodeList.size(); j++) {
				final PCQProteinNode proteinNode1 = proteinNodeList.get(i);

				final PCQProteinNode proteinNode2 = proteinNodeList.get(j);
				if (PCQUtils.shareAtLeastOnePeptideNode(proteinNode1, proteinNode2, false)
						|| (peptideAlignments != null && PCQUtils.shareAtLeastOnePeptideBySimilarity(proteinNode1,
								proteinNode2, peptideAlignments, false))) {
					final ProteinPair pair = new ProteinPair(proteinNode1, proteinNode2);
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
	 * Get all protein pairs, the regular ones and the discarded ones (containing a
	 * discarded protein node)
	 *
	 * @return
	 */
	public Set<ProteinPair> getAllProteinPairs() {
		final Set<ProteinPair> ret = new THashSet<ProteinPair>();
		ret.addAll(getProteinPairs());
		ret.addAll(getDiscardedProteinPairs());
		return ret;
	}

	@Override
	public String toString() {
		final StringBuilder sb = new StringBuilder();
		sb.append("PROTN[");
		for (final PCQProteinNode protein : PCQUtils.getSortedProteinNodesByAcc(proteinNodes)) {
			if (!"PROTN[".equals(sb.toString())) {
				sb.append(",");
			}
			sb.append(protein.getKey());
		}
		sb.append("]NTORP");
		final StringBuilder sb2 = new StringBuilder();
		sb2.append("PEPN[");
		for (final PCQPeptideNode peptide : PCQUtils.getSortedPeptideNodesBySequence(peptideNodes)) {
			if (!"PEPN[".equals(sb2.toString())) {
				sb2.append(",");
			}
			sb2.append(peptide.getKey());
		}
		sb2.append("]NPEP");
		final String string = sb.append("\t").append(sb2).toString();
		final StringBuilder sb3 = new StringBuilder();
		final List<QuantifiedProteinInterface> proteinList = new ArrayList<QuantifiedProteinInterface>();
		proteinList.addAll(getProteinSet());
		Collections.sort(proteinList, new Comparator<QuantifiedProteinInterface>() {
			@Override
			public int compare(QuantifiedProteinInterface o1, QuantifiedProteinInterface o2) {
				return o1.getKey().compareTo(o2.getKey());
			}
		});
		for (final QuantifiedProteinInterface protein : proteinList) {
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
		final QuantifiedPeptideInterface peptide = peptideNode.getQuantifiedPeptides().iterator().next();

		final Set<PCQPeptideNode> ret = new THashSet<PCQPeptideNode>();

		for (final AlignedPeptides alignment : alignmentSet.getAlignmentsForPeptide(peptide)) {
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
		for (final AlignedPeptides alignment : alignments) {
			if (alignment.getPeptideAligned(peptide1).equals(peptide2)) {
				return alignment.getAlignmentResult();
			}
		}

		return null;
	}

	public Set<String> getPeptideNodeKeys() {
		final Set<String> ret = new THashSet<String>();
		for (final PCQPeptideNode peptideNode : peptideNodes) {
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
		final Set<String> set = new THashSet<String>();
		for (final PCQPeptideNode peptideNode : getNonDiscardedPeptideNodes()) {
			for (final QuantifiedPeptideInterface peptide : peptideNode.getItemsInNode()) {
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
		final Set<String> set = new THashSet<String>();
		for (final PCQPeptideNode peptideNode : getPeptideNodes()) {
			if (peptideNode.isDiscarded()) {
				for (final QuantifiedPeptideInterface peptide : peptideNode.getItemsInNode()) {
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
		final Set<String> set = new THashSet<String>();
		for (final PCQProteinNode proteinNode : getNonDiscardedProteinNodes()) {
			for (final QuantifiedProteinInterface protein : proteinNode.getItemsInNode()) {
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

		for (final PCQPeptideNode peptideNode : peptideNodes) {
			if (peptideNode.getKey().equals(peptideNodeID)) {
				return peptideNode;
			}
		}

		return null;
	}

	public Set<PCQProteinNode> getNonDiscardedProteinNodes() {
		final Set<PCQProteinNode> ret = new THashSet<PCQProteinNode>();
		for (final PCQProteinNode proteinNode : getProteinNodes()) {
			if (!proteinNode.isDiscarded()) {
				ret.add(proteinNode);
			}
		}
		return ret;
	}

	public Set<PCQPeptideNode> getNonDiscardedPeptideNodes() {
		final Set<PCQPeptideNode> ret = new THashSet<PCQPeptideNode>();
		for (final PCQPeptideNode peptideNode : getPeptideNodes()) {
			if (!peptideNode.isDiscarded()) {
				ret.add(peptideNode);
			} else {
				// log.info("asdf");
			}
		}
		return ret;
	}

	public Set<QuantifiedProteinInterface> getNonDiscardedProteinSet() {
		final Set<QuantifiedProteinInterface> ret = new THashSet<QuantifiedProteinInterface>();
		final Set<QuantifiedProteinInterface> proteinSet = getProteinSet();
		for (final QuantifiedProteinInterface protein : proteinSet) {
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
				for (final QuantifiedPSMInterface psm : quantifiedPSMs) {
					final Iterator<QuantifiedProteinInterface> proteinIterator = psm.getQuantifiedProteins().iterator();
					while (proteinIterator.hasNext()) {
						final QuantifiedProteinInterface protein = proteinIterator.next();

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

	// @Override
	// public int hashCode() {
	// return Objects.hash(getClusterID());
	// }
}
