package edu.scripps.yates.pcq.model;

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

import org.apache.log4j.Logger;

import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.pcq.ProteinClusterQuantParameters;
import edu.scripps.yates.pcq.filter.PCQFilter;
import edu.scripps.yates.pcq.util.PCQUtils;
import edu.scripps.yates.pcq.xgmml.util.AlignedPeptides;
import edu.scripps.yates.pcq.xgmml.util.AlignmentSet;
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
	// proteinPairs with discarded proteins
	private final Set<ProteinPair> discardedProteinPairs = new HashSet<ProteinPair>();
	private final AlignmentSet alignmentSet = new AlignmentSet();
	private Map<String, QuantifiedPeptideInterface> peptideMap;

	private final HashMap<String, PCQProteinNode> proteinNodesByProteinKey = new HashMap<String, PCQProteinNode>();
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
	public void createNodes() {
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
		Set<String> set = new HashSet<String>();
		for (QuantifiedPeptideInterface pep : peptides) {
			if (set.contains(pep.getKey())) {
				log.info("error");
			} else {
				set.add(pep.getKey());
			}
		}
		if (peptides.size() > 1) {
			for (int i = 0; i < peptides.size(); i++) {
				QuantifiedPeptideInterface peptide1 = peptides.get(i);
				for (int j = i + 1; j < peptides.size(); j++) {
					QuantifiedPeptideInterface peptide2 = peptides.get(j);
					if (peptide1.getKey().startsWith("EFMDDS") || peptide2.getKey().startsWith("EFMDDS")) {
						log.info(peptide1.getKey() + "  " + peptide2.getKey());
					}
					if (getParams().isCollapseIndistinguishablePeptides()
							&& PCQUtils.peptidesShareAllProteins(peptide1, peptide2)) {

						PCQPeptideNode peptideNode = null;
						PCQPeptideNode peptideNode2 = null;
						// if the protein was already associated with some
						// protein
						// node
						if (peptideNodesByPeptideSequence.containsKey(peptide1.getKey())) {
							peptideNode = peptideNodesByPeptideSequence.get(peptide1.getKey());
						}
						if (peptideNodesByPeptideSequence.containsKey(peptide2.getKey())) {
							peptideNode2 = peptideNodesByPeptideSequence.get(peptide2.getKey());
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
						peptideNodesByPeptideSequence.put(peptide1.getKey(), peptideNode);
						peptideNodesByPeptideSequence.put(peptide2.getKey(), peptideNode);

					} else {
						// create a peptide node for each peptide separated
						// peptide node for psm1
						PCQPeptideNode peptideNode = null;
						if (peptideNodesByPeptideSequence.containsKey(peptide1.getKey())) {
							peptideNode = peptideNodesByPeptideSequence.get(peptide1.getKey());
						} else {
							peptideNode = new PCQPeptideNode(this, peptide1);
						}
						peptideNodes.add(peptideNode);
						peptideNodesByPeptideSequence.put(peptide1.getKey(), peptideNode);

						// protein node for protein2
						PCQPeptideNode peptideNode2 = null;
						if (peptideNodesByPeptideSequence.containsKey(peptide2.getKey())) {
							peptideNode2 = peptideNodesByPeptideSequence.get(peptide2.getKey());
						} else {
							peptideNode2 = new PCQPeptideNode(this, peptide2);
						}
						peptideNodes.add(peptideNode2);
						peptideNodesByPeptideSequence.put(peptide2.getKey(), peptideNode2);
					}
				}
			}
		} else {
			// only one peptide
			// create a peptide node for the peptide
			QuantifiedPeptideInterface peptide = peptides.iterator().next();
			PCQPeptideNode peptideNode = new PCQPeptideNode(this, peptide);
			peptideNodes.add(peptideNode);
			peptideNodesByPeptideSequence.put(peptide.getKey(), peptideNode);
		}

	}

	private void createProteinNodes() {
		// create a map to store proteins by accession
		Map<String, Set<QuantifiedProteinInterface>> proteinMap = new HashMap<String, Set<QuantifiedProteinInterface>>();
		for (QuantifiedProteinInterface protein : individualQuantifiedProteinSet) {
			if (proteinMap.containsKey(protein.getKey())) {
				proteinMap.get(protein.getKey()).add(protein);
			} else {
				Set<QuantifiedProteinInterface> set = new HashSet<QuantifiedProteinInterface>();
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
		Set<ProteinPair> ret = new HashSet<ProteinPair>();
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

		Set<PCQPeptideNode> ret = new HashSet<PCQPeptideNode>();

		for (AlignedPeptides alignment : alignmentSet.getAlignmentsForPeptide(peptide)) {
			ret.add(peptideNodesByPeptideSequence.get(alignment.getPeptide1().getKey()));
			ret.add(peptideNodesByPeptideSequence.get(alignment.getPeptide2().getKey()));
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
				if (peptide.getQuantifiedProteins().isEmpty()) {
					log.info(peptide);
				}
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

	public HashMap<String, PCQProteinNode> getProteinNodesByProteinKey() {
		return proteinNodesByProteinKey;

	}

	/**
	 * @return the peptideNodesByPeptideSequence
	 */
	public HashMap<String, PCQPeptideNode> getPeptideNodesByPeptideSequence() {
		return peptideNodesByPeptideSequence;
	}

}
