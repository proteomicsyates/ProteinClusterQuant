package edu.scripps.yates.pcq.filter;

import java.util.Iterator;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.pcq.model.PCQPeptideNode;
import edu.scripps.yates.pcq.model.PCQProteinNode;
import edu.scripps.yates.pcq.model.ProteinCluster;
import edu.scripps.yates.pcq.params.ProteinClusterQuantParameters;
import gnu.trove.set.hash.THashSet;

public abstract class PCQFilter {
	protected final static Logger log = Logger.getLogger(PCQFilter.class);
	private static Set<PCQPeptideNode> staticDiscardedPeptideNodesForStatistics = new THashSet<PCQPeptideNode>();
	private static Set<PCQProteinNode> staticDiscardedProteinNodesForStatistics = new THashSet<PCQProteinNode>();

	public void filter(ProteinCluster cluster) {
		final Set<PCQProteinNode> proteinNodes = cluster.getProteinNodes();
		final Iterator<PCQProteinNode> proteinNodesIterator = proteinNodes.iterator();
		Set<PCQProteinNode> discardedProteinNodes = new THashSet<PCQProteinNode>();
		// final int originalSize = cluster.getProteinNodes().size();
		// log.debug("Filtering " + originalSize + " protein nodes");
		while (proteinNodesIterator.hasNext()) {
			PCQProteinNode pcqProteinNode = proteinNodesIterator.next();
			if (!filterNonQuantifiedNodes() && !pcqProteinNode.isQuantified()) {
				continue;
			}
			final boolean valid = filter(pcqProteinNode);
			if (!valid) {
				pcqProteinNode.setDiscarded(true);
				staticDiscardedProteinNodesForStatistics.add(pcqProteinNode);
				if (ProteinClusterQuantParameters.getInstance().isRemoveFilteredNodes()) {
					proteinNodesIterator.remove();
					discardedProteinNodes.add(pcqProteinNode);
				}
			}
		}

		// this is only populated if params.isRemoveFilteredNodes is true
		for (PCQProteinNode pcqProteinNode : discardedProteinNodes) {
			// remove individual proteins from cluster
			final Set<QuantifiedProteinInterface> individualProteins = pcqProteinNode.getQuantifiedProteins();
			for (QuantifiedProteinInterface protein : individualProteins) {
				cluster.getProteinSet().remove(protein);
			}
			pcqProteinNode.removeProteinsFromPeptidesInNode();
		}

		// if (cluster.getProteinNodes().size() != originalSize) {
		// log.info(originalSize - cluster.getProteinNodes().size() + " protein
		// nodes where discarded");
		// }
		Set<PCQPeptideNode> discardedPeptideNodes = new THashSet<PCQPeptideNode>();
		final Iterator<PCQPeptideNode> peptideNodesIterator = cluster.getPeptideNodes().iterator();
		// final int originalSize2 = cluster.getPeptideNodes().size();
		// log.debug("Filtering " + originalSize2 + " peptide nodes");
		while (peptideNodesIterator.hasNext()) {
			PCQPeptideNode pcqPeptideNode = peptideNodesIterator.next();
			if (!filterNonQuantifiedNodes() && !pcqPeptideNode.isQuantified()) {
				continue;
			}
			final boolean valid = filter(pcqPeptideNode);
			if (!valid) {
				if (pcqPeptideNode.getKey().equals("ADRDEASPYAAMLAAQDVAEK")) {
					log.info(pcqPeptideNode);
				}
				pcqPeptideNode.setDiscarded(true);
				staticDiscardedPeptideNodesForStatistics.add(pcqPeptideNode);
				if (ProteinClusterQuantParameters.getInstance().isRemoveFilteredNodes()) {
					// remove peptide node from cluster
					peptideNodesIterator.remove();
					discardedPeptideNodes.add(pcqPeptideNode);
				}
			}
		}

		// this is only populated if params.isRemoveFilteredNodes() is true
		for (PCQPeptideNode peptideNode : discardedPeptideNodes) {
			// remove individual peptides from cluster
			final Iterator<QuantifiedPeptideInterface> peptidesFromPeptideNode = peptideNode.getQuantifiedPeptides()
					.iterator();
			while (peptidesFromPeptideNode.hasNext()) {
				QuantifiedPeptideInterface peptide = peptidesFromPeptideNode.next();
				final Iterator<QuantifiedPSMInterface> psmsFromPeptide = peptide.getQuantifiedPSMs().iterator();
				while (psmsFromPeptide.hasNext()) {
					final QuantifiedPSMInterface psm = psmsFromPeptide.next();
					final Iterator<QuantifiedProteinInterface> proteinsFromPSM = psm.getQuantifiedProteins().iterator();
					while (proteinsFromPSM.hasNext()) {
						final QuantifiedProteinInterface protein = proteinsFromPSM.next();
						// remove psm from its protein
						protein.getQuantifiedPSMs().remove(psm);
						// remove protein from psm
						proteinsFromPSM.remove();
						// if protein has no psm, remove protein from cluster
						// and remove it from its protein node
						if (protein.getQuantifiedPSMs().isEmpty()) {
							cluster.getProteinSet().remove(protein);
							final Iterator<PCQProteinNode> proteinNodesFromCluster = cluster.getProteinNodes()
									.iterator();
							while (proteinNodesFromCluster.hasNext()) {
								final PCQProteinNode proteinNode = proteinNodesFromCluster.next();
								// remove protein from protein node
								proteinNode.getItemsInNode().remove(protein);
								// remove the link between protein and protein
								// node in the cluster
								cluster.getProteinNodesByProteinKey().remove(protein.getKey());
								// remove peptide node from protein node
								proteinNode.getPeptideNodes().remove(peptideNode);
								// remove protein node from peptide node
								peptideNode.getProteinNodes().remove(proteinNode);
								// if protein node has no proteins, remove it
								// from cluster
								if (proteinNode.getItemsInNode().isEmpty()) {
									proteinNodesFromCluster.remove();
								}
							}

						}
					}
					if ("ADRDEASPYAAMLAAQDVAEK".equals(peptide.getKey())) {
						log.info(peptide);
					}
					// remove psm from its peptide
					psmsFromPeptide.remove();
					// remove peptide of the psm
					psm.setQuantifiedPeptide(null, false);
				}
				// remove peptide from peptide node
				peptidesFromPeptideNode.remove();
				// remove peptide node from its connected protein nodes and
				// viceversa
				final Iterator<PCQProteinNode> proteinNodesFromPeptideNode = peptideNode.getProteinNodes().iterator();
				while (proteinNodesFromPeptideNode.hasNext()) {
					PCQProteinNode proteinNode = proteinNodesFromPeptideNode.next();
					proteinNode.getPeptideNodes().remove(peptideNode);
					// if protein node has no peptide nodes, remove from cluster
					if (proteinNode.getPeptideNodes().isEmpty()) {
						cluster.getProteinNodes().remove(proteinNode);
					}
					proteinNodesFromPeptideNode.remove();
				}
				// remove peptide from cluster
				cluster.getPeptideSet().remove(peptide);
				// remove peptide node from cluster
				cluster.getPeptideNodes().remove(peptideNode);
				// remove link between peptide node and peptide sequence in the
				// cluster
				cluster.getPeptideNodesByPeptideSequence().remove(peptide.getKey());
			}
			// if
			// (ProteinClusterQuantParameters.getInstance().isRemoveFilteredNodes())
			// {
			// pcqPeptideNode.removePeptidesFromProteinsInNode();
			// }
		}

		// check if all the individual proteins have psms or not
		// if so, remove those proteins
		cluster.removeIndividualProteinsWithNoPSMs();
	}

	/**
	 * States wether this filter is also applied to NonQuantifiedNodes or not.
	 *
	 * @return
	 */
	protected abstract boolean filterNonQuantifiedNodes();

	/**
	 * Returns true if the {@link PCQProteinNode} pass the filter and false
	 * otherwise
	 *
	 * @param proteinNode
	 * @return
	 */
	protected abstract boolean filter(PCQProteinNode proteinNode);

	/**
	 * Returns true if the {@link PCQPeptideNode} pass the filetr and false
	 * otherwise
	 *
	 * @param peptideNode
	 * @return
	 */
	protected abstract boolean filter(PCQPeptideNode peptideNode);

	/**
	 * Get discardedPeptideNodes (for statistics)
	 *
	 * @return
	 */
	public static Set<PCQPeptideNode> getDiscardedPeptideNodes() {
		return staticDiscardedPeptideNodesForStatistics;
	}

	/**
	 * Get discardedProteinNodes (for statistics)
	 *
	 * @return
	 */
	public static Set<PCQProteinNode> getDiscardedProteinNodes() {
		return staticDiscardedProteinNodesForStatistics;
	}
}
