package edu.scripps.yates.pcq.filter;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.pcq.ProteinClusterQuantParameters;
import edu.scripps.yates.pcq.model.PCQPeptideNode;
import edu.scripps.yates.pcq.model.PCQProteinNode;
import edu.scripps.yates.pcq.model.ProteinCluster;

public abstract class PCQFilter {
	protected final static Logger log = Logger.getLogger(PCQFilter.class);
	private static Set<PCQPeptideNode> staticDiscardedPeptideNodesForStatistics = new HashSet<PCQPeptideNode>();
	private static Set<PCQProteinNode> staticDiscardedProteinNodesForStatistics = new HashSet<PCQProteinNode>();

	public void filter(ProteinCluster cluster) {
		final Iterator<PCQProteinNode> proteinNodesIterator = cluster.getProteinNodes().iterator();
		Set<PCQProteinNode> discardedProteinNodes = new HashSet<PCQProteinNode>();
		// final int originalSize = cluster.getProteinNodes().size();
		// log.debug("Filtering " + originalSize + " protein nodes");
		while (proteinNodesIterator.hasNext()) {
			PCQProteinNode pcqProteinNode = proteinNodesIterator.next();
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
		Set<PCQPeptideNode> discardedPeptideNodes = new HashSet<PCQPeptideNode>();
		final Iterator<PCQPeptideNode> peptideNodesIterator = cluster.getPeptideNodes().iterator();
		// final int originalSize2 = cluster.getPeptideNodes().size();
		// log.debug("Filtering " + originalSize2 + " peptide nodes");
		while (peptideNodesIterator.hasNext()) {
			PCQPeptideNode pcqPeptideNode = peptideNodesIterator.next();
			final boolean valid = filter(pcqPeptideNode);
			if (!valid) {
				pcqPeptideNode.setDiscarded(true);
				staticDiscardedPeptideNodesForStatistics.add(pcqPeptideNode);
				if (ProteinClusterQuantParameters.getInstance().isRemoveFilteredNodes()) {
					peptideNodesIterator.remove();
					discardedPeptideNodes.add(pcqPeptideNode);
				}
			}
		}

		// this is only populated if params.isRemoveFilteredNodes() is true
		for (PCQPeptideNode pcqPeptideNode : discardedPeptideNodes) {
			// remove individual peptides from cluster
			final Set<QuantifiedPeptideInterface> individualPeptides = pcqPeptideNode.getQuantifiedPeptides();
			for (QuantifiedPeptideInterface peptide : individualPeptides) {
				cluster.getPeptideSet().remove(peptide);
			}
			if (ProteinClusterQuantParameters.getInstance().isRemoveFilteredNodes()) {
				pcqPeptideNode.removePeptidesFromProteinsInNode();
			}
		}
	}

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
