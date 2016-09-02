package edu.scripps.yates.pcq.filter;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.pcq.model.PCQPeptideNode;
import edu.scripps.yates.pcq.model.PCQProteinNode;
import edu.scripps.yates.pcq.model.ProteinCluster;

public abstract class PCQFilter {
	protected final static Logger log = Logger.getLogger(PCQFilter.class);
	private static Set<PCQPeptideNode> staticDiscardedPeptideNodes = new HashSet<PCQPeptideNode>();
	private static Set<PCQProteinNode> staticDiscardedProteinNodes = new HashSet<PCQProteinNode>();

	public void filter(ProteinCluster cluster) {
		final Iterator<PCQProteinNode> proteinNodesIterator = cluster.getProteinNodes().iterator();
		Set<PCQProteinNode> discardedProteinNodes = new HashSet<PCQProteinNode>();
		// final int originalSize = cluster.getProteinNodes().size();
		// log.debug("Filtering " + originalSize + " protein nodes");
		while (proteinNodesIterator.hasNext()) {
			PCQProteinNode pcqProteinNode = proteinNodesIterator.next();
			final boolean valid = filter(pcqProteinNode);
			if (!valid) {
				proteinNodesIterator.remove();
				discardedProteinNodes.add(pcqProteinNode);
				staticDiscardedProteinNodes.add(pcqProteinNode);
			}
		}
		for (PCQProteinNode pcqProteinNode : discardedProteinNodes) {
			// remove individual proteins from cluster
			final Set<QuantifiedProteinInterface> individualProteins = pcqProteinNode.getIndividualProteins();
			for (QuantifiedProteinInterface protein : individualProteins) {
				cluster.getProteinSet().remove(protein);
				// set as invalid for not use them in the QuantAnalysis
				protein.setDiscarded(true);
			}
			pcqProteinNode.disconnectProteinsInNode();
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
				peptideNodesIterator.remove();
				discardedPeptideNodes.add(pcqPeptideNode);
				staticDiscardedPeptideNodes.add(pcqPeptideNode);
			}
		}
		// if (cluster.getPeptideNodes().size() != originalSize2) {
		// log.info(originalSize2 - cluster.getPeptideNodes().size() + " peptide
		// nodes where discarded");
		// }
		for (PCQPeptideNode pcqPeptideNode : discardedPeptideNodes) {
			// remove individual peptides from cluster
			final Set<QuantifiedPeptideInterface> individualPeptides = pcqPeptideNode.getIndividualPeptides();
			for (QuantifiedPeptideInterface peptide : individualPeptides) {
				cluster.getPeptideSet().remove(peptide);
				// set as invalid for not use them in the QuantAnalysis
				peptide.setDiscarded(true);
			}
			pcqPeptideNode.disconnectPeptidesInNode();
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

	public static Set<PCQPeptideNode> getDiscardedPeptideNodes() {
		return staticDiscardedPeptideNodes;
	}

	public static Set<PCQProteinNode> getDiscardedProteinNodes() {
		return staticDiscardedProteinNodes;
	}
}
