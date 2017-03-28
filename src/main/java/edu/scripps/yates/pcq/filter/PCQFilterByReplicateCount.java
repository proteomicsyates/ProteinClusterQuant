package edu.scripps.yates.pcq.filter;

import org.apache.log4j.Logger;

import edu.scripps.yates.pcq.model.PCQPeptideNode;
import edu.scripps.yates.pcq.model.PCQProteinNode;

public class PCQFilterByReplicateCount extends PCQFilter {
	private final int minReplicatesPerNode;
	private final static Logger log = Logger.getLogger(PCQFilterByReplicateCount.class);

	public PCQFilterByReplicateCount(int minReplicatesPerNode) {
		this.minReplicatesPerNode = minReplicatesPerNode;
		log.info("Creating filter by replicate count=" + minReplicatesPerNode);

	}

	@Override
	protected boolean filter(PCQProteinNode proteinNode) {
		return true;

	}

	@Override
	protected boolean filter(PCQPeptideNode peptideNode) {
		final int minReplicatesPerNode = peptideNode.getFileNames().size();

		if (minReplicatesPerNode >= this.minReplicatesPerNode) {
			return true;
		}
		return false;
	}

	@Override
	protected boolean filterNonQuantifiedNodes() {
		return false;
	}

}
