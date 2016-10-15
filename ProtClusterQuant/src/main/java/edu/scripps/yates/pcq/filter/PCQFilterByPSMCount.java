package edu.scripps.yates.pcq.filter;

import org.apache.log4j.Logger;

import edu.scripps.yates.pcq.model.PCQPeptideNode;
import edu.scripps.yates.pcq.model.PCQProteinNode;

public class PCQFilterByPSMCount extends PCQFilter {
	private final int minPSMCountPerNode;
	private final static Logger log = Logger.getLogger(PCQFilterByPSMCount.class);

	public PCQFilterByPSMCount(int minPSMCountPerNode) {
		this.minPSMCountPerNode = minPSMCountPerNode;
		log.info("Creating filter by psm count=" + minPSMCountPerNode);
	}

	@Override
	protected boolean filter(PCQProteinNode proteinNode) {
		return true;

	}

	@Override
	protected boolean filter(PCQPeptideNode peptideNode) {
		final int psmCount = peptideNode.getQuantifiedPSMs().size();
		if (psmCount >= minPSMCountPerNode) {
			return true;
		}
		return false;
	}

	@Override
	protected boolean filterNonQuantifiedNodes() {
		return false;
	}

}
