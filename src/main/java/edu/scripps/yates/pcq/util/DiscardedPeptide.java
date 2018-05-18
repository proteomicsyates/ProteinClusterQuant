package edu.scripps.yates.pcq.util;

import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.pcq.util.DiscardedPeptidesSet.DISCARD_REASON;

public class DiscardedPeptide {
	private final QuantifiedPeptideInterface peptide;
	private final DiscardedPeptidesSet.DISCARD_REASON reason;
	private String additionalDescription;

	public DiscardedPeptide(QuantifiedPeptideInterface peptide, DISCARD_REASON reason) {
		super();
		this.peptide = peptide;
		this.reason = reason;
	}

	public QuantifiedPeptideInterface getPeptide() {
		return peptide;
	}

	public DiscardedPeptidesSet.DISCARD_REASON getReason() {
		return reason;
	}

	public void setAdditionalDescription(String additionalDescription) {
		this.additionalDescription = additionalDescription;
	}

	public String getAdditionalDescription() {
		return additionalDescription;
	}
}
