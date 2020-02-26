package edu.scripps.yates.pcq.util;

import java.util.ArrayList;

import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;

public class DiscardedPeptidesSet extends ArrayList<DiscardedPeptide> {
	/**
	 * 
	 */
	private static final long serialVersionUID = -7836944653117423789L;

	public static enum DISCARD_REASON {
		AMBIGUOS_QUANT, CONTAINING_PTM, NOT_CONTAINING_PTM, PEPTIDE_WITH_NO_PROTEIN, DISCARDED_BY_FILTER,
		PEPTIDE_WITH_NO_QUANT_SITE
	};

	private static DiscardedPeptidesSet instance;

	private DiscardedPeptidesSet() {
		super();
	}

	public static DiscardedPeptidesSet getInstance() {
		if (instance == null) {
			instance = new DiscardedPeptidesSet();
		}
		return instance;
	}

	public void add(QuantifiedPeptideInterface peptide, DISCARD_REASON reason) {
		add(peptide, reason, null);
	}

	public void add(QuantifiedPeptideInterface peptide, DISCARD_REASON reason, String additionalDescription) {
		final DiscardedPeptide d = new DiscardedPeptide(peptide, reason);
		if (additionalDescription != null) {
			d.setAdditionalDescription(additionalDescription);
		}
		add(d);
	}

}
