package edu.scripps.yates.pcq.xgmml.util;

import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.utilities.alignment.nwalign.NWResult;

public class AlignedPeptides {
	private final NWResult alignmentResult;
	private final QuantifiedPeptideInterface peptide1;
	private final QuantifiedPeptideInterface peptide2;

	public AlignedPeptides(NWResult alignmentResult, QuantifiedPeptideInterface peptide1,
			QuantifiedPeptideInterface peptide2) {
		super();
		this.alignmentResult = alignmentResult;
		this.peptide1 = peptide1;
		this.peptide2 = peptide2;
	}

	/**
	 * @return the alignmentResult
	 */
	public NWResult getAlignmentResult() {
		return alignmentResult;
	}

	/**
	 * @return the peptide1
	 */
	public QuantifiedPeptideInterface getPeptide1() {
		return peptide1;
	}

	/**
	 * @return the peptide2
	 */
	public QuantifiedPeptideInterface getPeptide2() {
		return peptide2;
	}

	public QuantifiedPeptideInterface getPeptideAligned(QuantifiedPeptideInterface peptide) {
		if (peptide.equals(peptide1)) {
			return peptide2;
		} else if (peptide.equals(peptide2)) {
			return peptide1;
		}
		return null;
	}

	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "AlignedPeptides [alignmentResult=" + alignmentResult + ", peptide1=" + peptide1 + ", peptide2="
				+ peptide2 + "]";
	}
}
