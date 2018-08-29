package edu.scripps.yates.pcq.model;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.CensusRatio;
import edu.scripps.yates.census.read.model.QuantifiedPSM;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.dtaselectparser.util.DTASelectModification;
import edu.scripps.yates.dtaselectparser.util.DTASelectPSM;
import edu.scripps.yates.utilities.model.enums.AggregationLevel;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.sequence.PTMInPeptide;

public class NonQuantifiedPSM extends QuantifiedPSM {
	private final DTASelectPSM psm;

	public NonQuantifiedPSM(DTASelectPSM psm) throws NumberFormatException, IOException {
		super(psm.getSequence().getRawSequence(), null, null, Integer.valueOf(psm.getScan()), psm.getChargeState(),
				psm.getRawFileName(), false);
		this.psm = psm;
	}

	@Override
	public Float getCalcMHplus() {
		return psm.getCalcMh().floatValue();
	}

	@Override
	public Float getMHplus() {
		return psm.getMh().floatValue();
	}

	@Override
	public boolean containsPTMs() {
		return psm.getModifications() != null && !psm.getModifications().isEmpty();
	}

	@Override
	public List<PTMInPeptide> getPtms() {
		final List<PTMInPeptide> ret = new ArrayList<PTMInPeptide>();
		final List<DTASelectModification> modifications = psm.getModifications();
		if (modifications != null) {
			for (final DTASelectModification dtaSelectModification : modifications) {
				final PTMInPeptide ptm = new PTMInPeptide(dtaSelectModification.getModPosition(),
						dtaSelectModification.getAa(), getSequence(), dtaSelectModification.getModificationShift());

				ret.add(ptm);
			}
		}
		return ret;
	}

	@Override
	public Set<QuantRatio> getRatios() {
		return Collections.emptySet();
	}

	@Override
	public QuantRatio getConsensusRatio(QuantCondition quantConditionNumerator,
			QuantCondition quantConditionDenominator) {
		return CensusRatio.getNaNRatio(quantConditionNumerator, quantConditionDenominator, AggregationLevel.PSM,
				"RATIO");
	}

	@Override
	public Set<QuantRatio> getNonInfinityRatios() {
		return Collections.emptySet();
	}

	@Override
	public void addRatio(QuantRatio ratio) {
		throw new IllegalArgumentException("addRatio not available for " + getClass().getSimpleName());

	}

	@Override
	public double getMeanRatios(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator) {
		throw new IllegalArgumentException("getMeanRatios not available for " + getClass().getSimpleName());

	}

	@Override
	public double getSTDRatios(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator) {
		throw new IllegalArgumentException("getSTDRatios not available for " + getClass().getSimpleName());

	}

	@Override
	public Set<Amount> getAmounts() {
		return Collections.emptySet();
	}

	@Override
	public void addAmount(Amount amount) {
		throw new IllegalArgumentException("addAmount not available for " + getClass().getSimpleName());

	}

	@Override
	public Float getDeltaCN() {
		return psm.getDeltacn().floatValue();
	}

	@Override
	public Float getXcorr() {
		return psm.getXcorr().floatValue();
	}

	@Override
	public Float getDeltaMass() {
		if (getCalcMHplus() != null && getMHplus() != null) {
			return getCalcMHplus() - getMHplus();
		}
		return null;
	}

	@Override
	public Double getMaxPeak() {
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.scripps.yates.census.read.model.interfaces.QuantifiedItem#
	 * isQuantified()
	 */
	@Override
	public boolean isQuantified() {
		return false;
	}
}
