package edu.scripps.yates.pcq.model;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.CensusRatio;
import edu.scripps.yates.census.read.model.QuantifiedPSM;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.proteomicsmodel.PSM;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AggregationLevel;
import edu.scripps.yates.utilities.sequence.PTMInPeptide;

public class NonQuantifiedPSM extends QuantifiedPSM {
	private final PSM psm;

	public NonQuantifiedPSM(PSM psm) throws NumberFormatException, IOException {
		super(psm.getFullSequence(), null, null, Integer.valueOf(psm.getScanNumber()), psm.getChargeState(),
				psm.getMSRun().getRunId(), false);
		this.psm = psm;
	}

	@Override
	public Float getCalcMH() {

		return psm.getCalcMH();
	}

	@Override
	public Float getExperimentalMH() {

		return psm.getExperimentalMH();
	}

	@Override
	public boolean containsPTMs() {
		return psm.containsPTMs();
	}

	@Override
	public List<PTMInPeptide> getPTMsInPeptide() {
		return psm.getPTMsInPeptide();
	}

	@Override
	public Set<QuantRatio> getQuantRatios() {
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
	public boolean addQuantRatio(QuantRatio ratio) {
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
	public boolean addAmount(Amount amount) {
		throw new IllegalArgumentException("addAmount not available for " + getClass().getSimpleName());

	}

	@Override
	public Float getDeltaCn() {
		return psm.getDeltaCn();
	}

	@Override
	public Float getXCorr() {
		return psm.getXCorr();
	}

	@Override
	public Float getMassErrorPPM() {

		if (getCalcMH() != null && getExperimentalMH() != null) {
			return getCalcMH() - getExperimentalMH();
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
