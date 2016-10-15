package edu.scripps.yates.pcq.model;

import java.util.Collections;
import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.CensusRatio;
import edu.scripps.yates.census.read.model.QuantifiedProtein;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.dtaselectparser.util.DTASelectProtein;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;

public class NonQuantifiedProtein extends QuantifiedProtein {
	private final DTASelectProtein protein;

	public NonQuantifiedProtein(DTASelectProtein protein) {
		super(FastaParser.getACC(protein.getLocus()).getFirstelement());
		this.protein = protein;
	}

	@Override
	public Set<QuantRatio> getRatios() {
		return Collections.emptySet();
	}

	@Override
	public QuantRatio getConsensusRatio(QuantCondition quantConditionNumerator,
			QuantCondition quantConditionDenominator) {
		return CensusRatio.NAN_RATIO;
	}

	@Override
	public QuantRatio getConsensusRatio(QuantCondition quantConditionNumerator,
			QuantCondition quantConditionDenominator, String replicateName) {
		return CensusRatio.NAN_RATIO;
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
		throw new IllegalArgumentException("getSTDRatio not available for " + getClass().getSimpleName());
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
	public String getDescription() {
		if (description == null) {
			description = protein.getDescription();
		}
		return description;
	}

	@Override
	public Integer getLength() {
		return protein.getLength();
	}

	/*
	 * (non-Javadoc)
	 * @see edu.scripps.yates.census.read.model.interfaces.QuantifiedItem#
	 * isQuantified()
	 */
	@Override
	public boolean isQuantified() {
		return false;
	}

}
