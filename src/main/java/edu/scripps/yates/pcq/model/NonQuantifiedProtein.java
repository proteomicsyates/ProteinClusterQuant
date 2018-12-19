package edu.scripps.yates.pcq.model;

import java.util.Collections;
import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.CensusRatio;
import edu.scripps.yates.census.read.model.QuantifiedProtein;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.proteomicsmodel.Ratio;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AggregationLevel;

public class NonQuantifiedProtein extends QuantifiedProtein {
	private final Protein protein;

	public NonQuantifiedProtein(Protein protein, boolean ignoreTaxonomies) {
		super(protein.getAccession(), ignoreTaxonomies);
		mergeWithProtein(protein);
		this.protein = protein;
	}

	@Override
	public Set<QuantRatio> getQuantRatios() {
		return Collections.emptySet();
	}

	@Override
	public Set<Ratio> getRatios() {
		return Collections.emptySet();
	}

	@Override
	public QuantRatio getConsensusRatio(QuantCondition quantConditionNumerator,
			QuantCondition quantConditionDenominator) {
		return CensusRatio.getNaNRatio(quantConditionNumerator, quantConditionDenominator, AggregationLevel.PROTEIN,
				"RATIO");

	}

	@Override
	public Set<QuantRatio> getNonInfinityRatios() {
		return Collections.emptySet();
	}

	@Override
	public boolean addQuantRatio(QuantRatio ratio) {
		throw new IllegalArgumentException("addQuantRatio not available for " + getClass().getSimpleName());
	}

	@Override
	public boolean addRatio(Ratio ratio) {
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
	public boolean addAmount(Amount amount) {
		throw new IllegalArgumentException("addAmount not available for " + getClass().getSimpleName());

	}

	@Override
	public String getDescription() {
		if (super.getDescription() == null) {
			setDescription(protein.getDescription());
		}
		return super.getDescription();
	}

	@Override
	public Integer getLength() {
		if (super.getLength() == null) {
			setLength(protein.getLength());
		}
		return super.getLength();
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
