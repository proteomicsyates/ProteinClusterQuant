package edu.scripps.yates.pcq.model;

import java.util.Collections;
import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.CensusRatio;
import edu.scripps.yates.census.read.model.QuantifiedProtein;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.proteomicsmodel.MSRun;
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

	@Override
	public void mergeWithProtein(Protein protein) {

		if (protein.getCoverage() != null) {
			setCoverage(protein.getCoverage());
		}
		if (getPrimaryAccession() == null) {
			setPrimaryAccession(protein.getPrimaryAccession());
		}
		final String description = protein.getDescription();
		if (protein.getDescription() != null) {
			setDescription(description);
		}
		if (protein.getEmpai() != null) {
			setEmpai(protein.getEmpai());
		}
		if (protein.getLength() != null) {
			setLength(protein.getLength());
		}

		if (getMw() == null) {
			setMw(protein.getMw());
		}
		if (protein.getNsaf() != null) {
			setNsaf(protein.getNsaf());
		}
		if (getNsaf_norm() == null) {
			setNsaf_norm(protein.getNsaf_norm());
		}
		if (getPi() == null) {
			setPi(protein.getPi());
		}
		// protein may have to be grouped again since it contains new PSMs, so
		// we reset group
		setProteinGroup(null);
		// spectrumCount may not be accurate now with new psms
		setSpectrumCount(null);

//		if (protein.getPSMs() != null) {
//			for (final PSM psm : protein.getPSMs()) {
//				addPSM(psm, true);
//			}
//		}
//		if (protein.getPeptides() != null) {
//			for (final Peptide peptide2 : protein.getPeptides()) {
//				boolean found = false;
//				for (final Peptide peptide : getPeptides()) {
//					if (peptide.equals(peptide2)) {
//						found = true;
//						peptide.mergeWithPeptide(peptide2);
//					}
//				}
//				if (!found) {
//					addPeptide(peptide2, true);
//				}
//			}
//		}

		if (getSearchEngine() == null) {
			setSearchEngine(protein.getSearchEngine());
		}
		if (protein.getMSRuns() != null) {
			for (final MSRun msRun : protein.getMSRuns()) {
				addMSRun(msRun);
			}
		}
	}
}
