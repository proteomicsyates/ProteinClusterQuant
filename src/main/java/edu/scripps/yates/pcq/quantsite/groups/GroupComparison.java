package edu.scripps.yates.pcq.quantsite.groups;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import edu.scripps.yates.pcq.compare.model.MyTTest;
import edu.scripps.yates.pcq.compare.model.QuantifiedSite;
import edu.scripps.yates.pcq.compare.model.QuantifiedSiteSet;
import edu.scripps.yates.pcq.util.PCQUtils;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.maths.PValueCorrection;
import edu.scripps.yates.utilities.maths.PValueCorrectionResult;
import edu.scripps.yates.utilities.maths.PValueCorrectionType;
import edu.scripps.yates.utilities.maths.PValuesCollection;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import smile.math.Math;

public class GroupComparison {
	private final static Logger log = Logger.getLogger(GroupComparison.class);
	private final String comparisonID;
	private final String subsetString;
	private final List<String> subsetlist = new ArrayList<String>();
	private final String sampleListString;
	private final List<String> sampleList = new ArrayList<String>();
	private QuantifiedSiteSet quantSites;
	private final boolean subsetting;
	private final PValueCorrectionType pValueCorrectionMethod;
	private final Map<QuantifiedSite, TDoubleList> group1Ratios = new THashMap<QuantifiedSite, TDoubleList>();
	private final Map<QuantifiedSite, TDoubleList> group2Ratios = new THashMap<QuantifiedSite, TDoubleList>();
	private final double distributionAvg;
	private final double distributionSigma;
	private final int numberSigmas;
	private final boolean useMayorityRule;

	public GroupComparison(String comparisonID, String subsetString, String sampleListString, boolean subsetting,
			PValueCorrectionType pValueCorrectionMethod, double distributionAvg, double distributionSigma,
			int numberSigmas, boolean useMayorityRule) {
		this.comparisonID = comparisonID;
		this.subsetString = subsetString;
		if (subsetString.contains(",")) {
			for (final String subset : subsetString.split(",")) {
				subsetlist.add(subset);
			}
		} else {
			subsetlist.add(subsetString);
		}
		this.sampleListString = sampleListString;
		if (sampleListString.contains(",")) {
			for (final String sampleName : sampleListString.split(",")) {
				sampleList.add(sampleName);
			}
		} else {
			sampleList.add(sampleListString);
		}
		this.subsetting = subsetting;
		this.pValueCorrectionMethod = pValueCorrectionMethod;
		this.distributionAvg = distributionAvg;
		this.distributionSigma = distributionSigma;
		this.numberSigmas = numberSigmas;
		this.useMayorityRule = useMayorityRule;
	}

	public String getComparisonID() {
		return comparisonID;
	}

	public String getSubsetString() {
		return subsetString;
	}

	public List<String> getSubsetlist() {
		return subsetlist;
	}

	public String getSampleListString() {
		return sampleListString;
	}

	public List<String> getSampleList() {
		return sampleList;
	}

	public void setSites(QuantifiedSiteSet quantSites) {
		this.quantSites = quantSites;
	}

	private Iterator<QuantifiedSite> getQuantSiteIterator() {
		if (!subsetting) {
			return quantSites.iterator();
		} else {
			return new QuantSiteSubsetIterator(quantSites, subsetlist);
		}
	}

	/**
	 * It analyzes the quant sites and returns how many of them are significant
	 * after doing t-tests and p-value correction
	 * 
	 * @param qvalueThreshold
	 * @return
	 */
	public TObjectDoubleMap<QuantifiedSite> run(double qvalueThreshold) {
		final List<String> totalSamples = quantSites.getUniqueSampleNames();
		final Iterator<QuantifiedSite> iterator = getQuantSiteIterator();
		final TObjectDoubleMap<QuantifiedSite> pValuesByQuantSite = new TObjectDoubleHashMap<QuantifiedSite>();
		final List<QuantifiedSite> quantSitesCompared = new ArrayList<QuantifiedSite>();

		while (iterator.hasNext()) {

			final QuantifiedSite quantSite = iterator.next();
			if (quantSite.getNodeKey().equals("Q15029#K654")) {
				log.info("asdf");
			}
			quantSitesCompared.add(quantSite);
			group1Ratios.put(quantSite, new TDoubleArrayList());
			group2Ratios.put(quantSite, new TDoubleArrayList());

			for (int sampleIndex = 0; sampleIndex < totalSamples.size(); sampleIndex++) {
				final Double ratio = quantSite.getLog2Ratio(sampleIndex);
				if (ratio != null && !Double.isNaN(ratio)) {
					if (sampleList.contains(totalSamples.get(sampleIndex))) {
						group1Ratios.get(quantSite).add(ratio);
					} else {
						group2Ratios.get(quantSite).add(ratio);
					}
				}
			}
			final MyTTest ttest = performTTest(quantSite, group1Ratios.get(quantSite), group2Ratios.get(quantSite),
					distributionAvg, distributionSigma, useMayorityRule);
			pValuesByQuantSite.put(quantSite, ttest.getPValue());
		}

		// p-value correction
		final PValuesCollection<QuantifiedSite> pValueCollection = new PValuesCollection<QuantifiedSite>(
				pValuesByQuantSite);
		PValueCorrectionResult<QuantifiedSite> pAdjust = null;
		if (pValuesByQuantSite.size() > 0) {
			log.info("Adjusting " + pValuesByQuantSite.size() + " p-values using method '"
					+ pValueCorrectionMethod.name() + "'");
			pAdjust = PValueCorrection.pAdjust(pValueCollection, pValueCorrectionMethod);
		} else {
			pAdjust = new PValueCorrectionResult<QuantifiedSite>();
			pAdjust.setCorrectedPValues(pValueCollection);
			pAdjust.setOriginalPValues(pValueCollection);
		}

		return pAdjust.getCorrectedPValues().getPValues();

	}

	private MyTTest performTTest(QuantifiedSite quantSite, TDoubleList ratios1, TDoubleList ratios2,
			double distributionAvg, double distributionSigma, boolean useMayorityRule) {

		Pair<Double, Integer> averagePair1 = PCQUtils.averageTakingIntoAccountInfinitiesAndNans(ratios1,
				useMayorityRule);
		if (averagePair1 == null) {
			averagePair1 = new Pair<Double, Integer>(Double.NaN, 0);
		}
		final double mean1 = averagePair1.getFirstelement();
		Pair<Double, Integer> averagePair2 = PCQUtils.averageTakingIntoAccountInfinitiesAndNans(ratios2,
				useMayorityRule);
		if (averagePair2 == null) {
			averagePair2 = new Pair<Double, Integer>(Double.NaN, 0);
		}
		final double mean2 = averagePair2.getFirstelement();
		final double stdev1 = PCQUtils.stdevTakingIntoAccountInfinitiesAndNans(ratios1, useMayorityRule);
		final double var1 = Math.pow(stdev1, 2);
		final double stdev2 = PCQUtils.stdevTakingIntoAccountInfinitiesAndNans(ratios2, useMayorityRule);
		final double var2 = Math.pow(stdev2, 2);
		final int n1 = averagePair1.getSecondElement();
		final int n2 = averagePair2.getSecondElement();

		if (n1 < 2 || n2 < 2 || Double.isNaN(stdev1) || Double.isNaN(stdev2)) {
			return new MyTTest(Double.NaN, false);
		}
		if ((Double.isInfinite(mean1) && Double.isInfinite(mean2))) {
			if (mean1 == mean2) {
				return new MyTTest(1.0, false); // no significant
			} else {
				return new MyTTest(0.0, false); // significant
			}
		} else {
			if (Double.isInfinite(mean1) || Double.isInfinite(mean2)) {
				// only significant if the one that is FINITE is beyond n sigmas
				// of the distribution
				final double finite = Double.isFinite(mean1) ? mean1 : mean2;
				final double infinite = Double.isInfinite(mean1) ? mean1 : mean2;
				if (Double.POSITIVE_INFINITY == infinite) {
					if (finite < distributionAvg + numberSigmas * distributionSigma) {
						return new MyTTest(0.0, false);
					} else {
						return new MyTTest(1.0, false);
					}
				} else {
					// negative infinity
					if (finite > distributionAvg + numberSigmas * distributionSigma) {
						// significant
						return new MyTTest(0.0, false);
					} else {
						// no significant
						return new MyTTest(1.0, false);
					}
				}
			}
		}
		final edu.scripps.yates.utilities.maths.TTest pValue = edu.scripps.yates.utilities.maths.TTest.test(mean1, var1,
				n1, mean2, var2, n2, true);
		return new MyTTest(pValue, true);

	}

	public Double getGroup1AvgRatio(QuantifiedSite quantSite) {
		return Maths.mean(group1Ratios.get(quantSite));
	}

	public Double getGroup2AvgRatio(QuantifiedSite quantSite) {
		return Maths.mean(group2Ratios.get(quantSite));
	}

	public Double getGroup1Stdev(QuantifiedSite quantSite) {
		return Maths.stddev(group1Ratios.get(quantSite));
	}

	public Double getGroup2Stdev(QuantifiedSite quantSite) {
		return Maths.stddev(group2Ratios.get(quantSite));
	}

}
