package edu.scripps.yates.pcq.compare.model;

import edu.scripps.yates.utilities.maths.TTest;

public class MyTTest {
	private final TTest ttestResult;
	private final boolean useForPVAlueCorrection;
	private Double pvalue;
	private Double correctedPValue;

	public MyTTest(TTest ttestResult, boolean useForPVAlueCorrection) {
		this.ttestResult = ttestResult;
		this.useForPVAlueCorrection = useForPVAlueCorrection;
		pvalue = null;
	}

	public MyTTest(double pvalue, boolean useForPVAlueCorrection) {
		ttestResult = null;
		this.useForPVAlueCorrection = useForPVAlueCorrection;
		this.pvalue = pvalue;
	}

	public TTest getTtestResult() {
		return ttestResult;
	}

	public boolean isUseForPValueCorrection() {
		return useForPVAlueCorrection;
	}

	public double getPValue() {
		if (ttestResult != null) {
			return ttestResult.pvalue;
		} else {
			return pvalue;
		}
	}

	public void setPValue(double pValue) {
		if (ttestResult != null) {
			ttestResult.pvalue = pValue;
		} else {
			pvalue = pValue;
		}
	}

	public Double getCorrectedPValue() {
		return correctedPValue;
	}

	public void setCorrectedPValue(Double correctedPValue) {
		this.correctedPValue = correctedPValue;
	}
}
