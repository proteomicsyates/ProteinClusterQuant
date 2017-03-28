package edu.scripps.yates.pcq.compare;

import java.util.Set;

import edu.scripps.yates.utilities.venndata.MultipleValuesObject;

public class ContainsFDRAndRatioItem extends MultipleValuesObject {
	private final double fdr;
	private final double ratio;
	private final double variance;
	private int numPSMs;

	public ContainsFDRAndRatioItem(String name, String separator, double fdr, double ratio, double variance,
			int numPSMs) {
		super(name, separator);
		this.fdr = fdr;
		this.ratio = ratio;
		this.variance = variance;
		this.numPSMs = numPSMs;
	}

	@Override
	public Set<String> getKeys() {
		return super.getKeys();
	}

	/**
	 * @return the fdr
	 */
	public double getFdr() {
		return fdr;
	}

	/**
	 * @return the ratio
	 */
	public double getRatio() {
		return ratio;
	}

	/**
	 * @return the numPSMs
	 */
	public int getNumPSMs() {
		return numPSMs;
	}

	/**
	 * @param numPSMs
	 *            the numPSMs to set
	 */
	public void setNumPSMs(int numPSMs) {
		this.numPSMs = numPSMs;
	}

	/**
	 * @return the variance
	 */
	public double getVariance() {
		return variance;
	}

	@Override
	public String toString() {
		return "[" + getName() + "\t" + ratio + "\tFDR=" + fdr + "\tVar=" + variance + "\t#psms=" + numPSMs + "]";
	}
}
