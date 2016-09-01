package edu.scripps.yates.pcq.util;

import edu.scripps.yates.pcq.model.ProteinPair;

public class ProteinPairPValue {
	private final double pValue;
	private final ProteinPair proteinpair;

	public ProteinPairPValue(double pValue, ProteinPair proteinpair) {

		this.pValue = pValue;
		this.proteinpair = proteinpair;
	}

	public double getpValue() {
		return pValue;
	}

	public ProteinPair getProteinpair() {
		return proteinpair;
	}

}
