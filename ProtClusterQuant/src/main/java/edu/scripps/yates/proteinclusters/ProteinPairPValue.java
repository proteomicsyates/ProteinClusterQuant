package edu.scripps.yates.proteinclusters;

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
