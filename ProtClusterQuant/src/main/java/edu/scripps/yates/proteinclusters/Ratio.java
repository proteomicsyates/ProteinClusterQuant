package edu.scripps.yates.proteinclusters;

import java.util.HashMap;
import java.util.Map;

import edu.scripps.yates.census.analysis.QuantCondition;

public class Ratio {
	private final Map<QuantCondition, Double> ionCountMap = new HashMap<QuantCondition, Double>();

	public Ratio() {

	}

	public void addIonCount(QuantCondition cond, double peakCount) {
		if (ionCountMap.containsKey(cond)) {
			ionCountMap.put(cond, ionCountMap.get(cond) + peakCount);
		} else {
			ionCountMap.put(cond, peakCount);
		}
	}

	public Double getCountRatio(QuantCondition cond1, QuantCondition cond2) {
		if (ionCountMap.containsKey(cond1) && ionCountMap.containsKey(cond2)) {
			return ionCountMap.get(cond1) / ionCountMap.get(cond2);
		}
		return null;
	}

	public Double getLog2CountRatio(QuantCondition cond1, QuantCondition cond2) {
		final Double countRatio = getCountRatio(cond1, cond2);
		if (countRatio == null)
			return null;
		return Math.log(countRatio) / Math.log(2.0);
	}

	public double getIonCount(QuantCondition cond) {
		return ionCountMap.get(cond);
	}
}
