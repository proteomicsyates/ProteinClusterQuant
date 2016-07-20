package edu.scripps.yates.pcq;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.utilities.model.enums.AggregationLevel;
import edu.scripps.yates.utilities.model.enums.CombinationType;
import edu.scripps.yates.utilities.proteomicsmodel.Condition;
import edu.scripps.yates.utilities.proteomicsmodel.Score;

public class IonCountRatio implements QuantRatio {
	private final Map<QuantCondition, Double> ionCountMap = new HashMap<QuantCondition, Double>();
	private Score score;
	private CombinationType combinationType;
	private final AggregationLevel aggregationLevel;

	public IonCountRatio(AggregationLevel aggregationLevel) {
		this.aggregationLevel = aggregationLevel;
	}

	public void addIonCount(QuantCondition cond, double peakCount) {
		if (ionCountMap.containsKey(cond)) {
			ionCountMap.put(cond, ionCountMap.get(cond) + peakCount);
		} else {
			ionCountMap.put(cond, peakCount);
		}
	}

	@Override
	public Double getNonLogRatio(QuantCondition cond1, QuantCondition cond2) {
		if (ionCountMap.containsKey(cond1) && ionCountMap.containsKey(cond2)) {
			return ionCountMap.get(cond1) / ionCountMap.get(cond2);
		}
		return null;
	}

	@Override
	public Double getLog2Ratio(QuantCondition cond1, QuantCondition cond2) {
		final Double countRatio = getNonLogRatio(cond1, cond2);
		if (countRatio == null)
			return null;
		return Math.log(countRatio) / Math.log(2.0);
	}

	public Double getIonCount(QuantCondition cond) {
		return ionCountMap.get(cond);
	}

	@Override
	public double getValue() {
		throw new IllegalArgumentException("Not supported in a consensus isobaric count ratio");
	}

	@Override
	public Condition getCondition1() {
		List<QuantCondition> conditionList = getSortedConditionlist();
		if (conditionList.size() > 0)
			return conditionList.get(0);
		return null;
	}

	private List<QuantCondition> getSortedConditionlist() {
		final Set<QuantCondition> keySet = ionCountMap.keySet();
		List<QuantCondition> conditionList = new ArrayList<QuantCondition>();
		conditionList.addAll(keySet);
		Collections.sort(conditionList, new Comparator<QuantCondition>() {

			@Override
			public int compare(QuantCondition o1, QuantCondition o2) {
				return o1.getName().compareTo(o2.getName());
			}
		});
		return conditionList;
	}

	@Override
	public Condition getCondition2() {
		List<QuantCondition> conditionList = getSortedConditionlist();
		if (conditionList.size() > 1)
			return conditionList.get(1);
		return null;
	}

	@Override
	public String getDescription() {
		return "consensus ion count ratio";
	}

	@Override
	public Score getAssociatedConfidenceScore() {
		return score;
	}

	public void setAssociatedConfidenceScore(Score score) {
		this.score = score;
	}

	@Override
	public CombinationType getCombinationType() {
		return combinationType;
	}

	public void setCombinationType(CombinationType combinationType) {
		this.combinationType = combinationType;
	}

	@Override
	public AggregationLevel getAggregationLevel() {
		return aggregationLevel;
	}

	@Override
	public Double getLog2Ratio(QuantificationLabel labelNumerator, QuantificationLabel labelDenominator) {
		throw new IllegalArgumentException(
				"This is not supported. This is a consensus ratio that may come from ions counts with different labels");
	}

	@Override
	public Double getNonLogRatio(QuantificationLabel labelNumerator, QuantificationLabel labelDenominator) {
		throw new IllegalArgumentException(
				"This is not supported. This is a consensus ratio that may come from ions counts with different labels");

	}

	@Override
	public QuantificationLabel getLabel1() {
		throw new IllegalArgumentException(
				"This is not supported. This is a consensus ratio that may come from ions counts with different labels");

	}

	@Override
	public QuantificationLabel getLabel2() {
		throw new IllegalArgumentException(
				"This is not supported. This is a consensus ratio that may come from ions counts with different labels");

	}

	@Override
	public Double getLog2Ratio(String condition1Name, String condition2Name) {
		QuantCondition condition1 = getConditionByName(condition1Name);
		QuantCondition condition2 = getConditionByName(condition1Name);
		if (condition1 != null && condition2 != null) {
			return getLog2Ratio(condition1, condition2);
		}
		return null;
	}

	private QuantCondition getConditionByName(String conditionName) {
		final Set<QuantCondition> conditionSet = ionCountMap.keySet();
		for (QuantCondition quantCondition : conditionSet) {
			if (quantCondition.getName().equalsIgnoreCase(conditionName)) {
				return quantCondition;
			}
		}
		return null;
	}

	@Override
	public Double getNonLogRatio(String condition1Name, String condition2Name) {
		QuantCondition condition1 = getConditionByName(condition1Name);
		QuantCondition condition2 = getConditionByName(condition1Name);
		if (condition1 != null && condition2 != null) {
			return getNonLogRatio(condition1, condition2);
		}
		return null;
	}

	@Override
	public QuantCondition getQuantCondition1() {
		return (QuantCondition) getCondition1();
	}

	@Override
	public QuantCondition getQuantCondition2() {
		return (QuantCondition) getCondition2();
	}
}
