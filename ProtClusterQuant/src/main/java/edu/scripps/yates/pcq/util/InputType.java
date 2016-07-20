package edu.scripps.yates.pcq.util;

public enum InputType {

	CENSUS_CHRO, CENSUS_OUT, SEPARATED_VALUES;

	public static String getPossibleValues() {
		StringBuilder sb = new StringBuilder();
		for (InputType inputType : InputType.values()) {
			if (!"".equals(sb.toString())) {
				sb.append(", ");
			}
			sb.append(inputType.name());
		}
		return sb.toString();
	}
}
