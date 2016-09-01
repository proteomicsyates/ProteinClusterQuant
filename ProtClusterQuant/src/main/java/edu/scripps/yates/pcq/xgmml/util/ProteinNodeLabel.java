package edu.scripps.yates.pcq.xgmml.util;

public enum ProteinNodeLabel {
	ACC, ID, GENE;

	public static ProteinNodeLabel getFrom(String property) {
		for (ProteinNodeLabel label : ProteinNodeLabel.values()) {
			if (label.name().equalsIgnoreCase(property)) {
				return label;
			}
		}
		return null;
	}
}
