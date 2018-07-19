package edu.scripps.yates.pcq.xgmml.util;

public enum ProteinNodeLabel {
	ACC, ID, GENE;

	public static ProteinNodeLabel getFrom(String property) {
		for (final ProteinNodeLabel label : ProteinNodeLabel.values()) {
			if (label.name().equals(property.trim())) {
				return label;
			}
		}
		return null;
	}
}
