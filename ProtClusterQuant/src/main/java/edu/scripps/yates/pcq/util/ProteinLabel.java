package edu.scripps.yates.proteinclusters.util;

public enum ProteinLabel {
	ACC, ID;

	public static ProteinLabel getFrom(String property) {
		for (ProteinLabel label : ProteinLabel.values()) {
			if (label.name().equalsIgnoreCase(property)) {
				return label;
			}
		}
		return null;
	}
}
