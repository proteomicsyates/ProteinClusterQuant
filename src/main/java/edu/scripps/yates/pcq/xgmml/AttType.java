package edu.scripps.yates.pcq.xgmml;

public enum AttType {
	real(AttCyType.Double), integer(AttCyType.Integer), string(AttCyType.String);
	private final AttCyType associatedCyType;

	private AttType(AttCyType associatedCyType) {
		this.associatedCyType = associatedCyType;
	}

	public AttCyType getAssociatedCyType() {
		return associatedCyType;
	}
}
