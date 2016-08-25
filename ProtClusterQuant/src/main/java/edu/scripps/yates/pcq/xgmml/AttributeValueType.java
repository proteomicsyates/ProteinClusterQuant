package edu.scripps.yates.pcq.xgmml;

public class AttributeValueType {
	private final AttType type;
	private final Object value;

	public AttributeValueType(Object value, AttType type) {
		this.value = value;
		this.type = type;
	}

	public AttributeValueType(Object value) {
		this.value = value;
		if (isNumeric(value)) {
			if (isInteger(value)) {
				type = AttType.integer;
			} else {
				type = AttType.real;
			}
		} else {
			type = AttType.string;
		}
	}

	private boolean isNumeric(Object obj) {
		try {
			Double.parseDouble(obj.toString());
		} catch (NumberFormatException nfe) {
			return false;
		}
		return true;
	}

	private boolean isInteger(Object obj) {
		try {
			Integer.parseInt(obj.toString());
		} catch (NumberFormatException nfe) {
			return false;
		}
		return true;
	}

	/**
	 * @return the type
	 */
	public AttType getType() {
		return type;
	}

	/**
	 * @return the value
	 */
	public Object getValue() {
		return value;
	}
}
