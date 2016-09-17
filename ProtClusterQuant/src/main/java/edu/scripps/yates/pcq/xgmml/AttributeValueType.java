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
		type = evaluateType(value);

	}

	private AttType evaluateType(Object value2) {
		if (isNumeric(value)) {
			if (isInteger(value)) {
				return AttType.integer;
			} else {
				return AttType.real;
			}
		} else {

			return AttType.string;

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
