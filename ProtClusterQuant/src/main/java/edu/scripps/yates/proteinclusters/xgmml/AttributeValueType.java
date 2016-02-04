package edu.scripps.yates.proteinclusters.xgmml;

public class AttributeValueType {
	private final AttType type;
	private final Object value;

	public AttributeValueType(Object value, AttType type) {
		this.value = value;
		this.type = type;
	}

	public static void main(String[] ar) {
		final Double valueOf = Double.valueOf("-Infinity");
		System.out.println(Double.POSITIVE_INFINITY);
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
			double d = Double.parseDouble(obj.toString());
		} catch (NumberFormatException nfe) {
			return false;
		}
		return true;
	}

	private boolean isInteger(Object obj) {
		try {
			int d = Integer.parseInt(obj.toString());
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
