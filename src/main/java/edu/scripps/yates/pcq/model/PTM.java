package edu.scripps.yates.pcq.model;

import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.masses.AssignMass;
import edu.scripps.yates.utilities.sequence.PTMInProtein;

public class PTM {
	private final static Logger log = Logger.getLogger(PTM.class);
	private final Double massShift;
	private final char[] aas;

	public PTM(char aa) {
		this(null, aa);
	}

	public PTM(String aa) {
		this(null, aa);
	}

	public PTM(double massShift) {
		this(massShift, "");
	}

	public PTM(Double massShift, char aa) {
		if (!AssignMass.containsMass(aa)) {
			throw new IllegalArgumentException("'" + aa + "' not recognized.");
		}
		this.massShift = massShift;
		aas = new char[1];
		aas[0] = aa;
	}

	public PTM(Double massShift, String aasString) {
		this.massShift = massShift;
		aas = new char[aasString.length()];
		int i = 0;
		for (final char aa : aasString.toCharArray()) {
			if (!AssignMass.containsMass(aa)) {
				throw new IllegalArgumentException("'" + aa + "' not recognized.");
			}
			aas[i++] = aa;
		}
	}

	public PTM(Double massShift, char[] aas) {
		this.massShift = massShift;
		this.aas = aas;
		for (final char aa : aas) {
			if (!AssignMass.containsMass(aa)) {
				throw new IllegalArgumentException("'" + aa + "' not recognized.");
			}
		}
	}

	public double getMassShift() {
		return massShift;
	}

	public char[] getAas() {
		return aas;
	}

	public boolean isEquivalent(PTMInProtein ptmInProtein) {
		if (isAAValid(ptmInProtein.getAa()) && isMassShiftValid(ptmInProtein.getDeltaMass())) {
			return true;
		}
		return false;
	}

	private boolean isMassShiftValid(Double deltaMass) {
		if (massShift == null) {
			return true;
		}
		if (deltaMass == null) {
			return false;
		}
		if (Double.compare(deltaMass, massShift) == 0) {
			return true;
		}
		return false;
	}

	private boolean isAAValid(char aa) {
		if (aas.length == 0) {
			return true;
		}
		for (final char aa2 : aas) {
			if (aa == aa2) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Parse a string containing a PTM encoded on it.<br>
	 * The string will be as:
	 * <ul>
	 * <li>+34.5@K, meaning, ptms with a mass shift of +34.5 at K aminoacid</li>
	 * <li>K, meaning, any ptm at K aminoacid</li>
	 * <li>+34.5, meaning, any aminoacid with mass shift of +34.5</li>
	 * </ul>
	 * 
	 * @param string
	 * @return
	 */
	public static PTM parseString(String string) {
		try {
			string = string.trim();
			if (!"".equals(string)) {
				if (string.contains("@")) {
					final double massShift = Double.valueOf(string.substring(0, string.indexOf("@")));
					final String aasString = string.substring(string.indexOf("@") + 1);
					final char[] aas = new char[aasString.length()];
					int i = 0;
					for (final char c : aasString.toCharArray()) {
						aas[i++] = c;
					}
					return new PTM(massShift, aas);
				} else {
					// try as mass shift
					try {
						final double massShift = Double.valueOf(string);
						return new PTM(massShift);
					} catch (final NumberFormatException e) {
						// then it is just aminoacids
						return new PTM(string);
					}
				}
			}
		} catch (final Exception e) {
			log.error("Malformed PTM string: '" + string + "'");
		}
		return null;
	}

	@Override
	public String toString() {
		final StringBuilder sb = new StringBuilder();
		if (massShift != null) {
			sb.append(massShift);
		} else {
			sb.append("Any PTM");
		}
		sb.append(" at ");
		if (aas.length > 0) {
			for (final char c : aas) {
				sb.append(c + " ");
			}
		} else {
			sb.append("any AA");
		}
		return sb.toString();
	}

	public static void main(String[] args) {
		PTM ptm = PTM.parseString("+34.5@K");
		System.out.println(ptm);
		ptm = PTM.parseString("+79.96@PST");
		System.out.println(ptm);
		ptm = PTM.parseString("K");
		System.out.println(ptm);
		ptm = PTM.parseString("+79.96");
		System.out.println(ptm);
	}
}
