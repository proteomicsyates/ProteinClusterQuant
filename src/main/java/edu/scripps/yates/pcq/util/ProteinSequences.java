package edu.scripps.yates.pcq.util;

import java.util.HashMap;

/**
 * Class wrapper of a map of protein sequences, that will be needed to check
 * whether
 * 
 * @author salvador
 *
 */
public class ProteinSequences extends HashMap<String, String> {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	@Override
	public String put(String key, String value) {
		if (key.contains("P04434")) {
			System.out.println(key + "\t" + value);
		}
		return super.put(key, value);
	}

}
