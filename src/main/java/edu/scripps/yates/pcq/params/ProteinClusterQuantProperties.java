package edu.scripps.yates.pcq.params;

import java.util.Properties;

import org.apache.log4j.Logger;

public class ProteinClusterQuantProperties extends Properties {
	/**
	 *
	 */
	private static final long serialVersionUID = -3860309523271704035L;
	private final static Logger log = Logger.getLogger(ProteinClusterQuantProperties.class);

	public ProteinClusterQuantProperties() {
		super();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.util.Properties#getProperty(java.lang.String, java.lang.String)
	 */
	@Override
	public String getProperty(String key, String defaultValue) {
		final String property = super.getProperty(key, defaultValue);
		if (property == null || "".equals(property)) {
			log.warn("Parameter '" + key + "' is not found. Using default value: '" + defaultValue + "'");
			return defaultValue;
		} else {
			return property.trim();
		}
	}

	/**
	 * Gets a property by propertyKey. If mandatory=true, it will throw an exception
	 * if not found. If mandatory=false, it will return null if not found.
	 *
	 * @param propertyKey
	 * @param mandatory
	 * @return
	 */
	public String getProperty(String propertyKey, boolean mandatory) {
		String ret = getProperty(propertyKey);
		if ((ret == null || "".equals(ret)) && mandatory) {
			throw new IllegalArgumentException("Parameter '" + propertyKey + "' is not found");
		}
		if ("".equals(ret)) {
			ret = null;
		}
		return ret;
	}

	@Override
	public synchronized boolean containsKey(Object key) {
		if (super.containsKey(key)) {
			if ("".equals(getProperty(key.toString()).trim())) {
				return false;
			}
			return true;
		}
		return false;
	}

}
