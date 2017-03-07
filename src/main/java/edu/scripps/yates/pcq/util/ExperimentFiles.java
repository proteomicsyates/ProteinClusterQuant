package edu.scripps.yates.pcq.util;

import java.util.ArrayList;
import java.util.List;

public class ExperimentFiles {
	private final String experimentName;
	private final List<String> relicateFileNames = new ArrayList<String>();

	public ExperimentFiles(String experimentName) {
		this.experimentName = experimentName;
	}

	public void addReplicateFileName(String replicateFileName) {
		if (!relicateFileNames.contains(replicateFileName)) {
			relicateFileNames.add(replicateFileName);
		} else {
			throw new IllegalArgumentException("Replicate name " + replicateFileName + " is repeated");
		}
	}

	/**
	 * @return the experimentName
	 */
	public String getExperimentName() {
		return experimentName;
	}

	/**
	 * @return the relicateFileNames
	 */
	public List<String> getRelicateFileNames() {
		return relicateFileNames;
	}
}
