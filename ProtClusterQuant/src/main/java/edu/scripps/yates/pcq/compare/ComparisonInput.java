package edu.scripps.yates.pcq.compare;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class ComparisonInput {
	private final List<String> names = new ArrayList<String>();
	private final List<File> files = new ArrayList<File>();
	private double threshold = 0.05; // default value
	private ComparisonType comparisonType;

	enum ComparisonType {
		PEPTIDE_NODES, PROTEIN_NODES
	};

	public ComparisonInput(File comparisonInputFile) throws IOException {
		FileReader fileReader = new FileReader(comparisonInputFile);
		BufferedReader br = null;
		try {
			br = new BufferedReader(fileReader);
			String line = null;
			while ((line = br.readLine()) != null) {
				line = line.trim();
				if (line.startsWith("#") || "".equals(line)) {
					continue;
				}

				final String[] split = line.split("=");
				String name = split[0];
				if (name.equals("threshold")) {
					threshold = Double.valueOf(split[1]);
					System.out.println("Comparison FDR threshold: " + threshold);
				} else if (name.equals("comparison_type")) {
					comparisonType = ComparisonType.valueOf(split[1]);
					System.out.println("Comparison type: " + comparisonType);
				} else {
					File file = new File(split[1]);
					addName(name);
					addFile(file);
					System.out.println("Dataset: " + name + "\t" + file.getAbsolutePath());
				}
			}
		} finally {
			if (br != null) {
				br.close();
			}
		}
	}

	public List<File> getFiles() {
		return files;
	}

	public void addFile(File file) throws FileNotFoundException {
		if (!file.exists()) {
			throw new FileNotFoundException(file.getAbsolutePath() + " not found");
		}
		files.add(file);
	}

	public File getFile(int index) {
		if (files.size() > index)
			return files.get(index);
		return null;
	}

	public List<String> getNames() {
		return names;
	}

	public void addName(String name) {
		names.add(name);
	}

	public String getName(int index) {
		if (names.size() > index)
			return names.get(index);
		return null;
	}

	public double getThreshold() {
		return threshold;
	}

	/**
	 * @param threshold
	 *            the threshold to set
	 */
	public void setThreshold(double threshold) {
		this.threshold = threshold;
	}

	/**
	 * @return the comparisonType
	 */
	public ComparisonType getComparisonType() {
		return comparisonType;
	}

	/**
	 * @param comparisonType
	 *            the comparisonType to set
	 */
	public void setComparisonType(ComparisonType comparisonType) {
		this.comparisonType = comparisonType;
	}
}
