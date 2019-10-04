package edu.scripps.yates.pcq.quantsite;

public class Mutation {
	private final String geneName;
	private final String enst;
	private final String cellLine;
	private final String mutation;
	private final String mutationType;
	private final String mutationType2;

	public Mutation(String line) {
		final String[] split = line.split("\t");
		geneName = split[0];
		enst = split[1];
		cellLine = split[2];
		mutation = split[3];
		mutationType = split[4];
		if (split.length == 6) {
			mutationType2 = split[5];
		} else {
			mutationType2 = null;
		}
	}

	public String getGeneName() {
		return geneName;
	}

	public String getEnst() {
		return enst;
	}

	public String getCellLine() {
		return cellLine;
	}

	public String getMutation() {
		return mutation;
	}

	public String getMutationType() {
		return mutationType;
	}

	public String getMutationType2() {
		return mutationType2;
	}
}
