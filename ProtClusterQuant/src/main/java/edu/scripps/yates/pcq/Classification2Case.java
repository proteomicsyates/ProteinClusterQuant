package edu.scripps.yates.proteinclusters;

public enum Classification2Case {

	sign_shared(1, "Shared peptide is significantly different and is not shared by other protein", true), //
	sign_shared_third_protein(2, "Shared peptide is significant and is shared by other protein", true), //
	sign_unique(3,
			"Unique peptide is significantly different, or different between the two unique peptides is larger than threshold",
			true), //
			sign_unique_both(4,
					"Both unique peptides are significantly different and shared peptide is sgnificantly different in between",
					true), //
					error(5, "Error. This should not happen.", false), //
					no_difference(6,
							"There is not a significant differences. Difference between unique peptides is less than threshold",
							false), //
							unclassified(7, "Not classified", false);
	private final String explanation;
	private final boolean inconsistence;
	private final int caseID;

	private Classification2Case(int caseID, String explanation, boolean inconsistenceCase) {
		this.caseID = caseID;
		this.explanation = explanation;
		inconsistence = inconsistenceCase;
	}

	/**
	 * @return the inconsistence
	 */
	public boolean isInconsistence() {
		return inconsistence;
	}

	/**
	 * @return the explanation
	 */
	public String getExplanation() {
		return explanation;
	}

	/**
	 * @return the caseID
	 */
	public int getCaseID() {
		return caseID;
	}

	public static Classification2Case getByCaseID(int caseID) {
		final Classification2Case[] values = values();
		for (Classification2Case classification2Case : values) {
			if (caseID == classification2Case.getCaseID()) {
				return classification2Case;
			}
		}
		return null;
	}
}