package edu.scripps.yates.pcq;

public enum Classification1Case {
	CASE3(0, "Having ratios but not significantly different", false), //
	CASE1(1, "Mix of infinities. No ratios.", true), //
	CASE2(2, "Mix of ratios and infinities", true), //
	CASE3Sig(3, "Shared peptide is significantly different", true), //
	CASE4(4, "Non classified case", false), //
	CASE5(5, "Some missing value (not quantified)", false), //
	CASE3ExplainedMiddleSig(6, "Shared peptide is significantly different and it may be explained by a third protein",
			true), //
			CASE3WrongMeassurementMiddleSig(7,
					"Shared peptide is significantly different and it is not shared by other protein", true) //
					;
	private final String explanation;
	private final boolean inconsistence;
	private final int caseID;

	private Classification1Case(int caseID, String explanation, boolean inconsistenceCase) {
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

	public static Classification1Case getByCaseID(int caseID) {
		final Classification1Case[] values = values();
		for (Classification1Case classification1Case : values) {
			if (caseID == classification1Case.getCaseID()) {
				return classification1Case;
			}
		}
		return null;
	}

}
