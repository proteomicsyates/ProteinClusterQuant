package edu.scripps.yates.proteinclusters.util;

import java.awt.Color;
import java.io.File;
import java.util.Arrays;

public class ProteinClusterQuantParameters {

	private boolean significantProteinPairAnalysis;
	private boolean labelSwap;
	private boolean stdAsSignficanceCutoffOn;
	private double thresholdForSignificance;
	private boolean printOnlyFirstGene;
	private boolean ionsPerPeptideThresholdOn;
	private int ionsPerPeptideThreshold;
	private double iglewiczHoaglinTestThreshold;
	private boolean collapseIndistinguishableProteins;
	private boolean collapseIndistinguishablePeptides;
	private boolean makeAlignments;
	private boolean printKMeans;
	private char[] enzymeArray;
	private int missedCleavages;
	private boolean uniquePepOnly;
	private boolean getTax;
	private File uniprotReleasesFolder;
	private File inputFileFolder;
	private File outputFileFolder;
	private String[] inputFileNames;
	private String outputPrefix;
	private String outputSuffix;
	private String lightSpecies;
	private String heavySpecies;
	private String[] replicateIdentifiers;
	private ColorManager colorManager;
	private File fastaFile;
	private String decoyRegexp;
	private int finalAlignmentScore;
	private double sequenceIdentity;
	private int maxConsecutiveIdenticalAlignment;
	private File temporalOutputFolder;
	private boolean excludeUniquePeptides;
	private static ProteinClusterQuantParameters instance;
	private ProteinLabel proteinLabel;
	private int proteinNodeWidth;
	private int proteinNodeHeight;
	private int peptideNodeWidth;
	private int peptideNodeHeight;
	private Shape proteinNodeShape;
	private Shape peptideNodeShape;
	private double minimumRatioForColor;
	private double maximumRatioForColor;
	private boolean showCasesInEdges;
	private Color colorRatioMin;
	private Color colorRatioMax;
	private boolean remarkSignificantPeptides;

	private ProteinClusterQuantParameters() {

	}

	public static ProteinClusterQuantParameters getInstance() {
		if (instance == null) {
			instance = new ProteinClusterQuantParameters();
		}
		return instance;
	}

	/**
	 * @return the significantProteinPairAnalysis
	 */
	public boolean isSignificantProteinPairAnalysis() {
		return significantProteinPairAnalysis;
	}

	/**
	 * @return the labelSwap
	 */
	public boolean isLabelSwap() {
		return labelSwap;
	}

	/**
	 * @return the stdAsSignficanceCutoffOn
	 */
	public boolean isStdAsSignficanceCutoffOn() {
		return stdAsSignficanceCutoffOn;
	}

	/**
	 * @return the thresholdForSignificance
	 */
	public double getThresholdForSignificance() {
		return thresholdForSignificance;
	}

	/**
	 * @return the printOnlyFirstGene
	 */
	public boolean isPrintOnlyFirstGene() {
		return printOnlyFirstGene;
	}

	/**
	 * @return the psmPerPeptideThresholdOn
	 */
	public boolean isIonsPerPeptideThresholdOn() {
		return ionsPerPeptideThresholdOn;
	}

	/**
	 * @return the ionsPerPeptideThreshold
	 */
	public int getIonsPerPeptideThreshold() {
		return ionsPerPeptideThreshold;
	}

	/**
	 * @return the iglewiczHoaglinTestThreshold
	 */
	public double getIglewiczHoaglinTestThreshold() {
		return iglewiczHoaglinTestThreshold;
	}

	/**
	 * @return the collapseIndistinguishableProteins
	 */
	public boolean isCollapseIndistinguishableProteins() {
		return collapseIndistinguishableProteins;
	}

	/**
	 * @return the makeAlignments
	 */
	public boolean isMakeAlignments() {
		return makeAlignments;
	}

	/**
	 * @return the printKMeans
	 */
	public boolean isPrintKMeans() {
		return printKMeans;
	}

	/**
	 * @return the enzymeArray
	 */
	public char[] getEnzymeArray() {
		return enzymeArray;
	}

	/**
	 * @return the missedCleavages
	 */
	public int getMissedCleavages() {
		return missedCleavages;
	}

	/**
	 * @return the uniquePepOnly
	 */
	public boolean isUniquePepOnly() {
		return uniquePepOnly;
	}

	/**
	 * @return the getTax
	 */
	public boolean isGetTax() {
		return getTax;
	}

	/**
	 * @return the uniprotReleasesFolder
	 */
	public File getUniprotReleasesFolder() {
		return uniprotReleasesFolder;
	}

	/**
	 * @return the inputFileFolder
	 */
	public File getInputFileFolder() {
		return inputFileFolder;
	}

	/**
	 * @return the outputFileFolder
	 */
	public File getOutputFileFolder() {
		return outputFileFolder;
	}

	/**
	 * @return the inputFileNames
	 */
	public String[] getInputFileNames() {
		return inputFileNames;
	}

	/**
	 * @return the outputPrefix
	 */
	public String getOutputPrefix() {
		return outputPrefix;
	}

	/**
	 * @return the outputSuffix
	 */
	public String getOutputSuffix() {
		return outputSuffix;
	}

	/**
	 * @return the lightSpecies
	 */
	public String getLightSpecies() {
		return lightSpecies;
	}

	/**
	 * @return the heavySpecies
	 */
	public String getHeavySpecies() {
		return heavySpecies;
	}

	/**
	 * @return the replicateIdentifiers
	 */
	public String[] getReplicateIdentifiers() {
		return replicateIdentifiers;
	}

	/**
	 * @return the colorManager
	 */
	public ColorManager getColorManager() {
		return colorManager;
	}

	public void setSignificantProteinPairAnalysis(boolean significantProteinPairAnalysis) {
		this.significantProteinPairAnalysis = significantProteinPairAnalysis;
	}

	public void setLabelSwap(boolean labelSwap) {
		this.labelSwap = labelSwap;

	}

	public void setStdAsSignficanceCutoffOn(boolean stdAsSignficanceCutoffOn) {
		this.stdAsSignficanceCutoffOn = stdAsSignficanceCutoffOn;
	}

	public void setThresholdForSignificance(double thresholdForSignificance) {
		this.thresholdForSignificance = thresholdForSignificance;
	}

	public void setPrintOnlyFirstGene(boolean printOnlyFirstGene) {
		this.printOnlyFirstGene = printOnlyFirstGene;
	}

	public void setIonsPerPeptideThreshold(int ionsPerPeptideThreshold) {
		this.ionsPerPeptideThreshold = ionsPerPeptideThreshold;
	}

	public void setIglewiczHoaglinTestThreshold(double iglewiczHoaglinTestThreshold) {
		this.iglewiczHoaglinTestThreshold = iglewiczHoaglinTestThreshold;
	}

	public void setCollapseIndistinguishableProteins(boolean collapseIndistinguishableProteins) {
		this.collapseIndistinguishableProteins = collapseIndistinguishableProteins;
	}

	public void setMakeAlignments(boolean makeAlignments) {
		this.makeAlignments = makeAlignments;
	}

	public void setPrintKMeans(boolean printKMeans) {
		this.printKMeans = printKMeans;
	}

	public void setEnzymeArray(char[] enzymeArray) {
		this.enzymeArray = enzymeArray;
	}

	public void setMissedCleavages(int missedCleavages) {
		this.missedCleavages = missedCleavages;
	}

	public void setUniquePepOnly(boolean uniquePepOnly) {
		this.uniquePepOnly = uniquePepOnly;
	}

	public void setGetTax(boolean getTax) {
		this.getTax = getTax;
	}

	public void setUniprotReleasesFolder(File uniprotReleasesFolder) {
		this.uniprotReleasesFolder = uniprotReleasesFolder;
	}

	public void setInputFileFolder(File inputFileFolder) {
		this.inputFileFolder = inputFileFolder;
	}

	public void setOutputFileFolder(File outputFileFolder) {

		this.outputFileFolder = outputFileFolder;
	}

	public void setInputFileNames(String[] inputFileNames) {
		this.inputFileNames = inputFileNames;
	}

	public void setOutputPrefix(String outputPrefix) {
		this.outputPrefix = outputPrefix;
	}

	public void setOutputSuffix(String outputSuffix) {
		this.outputSuffix = outputSuffix;
	}

	public void setLightSpecies(String lightSpecies) {
		this.lightSpecies = lightSpecies;
	}

	public void setHeavySpecies(String heavySpecies) {
		this.heavySpecies = heavySpecies;
	}

	public void setReplicateIdentifiers(String[] replicateIdentifiers) {
		this.replicateIdentifiers = replicateIdentifiers;
	}

	public void setColorManager(ColorManager colorManager) {
		this.colorManager = colorManager;
	}

	public void setFastaFile(File fastaFile) {
		this.fastaFile = fastaFile;

	}

	/**
	 * @return the fastaFile
	 */
	public File getFastaFile() {
		return fastaFile;
	}

	public void setDecoyRegexp(String decoyRegexp) {
		this.decoyRegexp = decoyRegexp;
	}

	/**
	 * @return the decoyRegexp
	 */
	public String getDecoyRegexp() {
		return decoyRegexp;
	}

	/**
	 * @param ionsPerPeptideThresholdOn
	 *            the ionsPerPeptideThresholdOn to set
	 */
	public void setIonsPerPeptideThresholdOn(boolean ionsPerPeptideThresholdOn) {
		this.ionsPerPeptideThresholdOn = ionsPerPeptideThresholdOn;
	}

	/**
	 * @return the collapseIndistinguishablePeptides
	 */
	public boolean isCollapseIndistinguishablePeptides() {
		return collapseIndistinguishablePeptides;
	}

	/**
	 * @param collapseIndistinguishablePeptides
	 *            the collapseIndistinguishablePeptides to set
	 */
	public void setCollapseIndistinguishablePeptides(boolean collapseIndistinguishablePeptides) {
		this.collapseIndistinguishablePeptides = collapseIndistinguishablePeptides;
	}

	public int getFinalAlignmentScore() {
		return finalAlignmentScore;
	}

	/**
	 * @param finalAlignmentScore
	 *            the finalAlignmentScore to set
	 */
	public void setFinalAlignmentScore(int finalAlignmentScore) {
		this.finalAlignmentScore = finalAlignmentScore;
	}

	public double getSequenceIdentity() {
		return sequenceIdentity;
	}

	/**
	 * @param sequenceIdentity
	 *            the sequenceIdentity to set
	 */
	public void setSequenceIdentity(double sequenceIdentity) {
		this.sequenceIdentity = sequenceIdentity;
	}

	public int getMaxConsecutiveIdenticalAlignment() {
		return maxConsecutiveIdenticalAlignment;
	}

	/**
	 * @param maxConsecutiveIdenticalAlignment
	 *            the maxConsecutiveIdenticalAlignment to set
	 */
	public void setMaxConsecutiveIdenticalAlignment(int maxConsecutiveIdenticalAlignment) {
		this.maxConsecutiveIdenticalAlignment = maxConsecutiveIdenticalAlignment;
	}

	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */

	public void setTemporalOutputFolder(File temporalOutputFolder) {
		this.temporalOutputFolder = temporalOutputFolder;

	}

	/**
	 * @return the temporalOutputFolder
	 */
	public File getTemporalOutputFolder() {
		return temporalOutputFolder;
	}

	public boolean isExcludeUniquePeptides() {
		return excludeUniquePeptides;
	}

	public void setExcludeUniquePeptides(boolean excludeUniquePeptides) {
		this.excludeUniquePeptides = excludeUniquePeptides;
	}

	/**
	 * @return the proteinLabel
	 */
	public ProteinLabel getProteinLabel() {
		return proteinLabel;
	}

	/**
	 * @param proteinLabel
	 *            the proteinLabel to set
	 */
	public void setProteinLabel(ProteinLabel proteinLabel) {
		this.proteinLabel = proteinLabel;
	}

	/**
	 * @return the proteinNodeShape
	 */
	public Shape getProteinNodeShape() {
		return proteinNodeShape;
	}

	/**
	 * @param proteinNodeShape
	 *            the proteinNodeShape to set
	 */
	public void setProteinNodeShape(Shape proteinNodeShape) {
		this.proteinNodeShape = proteinNodeShape;
	}

	/**
	 * @return the peptideNodeShape
	 */
	public Shape getPeptideNodeShape() {
		return peptideNodeShape;
	}

	/**
	 * @param peptideNodeShape
	 *            the peptideNodeShape to set
	 */
	public void setPeptideNodeShape(Shape peptideNodeShape) {
		this.peptideNodeShape = peptideNodeShape;
	}

	/**
	 * @return the proteinNodeWidth
	 */
	public int getProteinNodeWidth() {
		return proteinNodeWidth;
	}

	/**
	 * @param proteinNodeWidth
	 *            the proteinNodeWidth to set
	 */
	public void setProteinNodeWidth(int proteinNodeWidth) {
		this.proteinNodeWidth = proteinNodeWidth;
	}

	/**
	 * @return the proteinNodeHeight
	 */
	public int getProteinNodeHeight() {
		return proteinNodeHeight;
	}

	/**
	 * @param proteinNodeHeight
	 *            the proteinNodeHeight to set
	 */
	public void setProteinNodeHeight(int proteinNodeHeight) {
		this.proteinNodeHeight = proteinNodeHeight;
	}

	/**
	 * @return the peptideNodeWidth
	 */
	public int getPeptideNodeWidth() {
		return peptideNodeWidth;
	}

	/**
	 * @param peptideNodeWidth
	 *            the peptideNodeWidth to set
	 */
	public void setPeptideNodeWidth(int peptideNodeWidth) {
		this.peptideNodeWidth = peptideNodeWidth;
	}

	/**
	 * @return the peptideNodeHeight
	 */
	public int getPeptideNodeHeight() {
		return peptideNodeHeight;
	}

	/**
	 * @param peptideNodeHeight
	 *            the peptideNodeHeight to set
	 */
	public void setPeptideNodeHeight(int peptideNodeHeight) {
		this.peptideNodeHeight = peptideNodeHeight;
	}

	/**
	 * @return the minimumRatioForColor
	 */
	public double getMinimumRatioForColor() {
		return minimumRatioForColor;
	}

	/**
	 * @param minimumRatioForColor
	 *            the minimumRatioForColor to set
	 */
	public void setMinimumRatioForColor(double minimumRatioForColor) {
		this.minimumRatioForColor = minimumRatioForColor;
	}

	/**
	 * @return the maximumRatioForColor
	 */
	public double getMaximumRatioForColor() {
		return maximumRatioForColor;
	}

	/**
	 * @param maximumRatioForColor
	 *            the maximumRatioForColor to set
	 */
	public void setMaximumRatioForColor(double maximumRatioForColor) {
		this.maximumRatioForColor = maximumRatioForColor;
	}

	/**
	 * @return the showCasesInEdges
	 */
	public boolean isShowCasesInEdges() {
		return showCasesInEdges;
	}

	/**
	 * @param showCasesInEdges
	 *            the showCasesInEdges to set
	 */
	public void setShowCasesInEdges(boolean showCasesInEdges) {
		this.showCasesInEdges = showCasesInEdges;
	}

	/**
	 * @return the colorRatioMax
	 */
	public Color getColorRatioMax() {
		return colorRatioMax;
	}

	/**
	 * @param colorRatioMax
	 *            the colorRatioMax to set
	 */
	public void setColorRatioMax(Color colorRatioMax) {
		this.colorRatioMax = colorRatioMax;
	}

	/**
	 * @return the colorRatioMin
	 */
	public Color getColorRatioMin() {
		return colorRatioMin;
	}

	/**
	 * @param colorRatioMin
	 *            the colorRatioMin to set
	 */
	public void setColorRatioMin(Color colorRatioMin) {
		this.colorRatioMin = colorRatioMin;
	}

	/**
	 * @return the remarkSignificantPeptides
	 */
	public boolean isRemarkSignificantPeptides() {
		return remarkSignificantPeptides;
	}

	/**
	 * @param remarkSignificantPeptides
	 *            the remarkSignificantPeptides to set
	 */
	public void setRemarkSignificantPeptides(boolean remarkSignificantPeptides) {
		this.remarkSignificantPeptides = remarkSignificantPeptides;
	}

	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "ProteinClusterQuantParameters [significantProteinPairAnalysis=" + significantProteinPairAnalysis
				+ ", labelSwap=" + labelSwap + ", stdAsSignficanceCutoffOn=" + stdAsSignficanceCutoffOn
				+ ", thresholdForSignificance=" + thresholdForSignificance + ", printOnlyFirstGene="
				+ printOnlyFirstGene + ", ionsPerPeptideThresholdOn=" + ionsPerPeptideThresholdOn
				+ ", ionsPerPeptideThreshold=" + ionsPerPeptideThreshold + ", iglewiczHoaglinTestThreshold="
				+ iglewiczHoaglinTestThreshold + ", collapseIndistinguishableProteins="
				+ collapseIndistinguishableProteins + ", collapseIndistinguishablePeptides="
				+ collapseIndistinguishablePeptides + ", makeAlignments=" + makeAlignments + ", printKMeans="
				+ printKMeans + ", enzymeArray=" + Arrays.toString(enzymeArray) + ", missedCleavages=" + missedCleavages
				+ ", uniquePepOnly=" + uniquePepOnly + ", getTax=" + getTax + ", uniprotReleasesFolder="
				+ uniprotReleasesFolder + ", inputFileFolder=" + inputFileFolder + ", outputFileFolder="
				+ outputFileFolder + ", inputFileNames=" + Arrays.toString(inputFileNames) + ", outputPrefix="
				+ outputPrefix + ", outputSuffix=" + outputSuffix + ", lightSpecies=" + lightSpecies + ", heavySpecies="
				+ heavySpecies + ", replicateIdentifiers=" + Arrays.toString(replicateIdentifiers) + ", colorManager="
				+ colorManager + ", fastaFile=" + fastaFile + ", decoyRegexp=" + decoyRegexp + ", finalAlignmentScore="
				+ finalAlignmentScore + ", sequenceIdentity=" + sequenceIdentity + ", maxConsecutiveIdenticalAlignment="
				+ maxConsecutiveIdenticalAlignment + ", temporalOutputFolder=" + temporalOutputFolder
				+ ", excludeUniquePeptides=" + excludeUniquePeptides + ", proteinLabel=" + proteinLabel
				+ ", proteinNodeWidth=" + proteinNodeWidth + ", proteinNodeHeight=" + proteinNodeHeight
				+ ", peptideNodeWidth=" + peptideNodeWidth + ", peptideNodeHeight=" + peptideNodeHeight
				+ ", proteinNodeShape=" + proteinNodeShape + ", peptideNodeShape=" + peptideNodeShape
				+ ", minimumRatioForColor=" + minimumRatioForColor + ", maximumRatioForColor=" + maximumRatioForColor
				+ ", showCasesInEdges=" + showCasesInEdges + ", colorRatioMin=" + colorRatioMin + ", colorRatioMax="
				+ colorRatioMax + ", remarkSignificantPeptides=" + remarkSignificantPeptides + "]";
	}

}
