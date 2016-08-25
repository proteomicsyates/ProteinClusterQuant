package edu.scripps.yates.pcq.util;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.scripps.yates.pcq.compare.ComparisonInput;
import edu.scripps.yates.pcq.filter.PCQFilter;
import edu.scripps.yates.pcq.filter.PCQFilterByIonCount;
import edu.scripps.yates.pcq.filter.PCQFilterByPSMCount;

public class ProteinClusterQuantParameters {

	private boolean labelSwap;
	private double thresholdForSignificance;
	private boolean printOnlyFirstGene;
	private boolean ionsPerPeptideNodeThresholdOn;
	private int ionsPerPeptideNodeThreshold;
	private boolean psmsPerPeptideThresholdOn;
	private int psmsPerPeptideThreshold;
	private double iglewiczHoaglinTestThreshold;
	private boolean collapseIndistinguishableProteins;
	private boolean collapseIndistinguishablePeptides;
	private boolean makeAlignments;
	private boolean printKMeans;
	private char[] enzymeArray;
	private int missedCleavages;
	private boolean uniquePepOnly;
	private File uniprotReleasesFolder;
	private File inputFileFolder;
	private File outputFileFolder;
	private final List<ExperimentFiles> inputFiles = new ArrayList<ExperimentFiles>();
	private String outputPrefix;
	private String outputSuffix;
	private String lightSpecies;
	private String heavySpecies;
	private ColorManager colorManager;
	private File fastaFile;
	private String decoyRegexp;
	private int finalAlignmentScore;
	private double sequenceIdentity;
	private int minConsecutiveIdenticalAlignment;
	private File temporalOutputFolder;
	private static ProteinClusterQuantParameters instance;
	private ProteinNodeLabel proteinLabel;
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
	private String mongoDBURI;
	private String mongoProtDBName;
	private String mongoSeqDBName;
	private String mongoMassDBName;
	private boolean ignoreNotFoundPeptidesInDB;
	private AnalysisInputType inputType;
	private Double outliersRemovalFDR;
	private Double significantFDRThreshold;
	private boolean performRatioIntegration;
	private List<PCQFilter> filters;
	private String uniprotVersion;
	private boolean onlyOneSpectrumPerChromatographicPeakAndPerSaltStep;
	private boolean skipSingletons;
	private Color colorNonRegulated;
	private boolean generateMiscellaneousFiles;
	private String separator;
	private boolean applyClassificationsByProteinPair;
	private boolean analysisRun;
	private boolean comparisonRun;
	private ComparisonInput comparisonInput;

	private ProteinClusterQuantParameters() {

	}

	public static ProteinClusterQuantParameters getInstance() {
		if (instance == null) {
			instance = new ProteinClusterQuantParameters();
		}
		return instance;
	}

	/**
	 * @return the labelSwap
	 */
	public boolean isLabelSwap() {
		return labelSwap;
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
	 * @return the psmPerPeptideNodeThresholdOn
	 */
	public boolean isIonsPerPeptideNodeThresholdOn() {
		return ionsPerPeptideNodeThresholdOn;
	}

	/**
	 * @return the ionsPerPeptideThreshold
	 */
	public int getIonsPerPeptideNodeThreshold() {
		return ionsPerPeptideNodeThreshold;
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
	public List<ExperimentFiles> getInputFileNames() {
		return inputFiles;
	}

	public List<String> getExperimentNames() {
		List<String> ret = new ArrayList<String>();
		final List<ExperimentFiles> inputFileNames = getInputFileNames();
		for (ExperimentFiles experimentFiles : inputFileNames) {
			ret.add(experimentFiles.getExperimentName());
		}
		return ret;
	}

	public Map<String, List<String>> getReplicateNamesByExperimentNameMap() {
		Map<String, List<String>> map = new HashMap<String, List<String>>();
		for (ExperimentFiles experimentFiles : getInputFileNames()) {
			List<String> replicateNames = new ArrayList<String>();
			replicateNames.addAll(experimentFiles.getRelicateFileNames());
			map.put(experimentFiles.getExperimentName(), replicateNames);
		}
		return map;
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
	 * @return the colorManager
	 */
	public ColorManager getColorManager() {
		return colorManager;
	}

	public void setLabelSwap(boolean labelSwap) {
		this.labelSwap = labelSwap;

	}

	public void setThresholdForSignificance(double thresholdForSignificance) {
		this.thresholdForSignificance = thresholdForSignificance;
	}

	public void setPrintOnlyFirstGene(boolean printOnlyFirstGene) {
		this.printOnlyFirstGene = printOnlyFirstGene;
	}

	public void setIonsPerPeptideNodeThreshold(int ionsPerPeptideNodeThreshold) {
		this.ionsPerPeptideNodeThreshold = ionsPerPeptideNodeThreshold;
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

	public void setUniprotReleasesFolder(File uniprotReleasesFolder) {
		this.uniprotReleasesFolder = uniprotReleasesFolder;
	}

	public void setInputFileFolder(File inputFileFolder) {
		this.inputFileFolder = inputFileFolder;
	}

	public void setOutputFileFolder(File outputFileFolder) {

		this.outputFileFolder = outputFileFolder;
	}

	public void addInputFileNames(ExperimentFiles experimentFiles) {
		inputFiles.add(experimentFiles);
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
	 * @param ionsPerPeptideNodeThresholdOn
	 *            the ionsPerPeptideThresholdOn to set
	 */
	public void setIonsPerPeptideNodeThresholdOn(boolean ionsPerPeptideNodeThresholdOn) {
		this.ionsPerPeptideNodeThresholdOn = ionsPerPeptideNodeThresholdOn;
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

	public int getMinConsecutiveIdenticalAlignment() {
		return minConsecutiveIdenticalAlignment;
	}

	/**
	 * @param minConsecutiveIdenticalAlignment
	 *            the maxConsecutiveIdenticalAlignment to set
	 */
	public void setMinConsecutiveIdenticalAlignment(int minConsecutiveIdenticalAlignment) {
		this.minConsecutiveIdenticalAlignment = minConsecutiveIdenticalAlignment;
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

	/**
	 * @return the proteinLabel
	 */
	public ProteinNodeLabel getProteinLabel() {
		return proteinLabel;
	}

	/**
	 * @param proteinLabel
	 *            the proteinLabel to set
	 */
	public void setProteinLabel(ProteinNodeLabel proteinLabel) {
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

	/**
	 * @return the mongoDBURI
	 */
	public String getMongoDBURI() {
		return mongoDBURI;
	}

	/**
	 * @param mongoDBURI
	 *            the mongoDBURI to set
	 */
	public void setMongoDBURI(String mongoDBURI) {
		this.mongoDBURI = mongoDBURI;
	}

	/**
	 * @return the mongoProtDBName
	 */
	public String getMongoProtDBName() {
		return mongoProtDBName;
	}

	/**
	 * @param mongoProtDBName
	 *            the mongoProtDBName to set
	 */
	public void setMongoProtDBName(String mongoProtDBName) {
		this.mongoProtDBName = mongoProtDBName;
	}

	/**
	 * @return the mongoSeqDBName
	 */
	public String getMongoSeqDBName() {
		return mongoSeqDBName;
	}

	/**
	 * @param mongoSeqDBName
	 *            the mongoSeqDBName to set
	 */
	public void setMongoSeqDBName(String mongoSeqDBName) {
		this.mongoSeqDBName = mongoSeqDBName;
	}

	/**
	 * @return the mongoMassDBName
	 */
	public String getMongoMassDBName() {
		return mongoMassDBName;
	}

	/**
	 * @param mongoMassDBName
	 *            the mongoMassDBName to set
	 */
	public void setMongoMassDBName(String mongoMassDBName) {
		this.mongoMassDBName = mongoMassDBName;
	}

	public boolean isIgnoreNotFoundPeptidesInDB() {
		return ignoreNotFoundPeptidesInDB;
	}

	/**
	 * @param ignoreNotFoundPeptidesInDB
	 *            the ignoreNotFoundPeptidesInDB to set
	 */
	public void setIgnoreNotFoundPeptidesInDB(boolean ignoreNotFoundPeptidesInDB) {
		this.ignoreNotFoundPeptidesInDB = ignoreNotFoundPeptidesInDB;
	}

	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "ProteinClusterQuantParameters [applyClassificationsByProteinPair=" + applyClassificationsByProteinPair
				+ ", labelSwap=" + labelSwap + ", thresholdForSignificance=" + thresholdForSignificance
				+ ", printOnlyFirstGene=" + printOnlyFirstGene + ", ionsPerPeptideNodeThresholdOn="
				+ ionsPerPeptideNodeThresholdOn + ", ionsPerPeptideThreshold=" + ionsPerPeptideNodeThreshold
				+ ", psmsPerPeptideThresholdOn=" + psmsPerPeptideThresholdOn + ", psmsPerPeptideThreshold="
				+ psmsPerPeptideThreshold + ", iglewiczHoaglinTestThreshold=" + iglewiczHoaglinTestThreshold
				+ ", collapseIndistinguishableProteins=" + collapseIndistinguishableProteins
				+ ", collapseIndistinguishablePeptides=" + collapseIndistinguishablePeptides + ", makeAlignments="
				+ makeAlignments + ", printKMeans=" + printKMeans + ", enzymeArray=" + Arrays.toString(enzymeArray)
				+ ", missedCleavages=" + missedCleavages + ", uniquePepOnly=" + uniquePepOnly
				+ ", uniprotReleasesFolder=" + uniprotReleasesFolder + ", uniprotVersion=" + uniprotVersion
				+ ", inputFileFolder=" + inputFileFolder + ", outputFileFolder=" + outputFileFolder + ", inputFiles="
				+ inputFiles + ", outputPrefix=" + outputPrefix + ", outputSuffix=" + outputSuffix + ", lightSpecies="
				+ lightSpecies + ", heavySpecies=" + heavySpecies + ", colorManager=" + colorManager + ", fastaFile="
				+ fastaFile + ", decoyRegexp=" + decoyRegexp + ", finalAlignmentScore=" + finalAlignmentScore
				+ ", sequenceIdentity=" + sequenceIdentity + ", maxConsecutiveIdenticalAlignment="
				+ minConsecutiveIdenticalAlignment + ", temporalOutputFolder=" + temporalOutputFolder
				+ ", proteinLabel=" + proteinLabel + ", proteinNodeWidth=" + proteinNodeWidth + ", proteinNodeHeight="
				+ proteinNodeHeight + ", peptideNodeWidth=" + peptideNodeWidth + ", peptideNodeHeight="
				+ peptideNodeHeight + ", proteinNodeShape=" + proteinNodeShape + ", peptideNodeShape="
				+ peptideNodeShape + ", minimumRatioForColor=" + minimumRatioForColor + ", maximumRatioForColor="
				+ maximumRatioForColor + ", showCasesInEdges=" + showCasesInEdges + ", colorRatioMin=" + colorRatioMin
				+ ", colorRatioMax=" + colorRatioMax + ", remarkSignificantPeptides=" + remarkSignificantPeptides
				+ ", mongoDBURI=" + mongoDBURI + ", mongoProtDBName=" + mongoProtDBName + ", mongoSeqDBName="
				+ mongoSeqDBName + ", mongoMassDBName=" + mongoMassDBName + ", ignoreNotFoundPeptidesInDB="
				+ ignoreNotFoundPeptidesInDB + ", inputType=" + inputType + ", outliersRemovalFDR=" + outliersRemovalFDR
				+ ", significantFDRThreshold=" + significantFDRThreshold + " ]";
	}

	public String[] getInputFileNamesArray() {
		int size = 0;
		for (ExperimentFiles experimentFileName : inputFiles) {
			size += experimentFileName.getRelicateFileNames().size();
		}
		String[] ret = new String[size];
		int index = 0;
		for (ExperimentFiles experimentFileName : inputFiles) {
			final List<String> list = experimentFileName.getRelicateFileNames();
			for (String fileName : list) {
				ret[index++] = fileName;
			}
		}
		return ret;
	}

	public AnalysisInputType getInputType() {
		return inputType;
	}

	/**
	 * @param inputType
	 *            the inputType to set
	 */
	public void setInputType(AnalysisInputType inputType) {
		this.inputType = inputType;
	}

	public List<PCQFilter> getFilters() {
		if (filters == null) {
			filters = new ArrayList<PCQFilter>();
			if (isIonsPerPeptideNodeThresholdOn()) {
				filters.add(new PCQFilterByIonCount(getIonsPerPeptideNodeThreshold()));
			}
			if (isPsmsPerPeptideThresholdOn()) {
				filters.add(new PCQFilterByPSMCount(getPsmsPerPeptideThreshold()));
			}
			// TODO add more filters when available
		}
		return filters;
	}

	public Double getOutliersRemovalFDR() {
		return outliersRemovalFDR;
	}

	public void setOutliersRemovalFDR(Double fdr) {
		outliersRemovalFDR = fdr;
	}

	/**
	 * @return the psmsPerPeptideThresholdOn
	 */
	public boolean isPsmsPerPeptideThresholdOn() {
		return psmsPerPeptideThresholdOn;
	}

	/**
	 * @param psmsPerPeptideThresholdOn
	 *            the psmsPerPeptideThresholdOn to set
	 */
	public void setPsmsPerPeptideThresholdOn(boolean psmsPerPeptideThresholdOn) {
		this.psmsPerPeptideThresholdOn = psmsPerPeptideThresholdOn;
	}

	/**
	 * @return the psmsPerPeptideThreshold
	 */
	public int getPsmsPerPeptideThreshold() {
		return psmsPerPeptideThreshold;
	}

	/**
	 * @param psmsPerPeptideThreshold
	 *            the psmsPerPeptideThreshold to set
	 */
	public void setPsmsPerPeptideThreshold(int psmsPerPeptideThreshold) {
		this.psmsPerPeptideThreshold = psmsPerPeptideThreshold;
	}

	public boolean isPerformRatioIntegration() {
		return performRatioIntegration;
	}

	/**
	 * @param performRatioIntegration
	 *            the performRatioIntegration to set
	 */
	public void setPerformRatioIntegration(boolean performRatioIntegration) {
		this.performRatioIntegration = performRatioIntegration;
	}

	/**
	 * @return the significantFDRThreshold
	 */
	public Double getSignificantFDRThreshold() {
		return significantFDRThreshold;
	}

	/**
	 * @param significantFDRThreshold
	 *            the significantFDRThreshold to set
	 */
	public void setSignificantFDRThreshold(Double significantFDRThreshold) {
		this.significantFDRThreshold = significantFDRThreshold;
	}

	public String getUniprotVersion() {
		return uniprotVersion;
	}

	/**
	 * @param uniprotVersion
	 *            the uniprotVersion to set
	 */
	public void setUniprotVersion(String uniprotVersion) {
		this.uniprotVersion = uniprotVersion;
	}

	/**
	 * @return the onlyOneSpectrumPerChromatographicPeakAndPerSaltStep
	 */
	public boolean isOnlyOneSpectrumPerChromatographicPeakAndPerSaltStep() {
		return onlyOneSpectrumPerChromatographicPeakAndPerSaltStep;
	}

	/**
	 * @param onlyOneSpectrumPerChromatographicPeakAndPerSaltStep
	 *            the onlyOneSpectrumPerChromatographicPeakAndPerSaltStep to set
	 */
	public void setOnlyOneSpectrumPerChromatographicPeakAndPerSaltStep(
			boolean onlyOneSpectrumPerChromatographicPeakAndPerSaltStep) {
		this.onlyOneSpectrumPerChromatographicPeakAndPerSaltStep = onlyOneSpectrumPerChromatographicPeakAndPerSaltStep;
	}

	public boolean isSkipSingletons() {
		return skipSingletons;
	}

	/**
	 * @param skipSingletons
	 *            the skipSingletons to set
	 */
	public void setSkipSingletons(boolean skipSingletons) {
		this.skipSingletons = skipSingletons;
	}

	public void setColorNonRegulated(Color colorNonRegulated) {
		this.colorNonRegulated = colorNonRegulated;

	}

	/**
	 * @return the colorNonRegulated
	 */
	public Color getColorNonRegulated() {
		return colorNonRegulated;
	}

	public boolean isGenerateMiscellaneousFiles() {
		return generateMiscellaneousFiles;
	}

	/**
	 * @param generateMiscellaneousFiles
	 *            the generateMiscellaneousFiles to set
	 */
	public void setGenerateMiscellaneousFiles(boolean generateMiscellaneousFiles) {
		this.generateMiscellaneousFiles = generateMiscellaneousFiles;
	}

	public String getSeparator() {
		return separator;
	}

	/**
	 * @param separator
	 *            the separator to set
	 */
	public void setSeparator(String separator) {
		this.separator = separator;
	}

	public boolean isApplyClassificationsByProteinPair() {
		return applyClassificationsByProteinPair;
	}

	/**
	 * @param applyClassificationsByProteinPair
	 *            the applyClassificationsByProteinPair to set
	 */
	public void setApplyClassificationsByProteinPair(boolean applyClassificationsByProteinPair) {
		this.applyClassificationsByProteinPair = applyClassificationsByProteinPair;
	}

	public boolean isAnalysisRun() {
		return analysisRun;
	}

	/**
	 * @param analysisRun
	 *            the isAnalysisRun to set
	 */
	public void setAnalysisRun(boolean analysisRun) {
		this.analysisRun = analysisRun;
	}

	public boolean isComparisonRun() {
		return comparisonRun;
	}

	/**
	 * @param comparisonRun
	 *            the comparisonRun to set
	 */
	private void setComparisonRun(boolean comparisonRun) {
		this.comparisonRun = comparisonRun;
	}

	public ComparisonInput getComparisonInput() {
		return comparisonInput;
	}

	/**
	 * @param comparisonInput
	 *            the comparisonInput to set
	 */
	public void setComparisonInput(ComparisonInput comparisonInput) {
		this.comparisonInput = comparisonInput;
		if (comparisonInput != null) {
			setComparisonRun(true);
		}
	}
}
