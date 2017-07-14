package edu.scripps.yates.pcq.params;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import edu.scripps.yates.census.analysis.QuantParameters;
import edu.scripps.yates.pcq.compare.ComparisonInput;
import edu.scripps.yates.pcq.filter.PCQFilter;
import edu.scripps.yates.pcq.filter.PCQFilterByIonCount;
import edu.scripps.yates.pcq.filter.PCQFilterByPSMCount;
import edu.scripps.yates.pcq.filter.PCQFilterByReplicateCount;
import edu.scripps.yates.pcq.sanxot.SanxotRunner;
import edu.scripps.yates.pcq.util.AnalysisInputType;
import edu.scripps.yates.pcq.util.ExperimentFiles;
import edu.scripps.yates.pcq.xgmml.util.ColorManager;
import edu.scripps.yates.pcq.xgmml.util.ProteinNodeLabel;
import edu.scripps.yates.pcq.xgmml.util.Shape;
import edu.scripps.yates.pcq.xgmml.util.UniprotAnnotationColumn;
import gnu.trove.map.hash.THashMap;

public class ProteinClusterQuantParameters {

	private boolean labelSwap;
	private double thresholdForSignificance;
	private boolean printOnlyFirstGene;
	private boolean ionsPerPeptideNodeThresholdOn;
	private int ionsPerPeptideNodeThreshold;
	private boolean psmsPerPeptideNodeThresholdOn;
	private int psmsPerPeptideNodeThreshold;
	private double iglewiczHoaglinTestThreshold;
	private boolean collapseIndistinguishableProteins;
	private boolean collapseIndistinguishablePeptides;
	private boolean makeAlignments;
	private char[] enzymeArray;
	private int missedCleavages;
	private boolean uniquePepOnly;
	private File uniprotReleasesFolder;
	private File inputFileFolder;
	private File outputFileFolder;
	private final List<ExperimentFiles> inputQuantificationFiles = new ArrayList<ExperimentFiles>();
	private final List<ExperimentFiles> inputIdentificationFiles = new ArrayList<ExperimentFiles>();
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
	private Double significantFDRThreshold;
	private boolean performRatioIntegration;
	private List<PCQFilter> filters;
	private String uniprotVersion;
	private boolean onlyOneSpectrumPerChromatographicPeakAndPerSaltStep;
	private boolean skipSingletons;
	private Color colorNonRegulated;
	private String separator = "\t";
	private boolean applyClassificationsByProteinPair;
	private boolean analysisRun;
	private boolean comparisonRun;
	private ComparisonInput comparisonInput;
	private int replicatesPerPeptideNodeThreshold;
	private boolean replicatesPerPeptideNodeThresholdOn;
	private final QuantParameters quantParameters;
	private boolean removeFilteredNodes;
	private boolean statisticalTestForProteinPairApplied;
	private boolean ignorePTMs;
	private boolean semiCleavage;
	private String peptideFilterRegexp;
	private List<UniprotAnnotationColumn> uniprotAnnotationColumns = new ArrayList<UniprotAnnotationColumn>();
	private char[] aaQuantified;

	private ProteinClusterQuantParameters() {
		quantParameters = new QuantParameters();
		quantParameters.setTimeout(SanxotRunner.DEFAULT_TIMEOUT);
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
	public List<ExperimentFiles> getInputQuantificationFileNames() {
		return inputQuantificationFiles;
	}

	public List<String> getExperimentNames() {
		List<String> ret = new ArrayList<String>();
		final List<ExperimentFiles> inputFileNames = getInputQuantificationFileNames();
		for (ExperimentFiles experimentFiles : inputFileNames) {
			ret.add(experimentFiles.getExperimentName());
		}
		return ret;
	}

	public Map<String, List<String>> getReplicateNamesByExperimentNameMap() {
		Map<String, List<String>> map = new THashMap<String, List<String>>();
		for (ExperimentFiles experimentFiles : getInputQuantificationFileNames()) {
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

	public void addQuantificationInputFileNames(ExperimentFiles experimentFiles) {
		inputQuantificationFiles.add(experimentFiles);
	}

	public void addIdentificationInputFileNames(ExperimentFiles experimentFiles) {
		inputIdentificationFiles.add(experimentFiles);
	}

	public void setOutputPrefix(String outputPrefix) {
		if (outputPrefix != null)
			this.outputPrefix = outputPrefix.trim();
	}

	public void setOutputSuffix(String outputSuffix) {
		if (outputSuffix != null)
			this.outputSuffix = outputSuffix.trim();
	}

	public void setLightSpecies(String lightSpecies) {
		if (lightSpecies != null)
			this.lightSpecies = lightSpecies.trim();
	}

	public void setHeavySpecies(String heavySpecies) {
		if (heavySpecies != null)
			this.heavySpecies = heavySpecies.trim();
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
		if (decoyRegexp != null)
			this.decoyRegexp = decoyRegexp.trim();
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
	 * 
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
		if (mongoDBURI != null)
			this.mongoDBURI = mongoDBURI.trim();
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
		if (mongoProtDBName != null)
			this.mongoProtDBName = mongoProtDBName.trim();
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
		if (mongoSeqDBName != null)
			this.mongoSeqDBName = mongoSeqDBName.trim();
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
		if (mongoMassDBName != null)
			this.mongoMassDBName = mongoMassDBName.trim();
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
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "ProteinClusterQuantParameters [applyClassificationsByProteinPair=" + applyClassificationsByProteinPair
				+ ", labelSwap=" + labelSwap + ", thresholdForSignificance=" + thresholdForSignificance
				+ ", printOnlyFirstGene=" + printOnlyFirstGene + ", ionsPerPeptideNodeThresholdOn="
				+ ionsPerPeptideNodeThresholdOn + ", ionsPerPeptideThreshold=" + ionsPerPeptideNodeThreshold
				+ ", psmsPerPeptideThresholdOn=" + psmsPerPeptideNodeThresholdOn + ", psmsPerPeptideThreshold="
				+ psmsPerPeptideNodeThreshold + ", iglewiczHoaglinTestThreshold=" + iglewiczHoaglinTestThreshold
				+ ", collapseIndistinguishableProteins=" + collapseIndistinguishableProteins
				+ ", collapseIndistinguishablePeptides=" + collapseIndistinguishablePeptides + ", makeAlignments="
				+ makeAlignments + ", enzymeArray=" + Arrays.toString(enzymeArray) + ", missedCleavages="
				+ missedCleavages + ", uniquePepOnly=" + uniquePepOnly + ", uniprotReleasesFolder="
				+ uniprotReleasesFolder + ", uniprotVersion=" + uniprotVersion + ", inputFileFolder=" + inputFileFolder
				+ ", outputFileFolder=" + outputFileFolder + ", inputQuantificationFiles=" + inputQuantificationFiles
				+ ", inputIdentificationFiles=" + inputIdentificationFiles + ", outputPrefix=" + outputPrefix
				+ ", outputSuffix=" + outputSuffix + ", lightSpecies=" + lightSpecies + ", heavySpecies=" + heavySpecies
				+ ", colorManager=" + colorManager + ", fastaFile=" + fastaFile + ", decoyRegexp=" + decoyRegexp
				+ ", finalAlignmentScore=" + finalAlignmentScore + ", sequenceIdentity=" + sequenceIdentity
				+ ", maxConsecutiveIdenticalAlignment=" + minConsecutiveIdenticalAlignment + ", temporalOutputFolder="
				+ temporalOutputFolder + ", proteinLabel=" + proteinLabel + ", proteinNodeWidth=" + proteinNodeWidth
				+ ", proteinNodeHeight=" + proteinNodeHeight + ", peptideNodeWidth=" + peptideNodeWidth
				+ ", peptideNodeHeight=" + peptideNodeHeight + ", proteinNodeShape=" + proteinNodeShape
				+ ", peptideNodeShape=" + peptideNodeShape + ", minimumRatioForColor=" + minimumRatioForColor
				+ ", maximumRatioForColor=" + maximumRatioForColor + ", showCasesInEdges=" + showCasesInEdges
				+ ", colorRatioMin=" + colorRatioMin + ", colorRatioMax=" + colorRatioMax
				+ ", remarkSignificantPeptides=" + remarkSignificantPeptides + ", mongoDBURI=" + mongoDBURI
				+ ", mongoProtDBName=" + mongoProtDBName + ", mongoSeqDBName=" + mongoSeqDBName + ", mongoMassDBName="
				+ mongoMassDBName + ", ignoreNotFoundPeptidesInDB=" + ignoreNotFoundPeptidesInDB + ", inputType="
				+ inputType + ", outliersRemovalFDR=" + quantParameters.getOutlierRemovalFDR()
				+ ", significantFDRThreshold=" + significantFDRThreshold + ", removeFilteredNodes="
				+ removeFilteredNodes + " ]";
	}

	public String[] getQuantInputFileNamesArray() {
		int size = 0;
		for (ExperimentFiles experimentFileName : inputQuantificationFiles) {
			size += experimentFileName.getRelicateFileNames().size();
		}
		String[] ret = new String[size];
		int index = 0;
		for (ExperimentFiles experimentFileName : inputQuantificationFiles) {
			final List<String> list = experimentFileName.getRelicateFileNames();
			for (String fileName : list) {
				ret[index++] = fileName.trim();
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
			if (isPsmsPerPeptideNodeThresholdOn()) {
				filters.add(new PCQFilterByPSMCount(getPsmsPerPeptideNodeThreshold()));
			}
			if (isReplicatesPerPeptideNodeThresholdOn()) {
				filters.add(new PCQFilterByReplicateCount(getReplicatesPerPeptideNodeThreshold()));
			}
			// TODO add more filters when available
		}
		return filters;
	}

	public void setOutliersRemovalFDR(Double fdr) {
		quantParameters.setOutlierRemovalFDR(fdr);
	}

	/**
	 * @return the psmsPerPeptideNodeThresholdOn
	 */
	public boolean isPsmsPerPeptideNodeThresholdOn() {
		return psmsPerPeptideNodeThresholdOn;
	}

	/**
	 * @param psmsPerPeptideNodeThresholdOn
	 *            the psmsPerPeptideThresholdOn to set
	 */
	public void setPsmsPerPeptideNodeThresholdOn(boolean psmsPerPeptideNodeThresholdOn) {
		this.psmsPerPeptideNodeThresholdOn = psmsPerPeptideNodeThresholdOn;
	}

	/**
	 * @return the replicatesPerPeptideNodeThresholdOn
	 */
	public boolean isReplicatesPerPeptideNodeThresholdOn() {
		return replicatesPerPeptideNodeThresholdOn;
	}

	/**
	 * @param psmsPerPeptideNodeThresholdOn
	 *            the psmsPerPeptideThresholdOn to set
	 */
	public void setReplicatesPerPeptideNodeThresholdOn(boolean replicatesPerPeptideNodeThresholdOn) {
		this.replicatesPerPeptideNodeThresholdOn = replicatesPerPeptideNodeThresholdOn;
	}

	/**
	 * @return the psmsPerPeptideNodeThreshold
	 */
	public int getPsmsPerPeptideNodeThreshold() {
		return psmsPerPeptideNodeThreshold;
	}

	/**
	 * @param psmsPerPeptideNodeThreshold
	 *            the psmsPerPeptideThreshold to set
	 */
	public void setPsmsPerPeptideNodeThreshold(int psmsPerPeptideNodeThreshold) {
		this.psmsPerPeptideNodeThreshold = psmsPerPeptideNodeThreshold;
	}

	/**
	 * @return the replicatesPerPeptideNodeThreshold
	 */
	public int getReplicatesPerPeptideNodeThreshold() {
		return replicatesPerPeptideNodeThreshold;
	}

	/**
	 * @param replicatesPerPeptideNodeThreshold
	 *            the replicatesPerPeptideNodeThreshold to set
	 */
	public void setReplicatesPerPeptideNodeThreshold(int replicatesPerPeptideNodeThreshold) {
		this.replicatesPerPeptideNodeThreshold = replicatesPerPeptideNodeThreshold;
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
		if (uniprotVersion != null)
			this.uniprotVersion = uniprotVersion.trim();
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

	public String getSeparator() {
		return separator;
	}

	/**
	 * @param separator
	 *            the separator to set
	 */
	public void setSeparator(String separator) {
		if (separator != null)
			this.separator = separator.trim();
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

	/**
	 * Whether PCQ performs an analysis or just compares networks
	 *
	 * @return
	 */
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

	/**
	 * @return the quantParameters
	 */
	public QuantParameters getQuantParameters() {
		return quantParameters;
	}

	public void setSanXotPath(File scriptPaths) {
		quantParameters.setSanxotScriptsFolder(scriptPaths);
	}

	public boolean isRemoveFilteredNodes() {
		return removeFilteredNodes;
	}

	/**
	 * @param removeFilteredNodes
	 *            the removeFilteredNodes to set
	 */
	public void setRemoveFilteredNodes(boolean removeFilteredNodes) {
		this.removeFilteredNodes = removeFilteredNodes;
	}

	public boolean isStatisticalTestForProteinPairApplied() {
		return statisticalTestForProteinPairApplied;
	}

	/**
	 * @param statisticalTestForProteinPairApplied
	 *            the statisticalTestForProteinPairApplied to set
	 */
	public void setStatisticalTestForProteinPairApplied(boolean statisticalTestForProteinPairApplied) {
		this.statisticalTestForProteinPairApplied = statisticalTestForProteinPairApplied;
	}

	/**
	 * @return the ignorePTMs
	 */
	public boolean isIgnorePTMs() {
		return ignorePTMs;
	}

	/**
	 * @param ignorePTMs
	 *            the ignorePTMs to set
	 */
	public void setIgnorePTMs(boolean ignorePTMs) {
		this.ignorePTMs = ignorePTMs;
	}

	public String[] getIdentificationInputFileNamesArray() {
		int size = 0;
		for (ExperimentFiles experimentFileName : inputIdentificationFiles) {
			size += experimentFileName.getRelicateFileNames().size();
		}
		String[] ret = new String[size];
		int index = 0;
		for (ExperimentFiles experimentFileName : inputIdentificationFiles) {
			final List<String> list = experimentFileName.getRelicateFileNames();
			for (String fileName : list) {
				ret[index++] = fileName.trim();
			}
		}
		return ret;
	}

	public boolean isSemiCleavage() {
		return semiCleavage;
	}

	/**
	 * @param semiCleavage
	 *            the semiCleavage to set
	 */
	public void setSemiCleavage(boolean semiCleavage) {
		this.semiCleavage = semiCleavage;
	}

	public String getPeptideFilterRegexp() {
		return peptideFilterRegexp;
	}

	/**
	 * @param peptideFilterRegexp
	 *            the peptideFilterRegexp to set
	 */
	public void setPeptideFilterRegexp(String peptideFilterRegexp) {
		this.peptideFilterRegexp = peptideFilterRegexp;
	}

	public List<UniprotAnnotationColumn> getUniprotAnnotationColumns() {
		return uniprotAnnotationColumns;
	}

	public void addUniprotAnnotationColumn(UniprotAnnotationColumn uniprotAnnotationColumn) {
		this.uniprotAnnotationColumns.add(uniprotAnnotationColumn);
	}

	public char[] getAaQuantified() {
		return aaQuantified;
	}

	public void setAaQuantified(char[] aaQuantified) {
		this.aaQuantified = aaQuantified;
	}

	public boolean isCollapseBySites() {
		return this.aaQuantified != null && this.aaQuantified.length > 0;
	}
}
