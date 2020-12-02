package edu.scripps.yates.pcq.params;

import java.awt.Color;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import org.apache.commons.io.FilenameUtils;

import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.pcq.model.IsobaricRatioType;
import edu.scripps.yates.pcq.util.AnalysisInputType;
import edu.scripps.yates.pcq.util.ExperimentFiles;
import edu.scripps.yates.pcq.xgmml.util.ColorManager;
import edu.scripps.yates.pcq.xgmml.util.ProteinNodeLabel;
import edu.scripps.yates.pcq.xgmml.util.Shape;
import edu.scripps.yates.pcq.xgmml.util.UniprotAnnotationColumn;
import edu.scripps.yates.utilities.colors.ColorGenerator;

public class PropertiesReader {
	private static final String PROPERTIES_FILE_NAME = "setup.properties";
	private final static org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(PropertiesReader.class);

	public static ProteinClusterQuantProperties readerProperties(boolean batchMode) throws IOException {
		final ClassLoader cl = PropertiesReader.class.getClassLoader();
		final InputStream inputStream = cl.getResourceAsStream(PROPERTIES_FILE_NAME);
		return readProperties(inputStream, batchMode);
	}

	public static ProteinClusterQuantProperties readProperties(File setupPropertiesFile, boolean batchMode)
			throws IOException {
		if (setupPropertiesFile == null || !setupPropertiesFile.exists()) {
			throw new IllegalArgumentException(
					"parameters file '" + setupPropertiesFile.getAbsolutePath() + "' is null or doesn't exist");
		}
		try {
			final InputStream inputStream = new FileInputStream(setupPropertiesFile);
			final ProteinClusterQuantProperties prop = readProperties(inputStream, batchMode);
			return prop;
		} catch (final IOException e) {
			throw e;
		} catch (final IllegalArgumentException e) {
			e.printStackTrace();
			throw new IllegalArgumentException("Error reading properties file '" + setupPropertiesFile.getAbsolutePath()
					+ "'. One common error is to write file paths with '\' (for windows). Use '\\' or '/' instead and try again");
		}
	}

	public static ProteinClusterQuantProperties readProperties(InputStream inputStream, boolean batchMode)
			throws IOException {

		if (inputStream == null) {
			throw new IllegalArgumentException("input stream is null");
		}

		final ProteinClusterQuantProperties prop = new ProteinClusterQuantProperties();
		prop.load(inputStream);
		readParametersFromProperties(prop, batchMode);
		return prop;

	}

	private static void readParametersFromProperties(ProteinClusterQuantProperties properties) throws IOException {
		readParametersFromProperties(properties, false);
	}

	private static void readParametersFromProperties(ProteinClusterQuantProperties properties, boolean batchMode)
			throws IOException {
		final ProteinClusterQuantParameters params = ProteinClusterQuantParameters.getInstance();

		if (properties.containsKey("inputType")) {
			try {
				final AnalysisInputType inputType = AnalysisInputType
						.valueOf(properties.getProperty("inputType", false));
				params.setAnalysisInputType(inputType);
			} catch (final NullPointerException e) {

			}
		}
		if (params.getAnalysisInputType() == null) {
			log.info("inputType not present");
		}

		if (properties.containsKey("quantChannels")) {
			try {
				final String inputType = properties.getProperty("quantChannels", false);
				if (inputType != null && !"".equals(inputType)) {
					if (inputType.contains("/")) {
						final String[] split = inputType.split("/");
						final QuantificationLabel numeratorLabel = QuantificationLabel.getByName(split[0]);
						if (numeratorLabel == null) {
							throw new IllegalArgumentException(
									"'quantChannels' label '" + split[0] + "' is not recognized. Possible values are "
											+ QuantificationLabel.getValuesString());
						}
						params.setNumeratorLabel(numeratorLabel);
						final QuantificationLabel denominatorLabel = QuantificationLabel.getByName(split[1]);
						if (denominatorLabel == null) {
							throw new IllegalArgumentException(
									"'quantChannels' label '" + split[1] + "' is not recognized. Possible values are "
											+ QuantificationLabel.getValuesString());
						}
						params.setDenominatorLabel(denominatorLabel);
					} else {
						throw new IllegalArgumentException(
								"'quantChannels' parameter is not well formed. Format LABEL1/LABEL2 where a LABEL can be "
										+ QuantificationLabel.getValuesString());
					}
				} else {
					params.setNumeratorLabel(QuantificationLabel.LIGHT);
					params.setDenominatorLabel(QuantificationLabel.HEAVY);
				}
			} catch (final Exception e) {
				throw new IllegalArgumentException("'quantChannels' parameter error");
			}
		} else {
			log.info("quantChannels not present. setting default values: L/H");
			params.setNumeratorLabel(QuantificationLabel.LIGHT);
			params.setDenominatorLabel(QuantificationLabel.HEAVY);
		}

		if (properties.containsKey("onlyOneSpectrumPerChromatographicPeakAndPerSaltStep")) {
			final boolean onlyOneSpectrumPerChromatographicPeakAndPerSaltStep = Boolean
					.valueOf(properties.getProperty("onlyOneSpectrumPerChromatographicPeakAndPerSaltStep", false));
			params.setOnlyOneSpectrumPerChromatographicPeakAndPerSaltStep(
					onlyOneSpectrumPerChromatographicPeakAndPerSaltStep);
		}
		if (properties.containsKey("skipSingletons")) {
			final boolean skipSingletons = Boolean.valueOf(properties.getProperty("skipSingletons", false));
			params.setSkipSingletons(skipSingletons);
		}

		final boolean applyClassificationsByProteinPair = Boolean
				.valueOf(properties.getProperty("applyClassificationsByProteinPairs", "false"));
		params.setApplyClassificationsByProteinPair(applyClassificationsByProteinPair);

		final boolean labelSwap = Boolean.valueOf(properties.getProperty("labelSwap", "false"));
		params.setLabelSwap(labelSwap);

		if (applyClassificationsByProteinPair) {
			final boolean statisticalTestForProteinPairApplied = Boolean
					.valueOf(properties.getProperty("statisticalTestForProteinPairApplied", "false"));
			params.setStatisticalTestForProteinPairApplied(statisticalTestForProteinPairApplied);
			// cutoff for significance
			final double thresholdForSignificance = Double
					.valueOf(properties.getProperty("thresholdForSignificance", "2"));
			params.setThresholdForSignificance(thresholdForSignificance);
			// threshold for the iglewiczHoaglin Test. A result wquals or greater
			// than that, would be considered as an outlier
			final double iglewiczHoaglinTestThreshold = Double
					.valueOf(properties.getProperty("iglewiczHoaglinTest", "3.5"));
			params.setIglewiczHoaglinTestThreshold(iglewiczHoaglinTestThreshold);
			// do we only count truly unique peptides as unique
			final boolean uniquePepOnly = Boolean.valueOf(properties.getProperty("uniquePepOnly", "true"));
			params.setUniquePepOnly(uniquePepOnly);
		}

		final boolean printOnlyFirstGene = Boolean.valueOf(properties.getProperty("printOnlyFirstGene", "true"));
		params.setPrintOnlyFirstGene(printOnlyFirstGene);
		// do we have threhsold for minimum ion counts per Peptide
		final boolean ionsPerPeptideNodeThresholdOn = Boolean
				.valueOf(properties.getProperty("ionsPerPeptideNodeThresholdOn", "false"));
		params.setIonsPerPeptideNodeThresholdOn(ionsPerPeptideNodeThresholdOn);
		// threshold for minimum ions per peptide node
		final int ionsPerPeptideNodeThreshold = Integer
				.valueOf(properties.getProperty("ionsPerPeptideNodeThreshold", "0"));
		params.setIonsPerPeptideNodeThreshold(ionsPerPeptideNodeThreshold);
		// threshold for minimum PSMs per Peptide
		final boolean psmsPerPeptideNodeThresholdOn = Boolean
				.valueOf(properties.getProperty("psmsPerPeptideNodeThresholdOn", "false"));
		params.setPsmsPerPeptideNodeThresholdOn(psmsPerPeptideNodeThresholdOn);
		// threshold for minimum PSMs per peptide node
		final int psmsPerPeptideNodeThreshold = Integer
				.valueOf(properties.getProperty("psmsPerPeptideNodeThreshold", "0"));
		params.setPsmsPerPeptideNodeThreshold(psmsPerPeptideNodeThreshold);
		// threshold for minimum replicates per Peptide Node
		final boolean replicatesPerPeptideNodeThresholdOn = Boolean
				.valueOf(properties.getProperty("replicatesPerPeptideNodeThresholdOn", "false"));
		params.setReplicatesPerPeptideNodeThresholdOn(replicatesPerPeptideNodeThresholdOn);
		// threshold for minimum replicates per peptide node
		final int replicatesPerPeptideNodeThreshold = Integer
				.valueOf(properties.getProperty("replicatesPerPeptideNodeThreshold", "0"));
		params.setReplicatesPerPeptideNodeThreshold(replicatesPerPeptideNodeThreshold);

		// are we collapsing the indisintuishable proteins
		final boolean collapseIndistinguishableProteins = Boolean
				.valueOf(properties.getProperty("collapseIndistinguishableProteins", "true"));
		params.setCollapseIndistinguishableProteins(collapseIndistinguishableProteins);
		// aminacids to quantify
		char[] aaQuantified = null;
		final String aaQuantifiedArrayString = properties.getProperty("collapsePeptidesBySites");
		if (aaQuantifiedArrayString != null && !"".equals(aaQuantifiedArrayString.trim())) {
			if (aaQuantifiedArrayString.contains(",")) {
				final String[] tmp = aaQuantifiedArrayString.split(",");
				aaQuantified = new char[tmp.length];
				for (int i = 0; i < tmp.length; i++) {
					aaQuantified[i] = tmp[i].trim().charAt(0);
				}
			} else {
				if (aaQuantifiedArrayString.length() > 1) {
					throw new FileNotFoundException(
							"Aminoacid array contains more than one character but is not separated by ','. Try something like 'K,R'");
				}
				aaQuantified = new char[1];
				aaQuantified[0] = aaQuantifiedArrayString.trim().charAt(0);
			}
		}
		params.setAaQuantified(aaQuantified);

		if (properties.containsKey("collapsePeptidesByPTMs")) {
			final String collapsePeptidesByPTMsString = properties.getProperty("collapsePeptidesByPTMs");
			try {
				params.setPTMQuantified(collapsePeptidesByPTMsString);
				log.info("Collapsing by " + params.getPTMsQuantified().size() + " ptm(s)");
			} catch (final Exception e) {
				throw new IllegalArgumentException("collapsePeptidesByPTMs is not recognized as '"
						+ collapsePeptidesByPTMsString + ". Posible values are TRUE or FALSE");
			}
		}
		if (params.isCollapseBySites() && params.isCollapseByPTMs()) {
			throw new IllegalArgumentException(
					"collapsePeptidesBySites and collapsePeptidesByPTMs cannot be used at the same time");
		}

		if (properties.containsKey("createProteinPTMStates")) {

			final boolean createProteinPTMStates = Boolean
					.valueOf(properties.getProperty("createProteinPTMStates", "false"));
			if (createProteinPTMStates && !params.isCollapseByPTMs()) {
				log.warn(
						"Found createProteinPTMStates=TRUE, but collapseByPTMs is not found or is FALSE, therefore, it will be ignored.");
				params.setCreateProteinPTMStates(false);
			} else {
				params.setCreateProteinPTMStates(createProteinPTMStates);
			}
		} else {
			if (params.isCollapseByPTMs()) {
				log.warn("createProteinPTMStates not found. Setting to default value FALSE");
			}
		}

		if (properties.containsKey("maxNumPTMsPerProtein")) {
			if (params.isCollapseByPTMs()) {
				final String maxNumPTMsPerProteinString = properties.getProperty("maxNumPTMsPerProtein");
				try {
					final Integer max = Integer.valueOf(maxNumPTMsPerProteinString);
					if (max < 0) {
						throw new Exception();
					}
					params.setMaxNumPTMsPerProtein(max);
				} catch (final Exception e) {
					throw new IllegalArgumentException("maxNumPTMsPerProtein is not recognized as '"
							+ maxNumPTMsPerProteinString + ". Posible values are positive integers");
				}
			} else {
				log.warn(
						"maxNumPTMsPerProtein parameter is provided but collapsePeptidesByPTMs is FALSE or not provided. maxNumPTMsPerProtein will be ignored.");
			}
		} else {
			if (params.isCollapseByPTMs()) {
				log.warn("collapsePeptidesByPTMs is TRUE but maxNumPTMsPerProtein is not provided. Default value of "
						+ params.getMaxNumPTMsPerProtein() + " will be used");
			}
		}

		if (properties.containsKey("isobaricRatioType")) {
			final String property = properties.getProperty("isobaricRatioType", false);
			try {
				final IsobaricRatioType isobaricRatioType = IsobaricRatioType.valueOf(property);
				params.setIsobaricRatioType(isobaricRatioType);
			} catch (final Exception e) {
				throw new IllegalArgumentException("isobaricRatioType is not recognized as '" + property
						+ ". Posible values are " + IsobaricRatioType.values());
			}
		}

		// are we collapsing the indistinguishable peptides
		final boolean collapseIndistinguishablePeptides = Boolean
				.valueOf(properties.getProperty("collapseIndistinguishablePeptides", "true"));
		params.setCollapseIndistinguishablePeptides(collapseIndistinguishablePeptides);

		if (params.isCollapseBySites() && !params.isCollapseIndistinguishablePeptides()) {
			throw new IllegalArgumentException(
					"collapsePeptidesBySites can only be used when collapseIndistinguishablePeptides is TRUE");
		}

		// determines if we align the peptides or not
		final boolean makeAlignments = Boolean.valueOf(properties.getProperty("makeAlignments", "false"));
		params.setMakeAlignments(makeAlignments);

		// DmDv, A, B, and E use 'K'
		// Human samples use 'K', 'R'
		final String enzymeArrayString = properties.getProperty("enzymeArray", "K,R");
		char[] enzymeArray = null;
		if (enzymeArrayString.contains(",")) {
			final String[] tmp = enzymeArrayString.split(",");
			enzymeArray = new char[tmp.length];
			for (int i = 0; i < tmp.length; i++) {
				enzymeArray[i] = tmp[i].trim().charAt(0);
			}
		} else {
			if (enzymeArrayString.length() > 1) {
				throw new FileNotFoundException(
						"Enzyme array contains more than one character but is not separated by ','. Try something like 'K,R'");
			}
			enzymeArray = new char[1];
			enzymeArray[0] = enzymeArrayString.trim().charAt(0);
		}
		params.setEnzymeArray(enzymeArray);
		// missedcleavages
		final int missedCleavages = Integer.valueOf(properties.getProperty("missedCleavages", "0"));
		params.setMissedCleavages(missedCleavages);

		final String peptideFilterRegexp = properties.getProperty("peptideFilterRegexp");
		if (peptideFilterRegexp != null) {
			// validate the regexp
			try {
				Pattern.compile(peptideFilterRegexp);
				params.setPeptideFilterRegexp(peptideFilterRegexp);
			} catch (final PatternSyntaxException e) {
				e.printStackTrace();
				throw e;
			}
		}
		final boolean useMayorityRule = Boolean
				.valueOf(properties.getProperty("useMajorityRulesForInfinities", "true"));
		params.setUseMayorityRule(useMayorityRule);

		final File uniprotReleasesFolder = new File(
				properties.getProperty("uniprotReleasesFolder", System.getProperty("user.dir")));
		params.setUniprotReleasesFolder(uniprotReleasesFolder);

		final String uniprotVersion = properties.getProperty("uniprotVersion", null);
		params.setUniprotVersion(uniprotVersion);

		// output prefix
		final String outputPrefix = properties.getProperty("outputPrefix", "");
		params.setOutputPrefix(outputPrefix);
		// output suffix
		final String outputSuffix = properties.getProperty("outputSuffix", "");
		params.setOutputSuffix(outputSuffix);
		// input file folder
		final String inputFilePathProperty = properties.getProperty("inputFilePath", System.getProperty("user.dir"));
		final File inputFileFolder = new File(FilenameUtils.separatorsToSystem(inputFilePathProperty));
		if (!inputFileFolder.exists()) {
			throw new FileNotFoundException(
					"Input file folder " + inputFileFolder.getAbsolutePath() + " doesn't exist");
		}
		params.setInputFileFolder(inputFileFolder);

		final String timeStamp = new SimpleDateFormat("yyyy-MM-dd").format(new Date());
		final String timeStamp2 = new SimpleDateFormat("yyyy-MM-dd-HH-mm-ss").format(new Date());
		final String outputFilePathProperty = properties.getProperty("outputFilePath", System.getProperty("user.dir"));
		params.setOutputFilePath(outputFilePathProperty);
		final File outputFileFolder = new File(outputFilePathProperty + File.separator + timeStamp + File.separator
				+ outputPrefix + "_" + outputSuffix + "_" + timeStamp2);
		final File temporalOutputFolder = new File(outputFileFolder.getAbsolutePath() + "_TEMP");
		params.setTemporalOutputFolder(temporalOutputFolder);

		params.setOutputFileFolder(outputFileFolder);

		// input files
		if (properties.containsKey("inputFiles")) {
			final String fileNamesString = properties.getProperty("inputFiles", false);
			parseInputFiles(fileNamesString, params);

		} else {
			log.info("Parameter 'inputFiles' not found. PCQ will not process quantitative values.");
		}

		// input files
		final String fileNamesString = properties.getProperty("inputIDFiles", false);
		if (fileNamesString != null) {
			parseInputIDFiles(fileNamesString, params);

		}
		if (params.getInputQuantificationFileNames().isEmpty()) {
			if (params.getIdentificationInputFileNamesArray().length > 0) {
				log.info("Working only with identification data, not quantitation.");
			} else {
				if (!batchMode) {
					throw new IllegalArgumentException(
							"One of the two 'inputFiles' or 'inputIDFiles' should be present and no empty");
				}
			}
		}

		// fasta file
		final String fastaFileProp = properties.getProperty("fastaFile", false);
		if (fastaFileProp != null && !"".equals(fastaFileProp)) {
			final File fastaFile = new File(fastaFileProp);
			if (!fastaFile.exists()) {
				throw new FileNotFoundException(fastaFile.getAbsolutePath() + " doesn't exist");
			}
			params.setFastaFile(fastaFile);
		}
		// ignore peptides not in indexed DB
		final boolean ignoreNotFoundPeptidesInDB = Boolean
				.valueOf(properties.getProperty("ignoreNotFoundPeptidesInDB", "false"));
		params.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);

		final String lightSpecies = properties.getProperty("lightSpecies", false);
		if (lightSpecies != null)
			params.setLightSpecies(lightSpecies);

		final String heavySpecies = properties.getProperty("heavySpecies", false);
		if (heavySpecies != null)
			params.setHeavySpecies(heavySpecies);

		// graphic options
		final ColorManager colorManager = ColorManager.getInstance();
		final String colorsByTax = properties.getProperty("colorsByTaxonomy");
		if (colorsByTax != null && !"".equals(colorsByTax)) {
			if (colorsByTax.contains(",")) {
				final String[] split = colorsByTax.split(",");
				for (int i = 0; i < split.length; i = i + 2) {
					final String tax = split[i];
					final String colorString = split[i + 1].trim();

					colorManager.addColorByTaxonomy(tax, colorString);
				}
			}
		}
		final String multiTaxonomyColor = properties.getProperty("colorMultiTaxonomy", "#ff66ff");
		if (multiTaxonomyColor != null && !"".equals(multiTaxonomyColor)) {
			colorManager.setMultiTaxonomyColor(multiTaxonomyColor);
		}
		// final String sharedNodeLabelColor =
		// properties.getProperty("shared_node_label_color", "#606060");
		// if (sharedNodeLabelColor != null && !"".equals(sharedNodeLabelColor))
		// {
		// colorManager.setSharedNodeLabelColor(sharedNodeLabelColor);
		// }
		if (makeAlignments) {
			final String alignedPeptideEdgecolor = properties.getProperty("colorAlignedPeptidesEdge", "#00ff00");
			if (alignedPeptideEdgecolor != null && !"".equals(alignedPeptideEdgecolor)) {
				colorManager.setAlignedPeptidesEdgeColor(alignedPeptideEdgecolor);
			}
		}
		params.setColorManager(colorManager);

		final String decoyRegexp = properties.getProperty("discardDecoys", null);
		params.setDecoyRegexp(decoyRegexp);

		if (makeAlignments) {
			final String finalAlignmentScore = properties.getProperty("finalAlignmentScore", "30");
			if (finalAlignmentScore != null) {
				params.setFinalAlignmentScore(Integer.valueOf(finalAlignmentScore));
			}
			final String sequenceIdentity = properties.getProperty("sequenceIdentity", "0.8");
			if (sequenceIdentity != null) {
				params.setSequenceIdentity(Double.valueOf(sequenceIdentity));
			}
			final String minConsecutiveIdenticalAlignment = properties.getProperty("minConsecutiveIdenticalAlignment",
					"6");
			if (minConsecutiveIdenticalAlignment != null) {
				params.setMinConsecutiveIdenticalAlignment(Integer.valueOf(minConsecutiveIdenticalAlignment));
			}
		}
		if (properties.containsKey("proteinLabel")) {
			final ProteinNodeLabel proteinLabel = ProteinNodeLabel
					.getFrom(properties.getProperty("proteinLabel", false));
			if (proteinLabel != null) {
				params.setProteinLabel(proteinLabel);
			} else {
				params.setProteinLabel(ProteinNodeLabel.ACC);
			}
		} else {
			params.setProteinLabel(ProteinNodeLabel.ACC);
		}

		final int proteinNodeWidth = Integer.valueOf(properties.getProperty("proteinNodeWidth", "70"));
		params.setProteinNodeWidth(proteinNodeWidth);
		final int proteinNodeHeigth = Integer.valueOf(properties.getProperty("proteinNodeHeight", "30"));
		params.setProteinNodeHeight(proteinNodeHeigth);
		final int peptideNodeWidth = Integer.valueOf(properties.getProperty("peptideNodeWidth", "70"));
		params.setPeptideNodeWidth(peptideNodeWidth);
		final int peptideNodeHeigth = Integer.valueOf(properties.getProperty("peptideNodeHeight", "30"));
		params.setPeptideNodeHeight(peptideNodeHeigth);

		final Shape proteinShape = Shape.valueOf(properties.getProperty("proteinNodeShape", "ELLIPSE"));
		params.setProteinNodeShape(proteinShape);
		final Shape peptideShape = Shape.valueOf(properties.getProperty("peptideNodeShape", "ROUNDRECT"));
		params.setPeptideNodeShape(peptideShape);

		final Color colorRatioMin = ColorGenerator.hex2Rgb(properties.getProperty("colorRatioMin", "#0000ff"));
		params.setColorRatioMin(colorRatioMin);

		final Color colorRatioMax = ColorGenerator.hex2Rgb(properties.getProperty("colorRatioMax", "#ff0000"));
		params.setColorRatioMax(colorRatioMax);

		final String colorNonRegulatedString = properties.getProperty("colorNonRegulatedPeptides", false);
		if (colorNonRegulatedString != null) {
			final Color colorNonRegulated = ColorGenerator.hex2Rgb(colorNonRegulatedString);
			params.setColorNonRegulated(colorNonRegulated);
		}

		final double minimumRatioForColor = Double.valueOf(properties.getProperty("minimumRatioForColor", "-10"));
		params.setMinimumRatioForColor(minimumRatioForColor);
		final double maximumRatioForColor = Double.valueOf(properties.getProperty("maximumRatioForColor", "10"));
		params.setMaximumRatioForColor(maximumRatioForColor);

		// show cases in edges
		final boolean showCasesInEdges = Boolean.valueOf(properties.getProperty("showCasesInEdges", "true"));
		params.setShowCasesInEdges(showCasesInEdges);

		// remark significant peptides
		final boolean remarkSignificantPeptides = Boolean
				.valueOf(properties.getProperty("remarkSignificantPeptides", "true"));
		params.setRemarkSignificantPeptides(remarkSignificantPeptides);

		// MONGO DB
		final String mongoDBHostURI = properties.getProperty("mongoDBHostURI");
		params.setMongoDBURI(mongoDBHostURI);

		final String mongoProtDBName = properties.getProperty("mongoProtDBName");
		params.setMongoProtDBName(mongoProtDBName);

		final String mongoMassDBName = properties.getProperty("mongoMassDBName");
		params.setMongoMassDBName(mongoMassDBName);

		final String mongoSeqDBName = properties.getProperty("mongoSeqDBName");
		params.setMongoSeqDBName(mongoSeqDBName);

		// SANXOT
		final boolean performRatioIntegration = Boolean
				.valueOf(properties.getProperty("performRatioIntegration", "false"));
		params.setPerformRatioIntegration(performRatioIntegration);
		if (performRatioIntegration) {
			// sanxot scripts paths
			final String sanxotPath = properties.getProperty("sanxotPath",
					System.getProperty("user.dir") + File.separator + "SanXot");
			final File scriptsPath = new File(sanxotPath);
			params.setSanXotPath(scriptsPath);
			try {
				if (properties.containsKey("outliersRemovalFDR")) {
					final double outlierRemovalFDR = Double
							.valueOf(properties.getProperty("outliersRemovalFDR", false));
					params.setOutliersRemovalFDR(outlierRemovalFDR);
				}
			} catch (final NumberFormatException e) {
				// do nothing
			}
		}
		try {
			if (properties.containsKey("significantFDRThreshold")) {
				final double significantFDRThreshold = Double
						.valueOf(properties.getProperty("significantFDRThreshold", false));
				params.setSignificantFDRThreshold(significantFDRThreshold);
			}
		} catch (final NumberFormatException e) {
			// do nothing
		}

		final boolean removeFilteredNodes = Boolean.valueOf(properties.getProperty("removeFilteredNodes", "true"));
		params.setRemoveFilteredNodes(removeFilteredNodes);

		if (properties.containsKey("ignorePTMs")) {
			final boolean ignorePTMs = Boolean.valueOf(properties.getProperty("ignorePTMs", "true"));
			params.setIgnorePTMs(ignorePTMs);
		} else {
			log.info("ignorePTMs param not found. Setting default value TRUE");
			params.setIgnorePTMs(true);
		}

		final boolean semiCleavage = Boolean.valueOf(properties.getProperty("semiCleavage", "false"));
		params.setSemiCleavage(semiCleavage);

		final String tripletsString = properties.getProperty("uniprot_xpath", false);
		if (tripletsString != null && !"".contentEquals(tripletsString)) {
			final List<String> triplets = parseTriplets(tripletsString);
			for (final String triplet : triplets) {
				final UniprotAnnotationColumn column = new UniprotAnnotationColumn(triplet);
				params.addUniprotAnnotationColumn(column);
			}
		}

		// look in uniprot for proteoforms
		final boolean lookForProteoforms = Boolean
				.valueOf(properties.getProperty("lookInUniprotForProteoforms", "false"));
		params.setLookForProteoforms(lookForProteoforms);

		// ignore taxonomies
		final boolean ignoreTaxonomies = Boolean.valueOf(properties.getProperty("ignoreTaxonomies", "false"));
		params.setIgnoreTaxonomies(ignoreTaxonomies);

		// ignore ACC format
		final boolean ignoreACCFormat = !Boolean.valueOf(properties.getProperty("recognizeACCFormat", "true"));
		params.setIgnoreACCFormat(ignoreACCFormat);

		// print PTM position in protein column
		final boolean printPTMPositionInProtein = Boolean
				.valueOf(properties.getProperty("printPTMPositionInProtein", "true"));
		params.setPrintPTMPositionInProtein(printPTMPositionInProtein);

		// print PTM position in protein column
		final boolean writePSEAQuantInputFiles = Boolean
				.valueOf(properties.getProperty("writePSEAQuantInputFiles", "true"));
		params.setWritePSEAQuantInputFiles(writePSEAQuantInputFiles);

		// check errors
		checkErrorsInParameters(params);
	}

	public static void parseInputIDFiles(String fileNamesString, ProteinClusterQuantParameters params) {
		if (fileNamesString.contains("|")) {
			final String[] tmp = fileNamesString.split("\\|");
			for (int i = 0; i < tmp.length; i++) {
				final ExperimentFiles experimentFiles = parseExperimentFileNames(tmp[i].trim());
				params.addIdentificationInputFileNames(experimentFiles);
			}
		} else {
			if (!"".equals(fileNamesString)) {
				final ExperimentFiles experimentFiles = parseExperimentFileNames(fileNamesString);
				params.addIdentificationInputFileNames(experimentFiles);
			}
		}
	}

	public static void parseInputFiles(String fileNamesString, ProteinClusterQuantParameters params) {
		if (fileNamesString != null && fileNamesString.contains("|")) {
			final String[] tmp = fileNamesString.split("\\|");
			for (int i = 0; i < tmp.length; i++) {
				final ExperimentFiles experimentFiles = parseExperimentFileNames(tmp[i].trim());
				params.addQuantificationInputFileNames(experimentFiles);
			}
		} else {
			if (fileNamesString != null && !"".equals(fileNamesString)) {
				final ExperimentFiles experimentFiles = parseExperimentFileNames(fileNamesString);
				params.addQuantificationInputFileNames(experimentFiles);
			} else {
				log.info("Parameter 'inputFiles' not found. PCQ will not process quantitative values.");
			}
		}
	}

	public static List<String> parseTriplets(String tripletsString) {
		final List<String> ret = new ArrayList<String>();
		final StringBuilder triplet = new StringBuilder();
		boolean insideBrackets = false;
		for (int i = 0; i < tripletsString.length(); i++) {
			final char charAt = tripletsString.charAt(i);
			if (charAt == '[') {
				insideBrackets = true;
				continue;
			}
			if (charAt == ']') {
				if (!"".equals(triplet.toString())) {
					ret.add(triplet.toString());
				}
				insideBrackets = false;
				continue;
			}
			if (insideBrackets) {
				triplet.append(charAt);
			}
		}
		return ret;
	}

	/**
	 * Parse a strign like:<br>
	 * experiment_name [ replicate_name_1, replicate_name_2] <br>
	 * and put it in the returning map
	 *
	 * @param trim
	 * @return
	 */
	public static ExperimentFiles parseExperimentFileNames(String string) {
		ExperimentFiles ret = null;
		try {
			final int firstBracketPosition = string.indexOf("[");
			final int secondBracketPosition = string.indexOf("]");
			final String experimentName = string.substring(0, firstBracketPosition).trim();
			ret = new ExperimentFiles(experimentName);
			final String replicatesString = string.substring(firstBracketPosition + 1, secondBracketPosition).trim();
			if (replicatesString.contains(",")) {
				final String[] split = replicatesString.split(",");
				for (final String replicateName : split) {
					ret.addReplicateFileName(replicateName.trim());
				}
			} else {
				ret.addReplicateFileName(replicatesString);
			}
		} catch (final Exception e) {
			if (e instanceof IllegalArgumentException) {
				throw (IllegalArgumentException) e;
			}
			throw new IllegalArgumentException("File name string is not well formed: '" + string
					+ "'\nTry something like: exp1[rep11_file.xml, rep12_file.xml] | exp2[rep21_file.xml, rep22_file.xml]");
		}
		return ret;
	}

	private static void checkErrorsInParameters(ProteinClusterQuantParameters params) {
		// not allow collapsing peptides and align peptides
		if (params.isCollapseIndistinguishablePeptides() && params.isMakeAlignments()) {
			throw new IllegalArgumentException(
					"ERROR in input parameters: It is not possible to collapse the peptides and make the alignments in the same network analysis");
		}
		if (params.isPerformRatioIntegration() && !params.getQuantParameters().getSanxotScriptsFolder().exists()) {
			throw new IllegalArgumentException(
					params.getQuantParameters().getSanxotScriptsFolder().getAbsolutePath() + " does not exist");
		}
		if (params.getAaQuantified() != null && params.getAaQuantified().length > 0
				&& params.isApplyClassificationsByProteinPair()) {
			throw new IllegalArgumentException(
					"ERROR in input parameters: It is not possible to perform a protein pair analysis on site specific quantification results");
		}
		if (params.getAnalysisInputType() == AnalysisInputType.CENSUS_CHRO && params.getIsobaricRatioType() == null) {
			throw new IllegalArgumentException(
					"ERROR in input parameters: You need to specify 'isobaricRatioType' parameter when using inputType=CENSUS_CHRO");
		}
		if (params.getAnalysisInputType() != AnalysisInputType.CENSUS_CHRO && params.getIsobaricRatioType() != null) {
			throw new IllegalArgumentException(
					"ERROR in input parameters: 'isobaricRatioType' parameter is not valid for inputType="
							+ params.getAnalysisInputType() + ". It is only valid for inputType="
							+ AnalysisInputType.CENSUS_CHRO);
		}

	}

}
