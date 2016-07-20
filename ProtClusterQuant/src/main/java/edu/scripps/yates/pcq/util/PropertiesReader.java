package edu.scripps.yates.pcq.util;

import java.awt.Color;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.text.SimpleDateFormat;
import java.util.Date;

public class PropertiesReader {
	private static final String PROPERTIES_FILE_NAME = "setup.properties";
	private final static org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(PropertiesReader.class);

	public static ProteinClusterQuantProperties readerProperties() {
		ClassLoader cl = PropertiesReader.class.getClassLoader();
		InputStream inputStream = cl.getResourceAsStream(PROPERTIES_FILE_NAME);
		return readProperties(inputStream);
	}

	public static ProteinClusterQuantProperties readProperties(File setupPropertiesFile) {
		if (setupPropertiesFile == null || !setupPropertiesFile.exists()) {
			throw new IllegalArgumentException("setup properties file not valid or null");
		}
		try {
			InputStream inputStream = new FileInputStream(setupPropertiesFile);
			ProteinClusterQuantProperties prop = readProperties(inputStream);
			return prop;
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(e);
		}
	}

	public static ProteinClusterQuantProperties readProperties(InputStream inputStream) {

		if (inputStream == null) {
			throw new IllegalArgumentException("input stream is null");
		}
		try {
			ProteinClusterQuantProperties prop = new ProteinClusterQuantProperties();
			prop.load(inputStream);
			readParametersFromProperties(prop);
			return prop;
		} catch (IOException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(e);
		}

	}

	private static void readParametersFromProperties(ProteinClusterQuantProperties properties)
			throws FileNotFoundException {
		final ProteinClusterQuantParameters params = ProteinClusterQuantParameters.getInstance();
		try {
			InputType inputType = InputType.valueOf(properties.getProperty("inputType", true));
			params.setInputType(inputType);
		} catch (Exception e) {
			throw new IllegalArgumentException(
					"'inputType' parameter can only have the following values: " + InputType.getPossibleValues());
		}

		if (properties.containsKey("onlyOneSpectrumPerChromatographicPeakAndPerSaltStep")) {
			boolean onlyOneSpectrumPerChromatographicPeakAndPerSaltStep = Boolean
					.valueOf(properties.getProperty("onlyOneSpectrumPerChromatographicPeakAndPerSaltStep", false));
			params.setOnlyOneSpectrumPerChromatographicPeakAndPerSaltStep(
					onlyOneSpectrumPerChromatographicPeakAndPerSaltStep);
		}
		if (properties.containsKey("skipSingletons")) {
			boolean skipSingletons = Boolean.valueOf(properties.getProperty("skipSingletons", false));
			params.setSkipSingletons(skipSingletons);
		}
		// CASIMIR
		// protein pair analysis based on significant value cutoff
		boolean significantProteinPairAnalysis = Boolean
				.valueOf(properties.getProperty("significantProteinPairAnalysis", "false"));
		params.setSignificantProteinPairAnalysis(significantProteinPairAnalysis);
		// CASIMIR
		// protein pair analysis based on significant value cutoff
		boolean labelSwap = Boolean.valueOf(properties.getProperty("labelSwap", "false"));
		params.setLabelSwap(labelSwap);
		// CASIMIR
		// Standard deviation of peptide measurements as significance cutoff
		boolean stdAsSignficanceCutoff = Boolean.valueOf(properties.getProperty("stdAsSignficanceCutoff", "false"));
		params.setStdAsSignficanceCutoffOn(stdAsSignficanceCutoff);
		// CASIMIR
		// cutoff for significance
		double thresholdForSignificance = Double.valueOf(properties.getProperty("thresholdForSignificance", "0.05"));
		params.setThresholdForSignificance(thresholdForSignificance);
		boolean printOnlyFirstGene = Boolean.valueOf(properties.getProperty("printOnlyFirstGene", "true"));
		params.setPrintOnlyFirstGene(printOnlyFirstGene);
		// do we have threhsold for minimum ion counts per Peptide
		boolean ionsPerPeptideNodeThresholdOn = Boolean
				.valueOf(properties.getProperty("ionsPerPeptideNodeThresholdOn", "false"));
		params.setIonsPerPeptideNodeThresholdOn(ionsPerPeptideNodeThresholdOn);
		// threshold for minimum ions per peptide node
		int ionsPerPeptideNodeThreshold = Integer.valueOf(properties.getProperty("ionsPerPeptideNodeThreshold", "8"));
		params.setIonsPerPeptideNodeThreshold(ionsPerPeptideNodeThreshold);
		// threshold for minimum PSMs per Peptide
		boolean psmsPerPeptideNodeThresholdOn = Boolean
				.valueOf(properties.getProperty("psmsPerPeptideNodeThresholdOn", "false"));
		params.setPsmsPerPeptideThresholdOn(psmsPerPeptideNodeThresholdOn);
		// threshold for minimum PSMs per peptide node
		int psmsPerPeptideNodeThreshold = Integer.valueOf(properties.getProperty("psmsPerPeptideNodeThreshold", "0"));
		params.setPsmsPerPeptideThreshold(psmsPerPeptideNodeThreshold);

		// threshold for the iglewiczHoaglin Test. A result wquals or greater
		// than that, would be considered as an outlier
		double iglewiczHoaglinTestThreshold = Double.valueOf(properties.getProperty("iglewiczHoaglinTest", "3.5"));
		params.setIglewiczHoaglinTestThreshold(iglewiczHoaglinTestThreshold);
		// are we collapsing the indisintuishable proteins
		boolean collapseIndistinguishableProteins = Boolean
				.valueOf(properties.getProperty("collapseIndistinguishableProteins", "true"));
		params.setCollapseIndistinguishableProteins(collapseIndistinguishableProteins);
		// are we collapsing the indisintuishable peptides
		boolean collapseIndistinguishablePeptides = Boolean
				.valueOf(properties.getProperty("collapseIndistinguishablePeptides", "true"));
		params.setCollapseIndistinguishablePeptides(collapseIndistinguishablePeptides);
		// determines if we align the peptides or not
		boolean makeAlignments = Boolean.valueOf(properties.getProperty("makeAlignments", "true"));
		params.setMakeAlignments(makeAlignments);
		// if prints out data for k-means clustering
		boolean printKMeans = Boolean.valueOf(properties.getProperty("k_means", "true"));
		params.setPrintKMeans(printKMeans);
		// DmDv, A, B, and E use 'K'
		// Human samples use 'K', 'R'
		final String enzymeArrayString = properties.getProperty("enzymeArray", "K,R");
		char[] enzymeArray = null;
		if (enzymeArrayString.contains(",")) {
			String[] tmp = enzymeArrayString.split(",");
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
		int missedCleavages = Integer.valueOf(properties.getProperty("missedcleavages", "1"));
		params.setMissedCleavages(missedCleavages);
		// do we only count truly unique peptides as unique
		boolean uniquePepOnly = Boolean.valueOf(properties.getProperty("uniquePepOnly", "true"));
		params.setUniquePepOnly(uniquePepOnly);
		// if want to get Taxonomy
		boolean getTax = Boolean.valueOf(properties.getProperty("getTax", "false"));
		params.setGetTax(getTax);
		File uniprotReleasesFolder = new File(
				properties.getProperty("uniprotReleasesFolder", System.getProperty("user.dir")));
		params.setUniprotReleasesFolder(uniprotReleasesFolder);

		String uniprotVersion = properties.getProperty("uniprotVersion", null);
		params.setUniprotVersion(uniprotVersion);

		// output prefix
		String outputPrefix = properties.getProperty("output_prefix", true);
		params.setOutputPrefix(outputPrefix);
		// output suffix
		String outputSuffix = properties.getProperty("output_sufix", true);
		params.setOutputSuffix(outputSuffix);
		// input file folder
		File inputFileFolder = new File(properties.getProperty("input_file_path", System.getProperty("user.dir")));
		if (!inputFileFolder.exists()) {
			throw new FileNotFoundException(
					"Input file folder " + inputFileFolder.getAbsolutePath() + " doesn't exist");
		}
		params.setInputFileFolder(inputFileFolder);

		String timeStamp = new SimpleDateFormat("yyyy-MM-dd-HH-mm-ss").format(new Date());
		File outputFileFolder = new File(properties.getProperty("output_file_path", System.getProperty("user.dir"))
				+ File.separator + new SimpleDateFormat("yyyy-MM-dd").format(new Date()) + File.separator + outputPrefix
				+ "_" + outputSuffix + "_" + timeStamp);
		File temporalOutputFolder = new File(outputFileFolder.getAbsolutePath() + "_TEMP");
		params.setTemporalOutputFolder(temporalOutputFolder);
		if (!temporalOutputFolder.exists()) {
			// create it
			log.info("Creating temporal output folder at: " + outputFileFolder.getAbsolutePath());
			temporalOutputFolder.mkdirs();
		}
		params.setOutputFileFolder(outputFileFolder);

		// input files
		String fileNamesString = properties.getProperty("input_files");
		if (fileNamesString.contains("|")) {
			String[] tmp = fileNamesString.split("\\|");
			for (int i = 0; i < tmp.length; i++) {
				ExperimentFiles experimentFiles = parseExperimentFileNames(tmp[i].trim());
				params.addInputFileNames(experimentFiles);
			}
		} else {
			ExperimentFiles experimentFiles = parseExperimentFileNames(fileNamesString);
			params.addInputFileNames(experimentFiles);
		}

		// fasta file
		final String fastaFileProp = properties.getProperty("fasta_file", false);
		if (fastaFileProp != null && !"".equals(fastaFileProp)) {
			File fastaFile = new File(fastaFileProp);
			if (!fastaFile.exists()) {
				throw new FileNotFoundException(fastaFile.getAbsolutePath() + " doesn't exist");
			}
			params.setFastaFile(fastaFile);
		}
		// ignore peptides not in indexed DB
		boolean ignoreNotFoundPeptidesInDB = Boolean
				.valueOf(properties.getProperty("ignoreNotFoundPeptidesInDB", "false"));
		params.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);

		String lightSpecies = properties.getProperty("light_species", false);
		if (lightSpecies != null)
			params.setLightSpecies(lightSpecies);

		String heavySpecies = properties.getProperty("heavy_species", false);
		if (heavySpecies != null)
			params.setHeavySpecies(heavySpecies);

		//////////////
		// PSEA QUANT
		// replicate_identifiers

		boolean excludeUniquePeptides = Boolean.valueOf(properties.getProperty("excludeUniquePeptides", "true"));
		params.setExcludeUniquePeptides(excludeUniquePeptides);
		/////////////////

		// graphic options
		ColorManager colorManager = ColorManager.getInstance();
		final String colorsByTax = properties.getProperty("colors_by_taxonomy");
		if (colorsByTax != null && !"".equals(colorsByTax)) {
			if (colorsByTax.contains(",")) {
				final String[] split = colorsByTax.split(",");
				for (int i = 0; i < split.length; i = i + 2) {
					String tax = split[i];
					String colorString = split[i + 1];

					colorManager.addColorByTaxonomy(tax, colorString);
				}
			}
		}
		final String multiTaxonomyColor = properties.getProperty("color_multi_taxonomy", "#ff66ff");
		if (multiTaxonomyColor != null && !"".equals(multiTaxonomyColor)) {
			colorManager.setMultiTaxonomyColor(multiTaxonomyColor);
		}
		final String sharedNodeLabelColor = properties.getProperty("shared_node_label_color", "#606060");
		if (sharedNodeLabelColor != null && !"".equals(sharedNodeLabelColor)) {
			colorManager.setSharedNodeLabelColor(sharedNodeLabelColor);
		}
		final String alignedPeptideEdgecolor = properties.getProperty("color_aligned_peptides_edge", "#00ff00");
		if (alignedPeptideEdgecolor != null && !"".equals(alignedPeptideEdgecolor)) {
			colorManager.setAlignedPeptidesEdgeColor(alignedPeptideEdgecolor);
		}
		params.setColorManager(colorManager);

		String decoyRegexp = properties.getProperty("discardDecoys", null);
		params.setDecoyRegexp(decoyRegexp);

		String finalAlignmentScore = properties.getProperty("final_alignment_score", "30");
		if (finalAlignmentScore != null)
			params.setFinalAlignmentScore(Integer.valueOf(finalAlignmentScore));

		String sequenceIdentity = properties.getProperty("sequence_identity", "0.8");
		if (sequenceIdentity != null)
			params.setSequenceIdentity(Double.valueOf(sequenceIdentity));

		String maxConsecutiveIdenticalAlignment = properties.getProperty("max_consecutive_identical_alignment", "6");
		if (maxConsecutiveIdenticalAlignment != null)
			params.setMaxConsecutiveIdenticalAlignment(Integer.valueOf(maxConsecutiveIdenticalAlignment));

		ProteinLabel proteinLabel = ProteinLabel.getFrom(properties.getProperty("protein_label", false));
		if (proteinLabel != null) {
			params.setProteinLabel(proteinLabel);
		} else {
			params.setProteinLabel(ProteinLabel.ACC);
		}

		int proteinNodeWidth = Integer.valueOf(properties.getProperty("protein_node_width", "70"));
		params.setProteinNodeWidth(proteinNodeWidth);
		int proteinNodeHeigth = Integer.valueOf(properties.getProperty("protein_node_height", "30"));
		params.setProteinNodeHeight(proteinNodeHeigth);
		int peptideNodeWidth = Integer.valueOf(properties.getProperty("peptide_node_width", "70"));
		params.setPeptideNodeWidth(peptideNodeWidth);
		int peptideNodeHeigth = Integer.valueOf(properties.getProperty("peptide_node_height", "30"));
		params.setPeptideNodeHeight(peptideNodeHeigth);

		Shape proteinShape = Shape.valueOf(properties.getProperty("protein_node_shape", "ELLIPSE"));
		params.setProteinNodeShape(proteinShape);
		Shape peptideShape = Shape.valueOf(properties.getProperty("peptide_node_shape", "ROUNDRECT"));
		params.setPeptideNodeShape(peptideShape);

		Color colorRatioMin = ColorManager.hex2Rgb(properties.getProperty("color_ratio_min", "#00ff00"));
		params.setColorRatioMin(colorRatioMin);

		Color colorRatioMax = ColorManager.hex2Rgb(properties.getProperty("color_ratio_max", "#ff0000"));
		params.setColorRatioMax(colorRatioMax);

		final String colorNonRegulatedString = properties.getProperty("color_non_regulated", false);
		if (colorNonRegulatedString != null) {
			Color colorNonRegulated = ColorManager.hex2Rgb(colorNonRegulatedString);
			params.setColorNonRegulated(colorNonRegulated);
		}

		double minimumRatioForColor = Double.valueOf(properties.getProperty("minimum_ratio_for_color", "-10"));
		params.setMinimumRatioForColor(minimumRatioForColor);
		double maximumRatioForColor = Double.valueOf(properties.getProperty("maximum_ratio_for_color", "10"));
		params.setMaximumRatioForColor(maximumRatioForColor);

		// show cases in edges
		boolean showCasesInEdges = Boolean.valueOf(properties.getProperty("show_cases_in_edges", "true"));
		params.setShowCasesInEdges(showCasesInEdges);

		// remark significant peptides
		boolean remarkSignificantPeptides = Boolean
				.valueOf(properties.getProperty("remark_significant_peptides", "true"));
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
		boolean performRatioIntegration = Boolean.valueOf(properties.getProperty("perform_ratio_integration", "false"));
		params.setPerformRatioIntegration(performRatioIntegration);
		try {
			if (properties.containsKey("outliers_removal_FDR")) {
				double outlierRemovalFDR = Double.valueOf(properties.getProperty("outliers_removal_FDR", false));
				params.setOutliersRemovalFDR(outlierRemovalFDR);
			}
		} catch (NumberFormatException e) {
			// do nothing
		}
		try {
			if (properties.containsKey("significant_FDR_threshold")) {
				double significantFDRThreshold = Double
						.valueOf(properties.getProperty("significant_FDR_threshold", false));
				params.setSignificantFDRThreshold(significantFDRThreshold);
			}
		} catch (NumberFormatException e) {
			// do nothing
		}

		// check errors
		checkErrorsInParameters(params);
	}

	/**
	 * Parse a strign like:<br>
	 * experiment_name [ replicate_name_1, replicate_name_2] <br>
	 * and put it in the returning map
	 *
	 * @param trim
	 * @return
	 */
	private static ExperimentFiles parseExperimentFileNames(String string) {
		ExperimentFiles ret = null;
		try {
			int firstBracketPosition = string.indexOf("[");
			int secondBracketPosition = string.indexOf("]");
			String experimentName = string.substring(0, firstBracketPosition).trim();
			ret = new ExperimentFiles(experimentName);
			final String replicatesString = string.substring(firstBracketPosition + 1, secondBracketPosition).trim();
			if (replicatesString.contains(",")) {
				final String[] split = replicatesString.split(",");
				for (String replicateName : split) {
					ret.addReplicateFileName(replicateName.trim());
				}
			} else {
				ret.addReplicateFileName(replicatesString);
			}
		} catch (Exception e) {
			if (e instanceof IllegalArgumentException) {
				throw e;
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

	}
}
