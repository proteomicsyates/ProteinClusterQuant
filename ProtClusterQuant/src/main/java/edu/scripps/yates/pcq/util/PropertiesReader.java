package edu.scripps.yates.proteinclusters.util;

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
		double thresholdForSignificance = Double.valueOf(properties.getProperty("thresholdForSignificance", "0"));
		params.setThresholdForSignificance(thresholdForSignificance);
		boolean printOnlyFirstGene = Boolean.valueOf(properties.getProperty("printOnlyFirstGene", "true"));
		params.setPrintOnlyFirstGene(printOnlyFirstGene);
		// do we have threhsold for minimum PSMs per Peptide
		boolean ionsPerPeptideThresholdOn = Boolean
				.valueOf(properties.getProperty("ionsPerPeptideThresholdOn", "false"));
		params.setIonsPerPeptideThresholdOn(ionsPerPeptideThresholdOn);
		// threshold for minimum ions per peptide
		int ionsPerPeptideThreshold = Integer.valueOf(properties.getProperty("ionsPerPeptideThreshold", "8"));
		params.setIonsPerPeptideThreshold(ionsPerPeptideThreshold);
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
		String[] inputFileNames = null;
		String fileNamesString = properties.getProperty("input_files");
		if (fileNamesString.contains(",")) {
			String[] tmp = fileNamesString.split(",");
			inputFileNames = new String[tmp.length];
			for (int i = 0; i < tmp.length; i++) {
				inputFileNames[i] = tmp[i].trim();
			}
		} else {
			inputFileNames = new String[1];
			inputFileNames[0] = fileNamesString.trim();
		}
		params.setInputFileNames(inputFileNames);
		final String fastaFileProp = properties.getProperty("fasta_file", false);
		if (fastaFileProp != null) {
			File fastaFile = new File(fastaFileProp);
			if (!fastaFile.exists()) {
				throw new FileNotFoundException(fastaFile.getAbsolutePath() + " doesn't exist");
			}
			params.setFastaFile(fastaFile);
		}

		String lightSpecies = properties.getProperty("light_species", false);
		if (lightSpecies != null)
			params.setLightSpecies(lightSpecies);

		String heavySpecies = properties.getProperty("heavy_species", false);
		if (heavySpecies != null)
			params.setHeavySpecies(heavySpecies);

		//////////////
		// PSEA QUANT
		// replicate_identifiers
		String replicateIdentifiersString = properties.getProperty("replicate_identifiers", false);
		if (replicateIdentifiersString != null) {
			String[] replicateIdentifiers = null;
			if (replicateIdentifiersString.contains(",")) {
				String[] tmp = replicateIdentifiersString.split(",");
				replicateIdentifiers = new String[tmp.length];
				for (int i = 0; i < tmp.length; i++) {
					replicateIdentifiers[i] = tmp[i].trim();
				}
			} else {
				replicateIdentifiers = new String[1];
				replicateIdentifiers[0] = replicateIdentifiersString.trim();
			}
			params.setReplicateIdentifiers(replicateIdentifiers);
		}
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
		// check errors
		checkErrorsInParameters(params);
	}

	private static void checkErrorsInParameters(ProteinClusterQuantParameters params) {
		// not allow collapsing peptides and align peptides
		if (params.isCollapseIndistinguishablePeptides() && params.isMakeAlignments()) {
			throw new IllegalArgumentException(
					"ERROR in input parameters: It is not possible to collapse the peptides and make the alignments in the same network analysis");
		}

	}
}
