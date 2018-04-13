
package edu.scripps.yates.pcq.util;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.annotations.uniprot.xml.Entry;
import edu.scripps.yates.annotations.uniprot.xml.GeneNameType;
import edu.scripps.yates.annotations.uniprot.xml.GeneType;
import edu.scripps.yates.annotations.uniprot.xml.OrganismNameType;
import edu.scripps.yates.annotations.uniprot.xml.OrganismType;
import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.CensusChroParser;
import edu.scripps.yates.census.read.CensusOutParser;
import edu.scripps.yates.census.read.SeparatedValuesParser;
import edu.scripps.yates.census.read.model.CensusRatio;
import edu.scripps.yates.census.read.model.Ion;
import edu.scripps.yates.census.read.model.IonCountRatio;
import edu.scripps.yates.census.read.model.IonSerie.IonSerieType;
import edu.scripps.yates.census.read.model.IsoRatio;
import edu.scripps.yates.census.read.model.IsobaricQuantifiedPSM;
import edu.scripps.yates.census.read.model.IsobaricQuantifiedPeptide;
import edu.scripps.yates.census.read.model.QuantifiedPSM;
import edu.scripps.yates.census.read.model.QuantifiedPeptide;
import edu.scripps.yates.census.read.model.RatioScore;
import edu.scripps.yates.census.read.model.interfaces.QuantParser;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.dbindex.DBIndexInterface;
import edu.scripps.yates.dbindex.io.DBIndexSearchParams;
import edu.scripps.yates.dbindex.io.DBIndexSearchParamsImpl;
import edu.scripps.yates.dbindex.util.IndexUtil;
import edu.scripps.yates.dbindex.util.PeptideFilterBySequence;
import edu.scripps.yates.dtaselectparser.DTASelectParser;
import edu.scripps.yates.pcq.model.IsobaricRatioType;
import edu.scripps.yates.pcq.model.PCQPeptideNode;
import edu.scripps.yates.pcq.model.PCQProteinNode;
import edu.scripps.yates.pcq.model.ProteinCluster;
import edu.scripps.yates.pcq.params.ProteinClusterQuantParameters;
import edu.scripps.yates.pcq.xgmml.util.AlignedPeptides;
import edu.scripps.yates.pcq.xgmml.util.AlignmentSet;
import edu.scripps.yates.utilities.alignment.nwalign.NWAlign;
import edu.scripps.yates.utilities.alignment.nwalign.NWResult;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.model.enums.AggregationLevel;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;
import edu.scripps.yates.utilities.sequence.PositionInPeptide;
import edu.scripps.yates.utilities.sequence.PositionInProtein;
import edu.scripps.yates.utilities.strings.StringUtils;
import edu.scripps.yates.utilities.util.Pair;
import edu.scripps.yates.utilities.util.StringPosition;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;
import gnu.trove.set.hash.TIntHashSet;

public class PCQUtils {

	public static final String PROTEIN_DESCRIPTION_SEPARATOR = "####";
	public static DecimalFormat df = new DecimalFormat("#.#");
	private static Map<String, DBIndexInterface> indexByFastaIndexKey = new THashMap<String, DBIndexInterface>();
	private static final Map<String, UniprotProteinLocalRetriever> uplrMap = new THashMap<String, UniprotProteinLocalRetriever>();
	private final static Logger log = Logger.getLogger(PCQUtils.class);
	public static final String PROTEIN_ACC_SEPARATOR = " ";
	public static final double factor = 1.2;
	public static final String FDR_CONFIDENCE_SCORE_NAME = "FDR";
	private static Map<String, QuantParser> quantParsersByFileNamesKey = new THashMap<String, QuantParser>();
	private static Map<String, DTASelectParser> dtaSelectParsersByFileNamesKey = new THashMap<String, DTASelectParser>();
	public static final String KEY_SEPARATOR = "-";

	public static AlignmentSet alignPeptides(List<QuantifiedPeptideInterface> peptideList, QuantCondition cond1,
			QuantCondition cond2, FileWriter alignmentLogFile) throws IOException {
		final AlignmentSet ret = new AlignmentSet();
		boolean someAlignmentPassThresholds = false;
		long totalAligments = 0;
		try {
			final int finalAlignmentScore = ProteinClusterQuantParameters.getInstance().getFinalAlignmentScore();
			int maxAlignmentScore = -Integer.MAX_VALUE;
			final double sequenceIdentity = ProteinClusterQuantParameters.getInstance().getSequenceIdentity();
			double maxSeqIdentity = -Double.MAX_VALUE;
			final int minConsecutiveIdenticalAlignment = ProteinClusterQuantParameters.getInstance()
					.getMinConsecutiveIdenticalAlignment();
			int maxConsecutiveIdenticalAlignment = -Integer.MAX_VALUE;
			log.info("Aligning " + peptideList.size() + " peptides between them");
			log.info("Using minScore=" + finalAlignmentScore + ", minSeqIdentity=" + sequenceIdentity
					+ ", minConsecutiveAlignment=" + minConsecutiveIdenticalAlignment);
			final ProgressCounter counter = new ProgressCounter(peptideList.size(),
					ProgressPrintingType.PERCENTAGE_STEPS, 0);
			for (int i = 0; i < peptideList.size(); i++) {
				counter.increment();
				final QuantifiedPeptideInterface pep1 = peptideList.get(i);
				final String printIfNecessary = counter.printIfNecessary();
				if (printIfNecessary != null && !"".equals(printIfNecessary)) {
					log.info("Aligning Peptides " + printIfNecessary);
				}

				for (int j = i + 1; j < peptideList.size(); j++) {
					totalAligments++;
					final QuantifiedPeptideInterface pep2 = peptideList.get(j);
					final NWResult alignment = NWAlign.needlemanWunsch(pep1.getSequence(), pep2.getSequence(), -11, -1);
					maxAlignmentScore = maxAlignmentScore < alignment.getFinalAlignmentScore()
							? alignment.getFinalAlignmentScore() : maxAlignmentScore;
					maxSeqIdentity = maxSeqIdentity < alignment.getSequenceIdentity() ? alignment.getSequenceIdentity()
							: maxSeqIdentity;
					maxConsecutiveIdenticalAlignment = maxConsecutiveIdenticalAlignment < alignment
							.getMaxConsecutiveIdenticalAlignment() ? alignment.getMaxConsecutiveIdenticalAlignment()
									: maxConsecutiveIdenticalAlignment;

					if ((alignment.getFinalAlignmentScore() >= finalAlignmentScore
							&& alignment.getSequenceIdentity() >= sequenceIdentity)
							&& alignment.getMaxConsecutiveIdenticalAlignment() >= minConsecutiveIdenticalAlignment) {
						someAlignmentPassThresholds = true;
						// store aligment
						ret.addAlignment(new AlignedPeptides(alignment, pep1, pep2));
						// print pep1 and pep2
						alignmentLogFile.append(pep1.getSequence() + "\t");
						alignmentLogFile.append(pep2.getSequence() + "\t");
						Double ratio1 = null;
						if (pep1 instanceof IsobaricQuantifiedPeptide) {
							ratio1 = ((IsobaricQuantifiedPeptide) pep1).getIonCountRatio(cond1, cond2)
									.getLog2Ratio(cond1, cond2);
						} else if (pep1 instanceof QuantifiedPeptide) {
							ratio1 = ((QuantifiedPeptide) pep1).getConsensusRatio(cond1, cond2).getLog2Ratio(cond1,
									cond2);
						}
						if (ratio1 == Double.NEGATIVE_INFINITY) {
							alignmentLogFile.append("NEG_INF" + "\t");
						} else if (ratio1 == Double.POSITIVE_INFINITY) {
							alignmentLogFile.append("POS_INF" + "\t");
						} else {
							alignmentLogFile.append(ratio1 + "\t");
						}
						Double ratio2 = null;
						if (pep2 instanceof IsobaricQuantifiedPeptide) {
							ratio2 = ((IsobaricQuantifiedPeptide) pep2).getIonCountRatio(cond1, cond2)
									.getLog2Ratio(cond1, cond2);
						} else if (pep2 instanceof QuantifiedPeptide) {
							ratio2 = ((QuantifiedPeptide) pep2).getConsensusRatio(cond1, cond2).getLog2Ratio(cond1,
									cond2);
						}
						if (ratio2 == Double.NEGATIVE_INFINITY) {
							alignmentLogFile.append("NEG_INF" + "\t");
						} else if (ratio2 == Double.POSITIVE_INFINITY) {
							alignmentLogFile.append("POS_INF" + "\t");
						} else {
							alignmentLogFile.append(ratio2 + "\t");
						}

						alignmentLogFile.append(alignment.getSequenceIdentity() + "\t");
						alignmentLogFile.append("1" + "\n");
						alignmentLogFile.flush();

						/*
						 * // prints p1 and pep2
						 * System.out.print(quantifiedProtein.getAccession() +
						 * "\t"); System.out.print(pep2.getSequence() + "\t");
						 * System.out.print("1" + "\t"); System.out.print("0" +
						 * "\t"); log.info("1");
						 */

						// }
					}

				}
			}
			if (!someAlignmentPassThresholds) {
				log.info("None of the aligments passed the thresholds: minAlignmentScore=" + finalAlignmentScore
						+ ", minSeqIdentity=" + sequenceIdentity + ", minConsecutiveAlignment="
						+ minConsecutiveIdenticalAlignment);
				log.info("Better values were: minAlignmentScore=" + maxAlignmentScore + ", minSeqIdentity="
						+ maxSeqIdentity + ", minConsecutiveAlignment=" + maxConsecutiveIdenticalAlignment);
			}
			return ret;
		} finally {
			if (ret.getNumAligments() > 0) {
				log.info(ret.getNumAligments() + " aligments passed the threshold out of " + totalAligments);
			}
			if (alignmentLogFile != null) {
				alignmentLogFile.close();
			}
		}

	}

	// merges the clusters if there is a similar peptide pair between the two
	// clusters
	public static ProteinCluster mergeClusters(ProteinCluster cluster, ProteinCluster cluster2) {
		for (final QuantifiedPeptideInterface peptide : cluster2.getPeptideSet()) {
			cluster.addIndividualQuantifiedPeptide(peptide);
		}

		for (final QuantifiedProteinInterface protein : cluster2.getProteinSet()) {
			cluster.addIndividualQuantifiedProtein(protein);
		}
		return cluster;
	}

	public static CensusChroParser getCensusChroParser(File fastaFile, File inputFileFolder, String[] fileNames,
			Map<QuantCondition, QuantificationLabel> labelsByConditions, QuantificationLabel numeratorLabel,
			QuantificationLabel denominatorLabel, char[] enzymeArray, int missedCleavages, boolean semiCleavage,
			File uniprotReleasesFolder, String uniprotVersion, String decoyRegexp, boolean ignoreNotFoundPeptidesInDB,
			boolean distinguishModifiedPeptides, String peptideFilterRegexp, char[] quantifiedAAs,
			Boolean lookForProteoforms) throws FileNotFoundException {

		final List<Map<QuantCondition, QuantificationLabel>> list = new ArrayList<Map<QuantCondition, QuantificationLabel>>();
		for (int i = 0; i < fileNames.length; i++) {
			list.add(labelsByConditions);
		}
		return getCensusChroParser(fastaFile, inputFileFolder, fileNames, list, numeratorLabel, denominatorLabel,
				enzymeArray, missedCleavages, semiCleavage, uniprotReleasesFolder, uniprotVersion, decoyRegexp,
				ignoreNotFoundPeptidesInDB, distinguishModifiedPeptides, peptideFilterRegexp, quantifiedAAs,
				lookForProteoforms);
	}

	public static CensusChroParser getCensusChroParser(File fastaFile, File inputFilefolder, String[] fileNames,
			List<Map<QuantCondition, QuantificationLabel>> labelsByConditions, QuantificationLabel numeratorLabel,
			QuantificationLabel denominatorLabel, char[] enzymeArray, int missedCleavages, boolean semiCleavage,
			File uniprotReleasesFolder, String uniprotVersion, String decoyRegexp, boolean ignoreNotFoundPeptidesInDB,
			boolean distinguishModifiedPeptides, String peptideFilterRegexp, char[] quantifiedAAs,
			Boolean lookForProteoforms) throws FileNotFoundException {

		final List<RemoteSSHFileReference> xmlFiles = new ArrayList<RemoteSSHFileReference>();
		// Set parser (6 files) to peptides
		for (final String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.add(new RemoteSSHFileReference(inputXmlFile));
		}
		final String fileNamesKey = getFileNamesKey(fileNames);
		if (quantParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return (CensusChroParser) quantParsersByFileNamesKey.get(fileNamesKey);
		}

		final CensusChroParser parser = new CensusChroParser(xmlFiles, labelsByConditions, numeratorLabel,
				denominatorLabel);
		try {
			parser.setRetrieveFastaIsoforms(ProteinClusterQuantParameters.getInstance().isCollapseBySites());
			parser.addIonExclusion(IonSerieType.B, 1);
			parser.addIonExclusion(IonSerieType.Y, 1);
			parser.setDecoyPattern(decoyRegexp);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
			if (quantifiedAAs != null) {
				for (final char quantifiedAA : quantifiedAAs) {
					parser.addQuantifiedAA(quantifiedAA);
				}
			}

			final DBIndexInterface fastaDBIndex = getFastaDBIndex(fastaFile, enzymeArray, missedCleavages, semiCleavage,
					peptideFilterRegexp, uniprotReleasesFolder, lookForProteoforms);
			parser.setDbIndex(fastaDBIndex);
			parser.enableProteinMergingBySecondaryAccessions(
					getUniprotProteinLocalRetrieverByFolder(uniprotReleasesFolder), uniprotVersion);

			return parser;
		} finally {
			addQuantParserToStaticMap(fileNamesKey, parser);
		}
	}

	private static DTASelectParser getDTASelectParser(File fastaFile, File inputFilefolder, String[] fileNames,
			char[] enzymeArray, int missedCleavages, boolean semiCleavage, File uniprotReleasesFolder,
			String uniprotVersion, String decoyRegexp, boolean ignoreNotFoundPeptidesInDB, String peptideFilterRegexp,
			Boolean lookForProteoforms) throws FileNotFoundException {

		// Set parser (6 files) to peptides
		final Map<String, RemoteSSHFileReference> xmlFiles = new THashMap<String, RemoteSSHFileReference>();
		for (final String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.put(fileName, new RemoteSSHFileReference(inputXmlFile));
		}
		final String fileNamesKey = getFileNamesKey(fileNames);
		if (dtaSelectParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return dtaSelectParsersByFileNamesKey.get(fileNamesKey);
		}

		final DTASelectParser parser = new DTASelectParser(xmlFiles);
		try {

			parser.setDecoyPattern(decoyRegexp);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
			final DBIndexInterface fastaDBIndex = getFastaDBIndex(fastaFile, enzymeArray, missedCleavages, semiCleavage,
					peptideFilterRegexp, uniprotReleasesFolder, lookForProteoforms);
			parser.setDbIndex(fastaDBIndex);
			parser.enableProteinMergingBySecondaryAccessions(
					getUniprotProteinLocalRetrieverByFolder(uniprotReleasesFolder), uniprotVersion);

			return parser;
		} finally {
			addDTASelectParserToStaticMap(fileNamesKey, parser);
		}
	}

	private static String getFileNamesKey(String[] xmlFiles) {
		final StringBuilder sb = new StringBuilder();
		for (final String key : xmlFiles) {
			sb.append(key);
		}
		return sb.toString();
	}

	private static CensusOutParser getCensusOutParser(File fastaFile, File inputFilefolder, String[] fileNames,
			List<Map<QuantCondition, QuantificationLabel>> labelsByConditions, QuantificationLabel numeratorLabel,
			QuantificationLabel denominatorLabel, char[] enzymeArray, int missedCleavages, boolean semiCleavage,
			File uniprotReleasesFolder, String uniprotVersion, String decoyRegexp, boolean ignoreNotFoundPeptidesInDB,
			boolean onlyOneSpectrumPerChromatographicPeakAndPerSaltStep, boolean skipSingletons,
			boolean distinguishModifiedPeptides, String peptideFilterRegexp, char[] quantifiedAAs,
			Boolean lookForProteoforms) throws FileNotFoundException {

		// Set parser (6 files) to peptides
		final List<RemoteSSHFileReference> xmlFiles = new ArrayList<RemoteSSHFileReference>();

		for (final String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.add(new RemoteSSHFileReference(inputXmlFile));
		}
		final String fileNamesKey = getFileNamesKey(fileNames);
		if (quantParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return (CensusOutParser) quantParsersByFileNamesKey.get(fileNamesKey);
		}
		final CensusOutParser parser = new CensusOutParser(xmlFiles, labelsByConditions, numeratorLabel,
				denominatorLabel);
		try {
			final boolean collapseBySites = ProteinClusterQuantParameters.getInstance().isCollapseBySites();
			parser.setRetrieveFastaIsoforms(collapseBySites);
			parser.setDecoyPattern(decoyRegexp);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
			parser.setOnlyOneSpectrumPerChromatographicPeakAndPerSaltStep(
					onlyOneSpectrumPerChromatographicPeakAndPerSaltStep);
			parser.setSkipSingletons(skipSingletons);
			final DBIndexInterface fastaDBIndex = getFastaDBIndex(fastaFile, enzymeArray, missedCleavages, semiCleavage,
					peptideFilterRegexp, uniprotReleasesFolder, lookForProteoforms);
			parser.setDbIndex(fastaDBIndex);
			parser.enableProteinMergingBySecondaryAccessions(
					getUniprotProteinLocalRetrieverByFolder(uniprotReleasesFolder), uniprotVersion);
			if (quantifiedAAs != null) {
				for (final char quantifiedAA : quantifiedAAs) {
					parser.addQuantifiedAA(quantifiedAA);
				}
			}
			return parser;
		} finally {
			addQuantParserToStaticMap(fileNamesKey, parser);

		}
	}

	private static SeparatedValuesParser getSeparatedValuesParser(File fastaFile, File inputFilefolder,
			String[] fileNames, String separator, List<Map<QuantCondition, QuantificationLabel>> labelsByConditions,
			QuantificationLabel numeratorLabel, QuantificationLabel denominatorLabel, char[] enzymeArray,
			int missedCleavages, boolean semiCleavage, File uniprotReleasesFolder, String uniprotVersion,
			String decoyRegexp, boolean ignoreNotFoundPeptidesInDB, boolean distinguishModifiedPeptides,
			String peptideFilterRegexp, char[] quantifiedAAs, Boolean lookForProteoforms) throws FileNotFoundException {
		// Set parser (6 files) to peptides
		final List<RemoteSSHFileReference> xmlFiles = new ArrayList<RemoteSSHFileReference>();

		for (final String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.add(new RemoteSSHFileReference(inputXmlFile));
		}
		final String fileNamesKey = getFileNamesKey(fileNames);
		if (quantParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return (SeparatedValuesParser) quantParsersByFileNamesKey.get(fileNamesKey);
		}
		final SeparatedValuesParser parser = new SeparatedValuesParser(xmlFiles, separator, labelsByConditions,
				numeratorLabel, denominatorLabel);
		try {
			parser.setRetrieveFastaIsoforms(ProteinClusterQuantParameters.getInstance().isCollapseBySites());
			parser.setDecoyPattern(decoyRegexp);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
			final DBIndexInterface fastaDBIndex = getFastaDBIndex(fastaFile, enzymeArray, missedCleavages, semiCleavage,
					peptideFilterRegexp, uniprotReleasesFolder, lookForProteoforms);
			parser.setDbIndex(fastaDBIndex);
			parser.enableProteinMergingBySecondaryAccessions(
					getUniprotProteinLocalRetrieverByFolder(uniprotReleasesFolder), uniprotVersion);
			if (quantifiedAAs != null) {
				for (final char quantifiedAA : quantifiedAAs) {
					parser.addQuantifiedAA(quantifiedAA);
				}
			}
			return parser;
		} finally {
			addQuantParserToStaticMap(fileNamesKey, parser);
		}
	}

	private static void addQuantParserToStaticMap(String fileNamesKey, QuantParser parser) {
		quantParsersByFileNamesKey.put(fileNamesKey, parser);
		log.info(quantParsersByFileNamesKey.size() + " parsers stored.");
	}

	private static void addDTASelectParserToStaticMap(String fileNamesKey, DTASelectParser parser) {
		dtaSelectParsersByFileNamesKey.put(fileNamesKey, parser);
		log.info(dtaSelectParsersByFileNamesKey.size() + " parsers stored.");
	}

	public static DBIndexInterface getFastaDBIndex(File fastaFile, char[] enzymeArray, int missedCleavages,
			boolean semicleavage, String peptideFilterRegexp, File uniprotReleasesFolder, Boolean lookForProteoforms) {
		if (fastaFile != null) {

			final DBIndexSearchParams defaultDBIndexParams = DBIndexInterface.getDefaultDBIndexParams(fastaFile);
			final String fastaIndexKey = IndexUtil.createFullIndexFileName(defaultDBIndexParams);
			if (indexByFastaIndexKey.containsKey(fastaIndexKey)) {
				return indexByFastaIndexKey.get(fastaIndexKey);
			}
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeArr(enzymeArray, missedCleavages, semicleavage);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeOffset(0);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeNocutResidues("");
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setH2OPlusProtonAdded(true);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setUniprotReleasesFolder(uniprotReleasesFolder);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setLookProteoforms(lookForProteoforms);
			// if looking for proteoforms, not use in memory
			final boolean inMemoryIndex = !lookForProteoforms;
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setInMemoryIndex(inMemoryIndex);
			if (peptideFilterRegexp != null) {
				((DBIndexSearchParamsImpl) defaultDBIndexParams)
						.setPeptideFilter(new PeptideFilterBySequence(peptideFilterRegexp));
			}

			final DBIndexInterface dbIndex = new DBIndexInterface(defaultDBIndexParams);
			indexByFastaIndexKey.put(fastaIndexKey, dbIndex);
			return dbIndex;
		}
		return null;
	}

	private static CensusChroParser getCensusChroParserUsingMongoDBIndex(String mongoDBURI, String mongoMassDBName,
			String mongoSeqDBName, String mongoProtDBName, File inputFilefolder, String[] fileNames,
			List<Map<QuantCondition, QuantificationLabel>> labelsByConditions, QuantificationLabel numeratorLabel,
			QuantificationLabel denominatorLabel, File uniprotReleasesFolder, String uniprotVersion, String decoyRegexp,
			boolean ignoreNotFoundPeptidesInDB, boolean distinguishModifiedPeptides, String peptideFilterRegexp,
			char[] quantifiedAAs) throws FileNotFoundException {
		// Set parser (6 files) to peptides
		final List<RemoteSSHFileReference> xmlFiles = new ArrayList<RemoteSSHFileReference>();

		for (final String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.add(new RemoteSSHFileReference(inputXmlFile));
		}
		final String fileNamesKey = getFileNamesKey(fileNames);
		if (quantParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return (CensusChroParser) quantParsersByFileNamesKey.get(fileNamesKey);
		}
		final CensusChroParser parser = new CensusChroParser(xmlFiles, labelsByConditions, numeratorLabel,
				denominatorLabel);
		try {
			parser.setRetrieveFastaIsoforms(ProteinClusterQuantParameters.getInstance().isCollapseBySites());
			parser.addIonExclusion(IonSerieType.B, 1);
			parser.addIonExclusion(IonSerieType.Y, 1);
			parser.setDecoyPattern(decoyRegexp);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
			if (quantifiedAAs != null) {
				for (final char quantifiedAA : quantifiedAAs) {
					parser.addQuantifiedAA(quantifiedAA);
				}
			}

			final DBIndexInterface mongoDBIndex = getMongoDBIndex(mongoDBURI, mongoMassDBName, mongoSeqDBName,
					mongoProtDBName, peptideFilterRegexp);
			parser.setDbIndex(mongoDBIndex);
			parser.enableProteinMergingBySecondaryAccessions(
					getUniprotProteinLocalRetrieverByFolder(uniprotReleasesFolder), uniprotVersion);

			return parser;
		} finally {
			addQuantParserToStaticMap(fileNamesKey, parser);

		}
	}

	public static UniprotProteinLocalRetriever getUniprotProteinLocalRetrieverByFolder(File uniprotReleasesFolder) {
		if (uniprotReleasesFolder == null) {
			return null;
		}
		if (!uplrMap.containsKey(uniprotReleasesFolder.getAbsolutePath())) {
			final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
			uplrMap.put(uniprotReleasesFolder.getAbsolutePath(), uplr);
		}
		return uplrMap.get(uniprotReleasesFolder.getAbsolutePath());
	}

	private static CensusOutParser getCensusOutParserUsingMongoDBIndex(String mongoDBURI, String mongoMassDBName,
			String mongoSeqDBName, String mongoProtDBName, File inputFilefolder, String[] fileNames,
			List<Map<QuantCondition, QuantificationLabel>> labelsByConditions, QuantificationLabel numeratorLabel,
			QuantificationLabel denominatorLabel, File uniprotReleasesFolder, String uniprotVersion, String decoyRegexp,
			boolean ignoreNotFoundPeptidesInDB, boolean onlyOneSpectrumPerChromatographicPeakAndPerSaltStep,
			boolean skipSingletons, boolean distinguishModifiedPeptides, String peptideFilterRegexp,
			char[] quantifiedAAs) throws FileNotFoundException {
		// Set parser (6 files) to peptides
		final List<RemoteSSHFileReference> xmlFiles = new ArrayList<RemoteSSHFileReference>();

		for (final String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.add(new RemoteSSHFileReference(inputXmlFile));
		}
		final String fileNamesKey = getFileNamesKey(fileNames);
		if (quantParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return (CensusOutParser) quantParsersByFileNamesKey.get(fileNamesKey);
		}
		final CensusOutParser parser = new CensusOutParser(xmlFiles, labelsByConditions, numeratorLabel,
				denominatorLabel);
		try {
			parser.setRetrieveFastaIsoforms(ProteinClusterQuantParameters.getInstance().isCollapseBySites());
			parser.setDecoyPattern(decoyRegexp);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
			parser.setOnlyOneSpectrumPerChromatographicPeakAndPerSaltStep(
					onlyOneSpectrumPerChromatographicPeakAndPerSaltStep);
			parser.setSkipSingletons(skipSingletons);
			final DBIndexInterface dbIndex = getMongoDBIndex(mongoDBURI, mongoMassDBName, mongoSeqDBName,
					mongoProtDBName, peptideFilterRegexp);
			parser.setDbIndex(dbIndex);
			parser.enableProteinMergingBySecondaryAccessions(
					getUniprotProteinLocalRetrieverByFolder(uniprotReleasesFolder), uniprotVersion);
			if (quantifiedAAs != null) {
				for (final char quantifiedAA : quantifiedAAs) {
					parser.addQuantifiedAA(quantifiedAA);
				}
			}
			return parser;
		} finally {
			addQuantParserToStaticMap(fileNamesKey, parser);

		}
	}

	private static DTASelectParser getDTASelectParserUsingMongoDBIndex(String mongoDBURI, String mongoMassDBName,
			String mongoSeqDBName, String mongoProtDBName, File inputFilefolder, String[] fileNames,
			File uniprotReleasesFolder, String uniprotVersion, String decoyRegexp, boolean ignoreNotFoundPeptidesInDB,
			String peptideFilterRegexp) throws FileNotFoundException {
		// Set parser (6 files) to peptides
		final Map<String, RemoteSSHFileReference> xmlFiles = new THashMap<String, RemoteSSHFileReference>();

		for (final String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.put(fileName, new RemoteSSHFileReference(inputXmlFile));
		}
		final String fileNamesKey = getFileNamesKey(fileNames);
		if (dtaSelectParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return dtaSelectParsersByFileNamesKey.get(fileNamesKey);
		}
		final DTASelectParser parser = new DTASelectParser(xmlFiles);
		try {
			parser.setDecoyPattern(decoyRegexp);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
			final DBIndexInterface dbIndex = getMongoDBIndex(mongoDBURI, mongoMassDBName, mongoSeqDBName,
					mongoProtDBName, peptideFilterRegexp);
			parser.setDbIndex(dbIndex);
			parser.enableProteinMergingBySecondaryAccessions(
					getUniprotProteinLocalRetrieverByFolder(uniprotReleasesFolder), uniprotVersion);
			return parser;
		} finally {
			addDTASelectParserToStaticMap(fileNamesKey, parser);

		}
	}

	private static SeparatedValuesParser getSeparatedValuesParserUsingMongoDBIndex(String mongoDBURI,
			String mongoMassDBName, String mongoSeqDBName, String mongoProtDBName, File inputFilefolder,
			String[] fileNames, String separator, List<Map<QuantCondition, QuantificationLabel>> labelsByConditions,
			QuantificationLabel numeratorLabel, QuantificationLabel denominatorLabel, File uniprotReleasesFolder,
			String uniprotVersion, String decoyRegexp, boolean ignoreNotFoundPeptidesInDB,
			boolean distinguishModifiedPeptides, String peptideFilterRegexp, char[] quantifiedAAs)
			throws FileNotFoundException {
		// Set parser (6 files) to peptides
		final List<RemoteSSHFileReference> xmlFiles = new ArrayList<RemoteSSHFileReference>();

		for (final String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.add(new RemoteSSHFileReference(inputXmlFile));
		}
		final String fileNamesKey = getFileNamesKey(fileNames);
		if (quantParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return (SeparatedValuesParser) quantParsersByFileNamesKey.get(fileNamesKey);
		}
		final SeparatedValuesParser parser = new SeparatedValuesParser(xmlFiles, separator, labelsByConditions,
				numeratorLabel, denominatorLabel);
		try {
			parser.setRetrieveFastaIsoforms(ProteinClusterQuantParameters.getInstance().isCollapseBySites());
			parser.setDecoyPattern(decoyRegexp);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
			final DBIndexInterface dbIndex = getMongoDBIndex(mongoDBURI, mongoMassDBName, mongoSeqDBName,
					mongoProtDBName, peptideFilterRegexp);
			parser.setDbIndex(dbIndex);
			parser.enableProteinMergingBySecondaryAccessions(
					getUniprotProteinLocalRetrieverByFolder(uniprotReleasesFolder), uniprotVersion);
			if (quantifiedAAs != null) {
				for (final char quantifiedAA : quantifiedAAs) {
					parser.addQuantifiedAA(quantifiedAA);
				}
			}
			return parser;
		} finally {
			addQuantParserToStaticMap(fileNamesKey, parser);
		}
	}

	private static DBIndexInterface getMongoDBIndex(String mongoDBURI, String mongoMassDBName, String mongoSeqDBName,
			String mongoProtDBName, String peptideFilterRegexp) {

		if (mongoDBURI != null) {
			log.info("Using mongoDB at: " + mongoDBURI + " massDBName:" + mongoMassDBName + " seqDBName:"
					+ mongoSeqDBName + " protDBName:" + mongoProtDBName);

			PeptideFilterBySequence peptideFilter = null;
			if (peptideFilterRegexp != null) {
				peptideFilter = new PeptideFilterBySequence(peptideFilterRegexp);
			}
			final DBIndexSearchParamsImpl params = new DBIndexSearchParamsImpl(mongoDBURI, mongoMassDBName,
					mongoSeqDBName, mongoProtDBName, peptideFilter);
			final DBIndexInterface dbIndex = new DBIndexInterface(params);
			return dbIndex;
		}
		return null;
	}

	public static boolean shareAtLeastOnePeptideNode(PCQProteinNode proteinNode1, PCQProteinNode proteinNode2,
			boolean skipDiscardedPeptideNodes) {

		final Set<PCQPeptideNode> peptideNodes1 = proteinNode1.getPeptideNodes();
		for (final PCQPeptideNode peptideNode : peptideNodes1) {
			if (skipDiscardedPeptideNodes && peptideNode.isDiscarded()) {
				continue;
			}
			if (proteinNode2.getPeptideNodes().contains(peptideNode)) {
				return true;
			}
		}
		return false;
	}

	public static boolean shareAtLeastOnePeptideBySimilarity(PCQProteinNode proteinNode1, PCQProteinNode proteinNode2,
			AlignmentSet peptideAlignments, boolean skipDiscardedPeptideNodes) {
		final Set<PCQPeptideNode> peptideNodes1 = proteinNode1.getPeptideNodes();

		for (final PCQPeptideNode peptideNode1 : peptideNodes1) {
			if (skipDiscardedPeptideNodes && peptideNode1.isDiscarded()) {
				continue;
			}
			final QuantifiedPeptideInterface peptide1 = peptideNode1.getItemsInNode().iterator().next();
			if (peptideAlignments != null) {
				final Set<AlignedPeptides> alignments = peptideAlignments.getAlignmentsForPeptide(peptide1);

				for (final AlignedPeptides alignment : alignments) {
					final QuantifiedPeptideInterface peptide2 = alignment.getPeptideAligned(peptide1);
					if (peptide2 != null) {

						if (proteinNodeContainsPeptideSequence(proteinNode2, peptide2.getSequence())) {
							return true;
						}
					}
				}

			}
		}
		return false;
	}

	private static Boolean proteinNodeContainsPeptideSequence(PCQProteinNode proteinNode, String peptideSeq) {
		final Set<QuantifiedPeptideInterface> peptides = proteinNode.getQuantifiedPeptides();
		for (final QuantifiedPeptideInterface peptide : peptides) {
			if (peptide.getSequence().equals(peptideSeq)) {
				return true;
			}
		}
		return false;
	}

	public static Set<PCQPeptideNode> getSharedPeptideNodeSet(PCQProteinNode proteinNode1, PCQProteinNode proteinNode2,
			boolean onlysharedByThisToProteins, boolean skipDiscarded) {
		final Set<PCQPeptideNode> peptideNodes1 = proteinNode1.getPeptideNodes();
		final Set<PCQPeptideNode> peptideNodes2 = proteinNode2.getPeptideNodes();
		final Set<PCQPeptideNode> ret = new THashSet<PCQPeptideNode>();
		for (final PCQPeptideNode peptideNode1 : peptideNodes1) {
			if (skipDiscarded && peptideNode1.isDiscarded()) {
				continue;
			}
			if (peptideNodes2.contains(peptideNode1)) {
				if (!onlysharedByThisToProteins) {
					ret.add(peptideNode1);
				} else {
					// only is valid if the peptide node only contains 2
					// protein nodes
					if (peptideNode1.getProteinNodes().size() == 2) {
						ret.add(peptideNode1);
					}
				}
			}
		}
		return ret;
	}

	/**
	 * Get the list of individual ConsensusRatios of the individual shared
	 * {@link QuantifiedPeptideInterface}s between two {@link PCQProteinNode}s
	 *
	 * @param proteinNode1
	 * @param proteinNode2
	 * @param cond1
	 * @param cond2
	 * @return
	 */
	public static List<Double> getSharedIndividualPeptideRatioValues(PCQProteinNode proteinNode1,
			PCQProteinNode proteinNode2, QuantCondition cond1, QuantCondition cond2,
			boolean onlySharedByTheseTwoProteins, boolean skipDiscarded) {
		final Set<PCQPeptideNode> sharedPeptideNodes = getSharedPeptideNodeSet(proteinNode1, proteinNode2,
				onlySharedByTheseTwoProteins, skipDiscarded);

		final List<Double> ratioValues = new ArrayList<Double>();

		for (final PCQPeptideNode peptideNode : sharedPeptideNodes) {
			if (skipDiscarded && peptideNode.isDiscarded()) {
				continue;
			}
			for (final QuantifiedPeptideInterface peptide : peptideNode.getQuantifiedPeptides()) {
				final QuantRatio consensusRatio = peptide.getConsensusRatio(cond1, cond2);
				if (consensusRatio != null) {
					ratioValues.add(consensusRatio.getLog2Ratio(cond1, cond2));
				}
			}
		}
		return ratioValues;
	}

	/**
	 * Get the list of individual ConsensusRatios of the individual shared
	 * {@link QuantifiedPSMInterface}s between two {@link PCQProteinNode}s
	 *
	 * @param proteinNode1
	 * @param proteinNode2
	 * @param cond1
	 * @param cond2
	 * @return
	 */
	public static List<Double> getSharedIndividualPSMRatioValues(PCQProteinNode proteinNode1,
			PCQProteinNode proteinNode2, QuantCondition cond1, QuantCondition cond2,
			boolean onlySharedByTheseTwoProteins, boolean skipDiscarded) {
		final Set<PCQPeptideNode> sharedPeptideNodes = getSharedPeptideNodeSet(proteinNode1, proteinNode2,
				onlySharedByTheseTwoProteins, skipDiscarded);

		final List<Double> ratioValues = new ArrayList<Double>();

		for (final PCQPeptideNode peptideNode : sharedPeptideNodes) {
			if (skipDiscarded && peptideNode.isDiscarded()) {
				continue;
			}
			for (final QuantifiedPSMInterface psm : peptideNode.getQuantifiedPSMs()) {
				final QuantRatio consensusRatio = psm.getConsensusRatio(cond1, cond2);
				if (consensusRatio != null) {
					ratioValues.add(consensusRatio.getLog2Ratio(cond1, cond2));
				}
			}
		}
		return ratioValues;
	}

	/**
	 * Get the shared peptides between protein1 and protein2. Each peptide will
	 * be mapped to the proteins that are mapped.
	 *
	 * @param proteinNode1
	 * @param proteinNode2
	 * @param onlySharedByTheseTwoProteins
	 * @return
	 */
	public static Map<String, Set<PCQPeptideNode>> getSharedPeptideNodesByProteinNode(PCQProteinNode proteinNode1,
			PCQProteinNode proteinNode2, boolean onlySharedByTheseTwoProteins, boolean skipDiscarded) {

		final Set<PCQPeptideNode> peptidesNodes1 = proteinNode1.getPeptideNodes();
		final Set<PCQPeptideNode> peptidesNodes2 = proteinNode2.getPeptideNodes();
		final Set<PCQPeptideNode> totalPeptideNodes = new THashSet<PCQPeptideNode>();
		totalPeptideNodes.addAll(peptidesNodes1);
		totalPeptideNodes.addAll(peptidesNodes2);
		final Map<String, Set<PCQPeptideNode>> map = new THashMap<String, Set<PCQPeptideNode>>();
		for (final PCQPeptideNode peptideNode : totalPeptideNodes) {
			if (skipDiscarded && peptideNode.isDiscarded()) {
				continue;
			}
			final Set<PCQProteinNode> proteinNodes = peptideNode.getProteinNodes();
			if (proteinNodes.contains(proteinNode1) && proteinNodes.contains(proteinNode2)) {
				boolean include = false;
				if (onlySharedByTheseTwoProteins) {
					if (proteinNodes.size() == 2) {
						include = true;
					}
				} else {
					include = true;
				}
				if (include) {
					// peptide shared by protein1 and protein2
					for (final PCQProteinNode proteinNode : proteinNodes) {
						if (skipDiscarded && proteinNode.isDiscarded()) {
							continue;
						}
						final String proteinAccKey = proteinNode.getKey();
						if (map.containsKey(proteinAccKey)) {
							map.get(proteinAccKey).add(peptideNode);
						} else {
							final Set<PCQPeptideNode> set = new THashSet<PCQPeptideNode>();
							set.add(peptideNode);
							map.put(proteinAccKey, set);
						}
					}

				}
			}
		}
		return map;
	}

	/**
	 * Get the shared peptides between proteinNode1 and proteinNode2. It is
	 * returned in a Map, in which each entry is a set of peptideNodes that
	 * share the same proteinNodes.<br>
	 * For example, one set will be the peptideNodes shared ONLY by proteinNode1
	 * and proteinNode2. Other set would be the peptideNodes shared by
	 * proteinNode1, proteinNode2 and a third proteinNode3. Consequently,
	 * another set would be the peptideNodes shared by proteinNode1,
	 * proteinNode2, proteinNode3 and a fourth proteinNode4, for example. And so
	 * on.<br>
	 * For each peptideNode set the key of the map will be the protein
	 * accessions where they belong sorted alphabetically.
	 *
	 * @param proteinNode1
	 * @param proteinNode2
	 * @param onlySharedByTheseTwoProteins
	 * @return
	 */
	public static Map<String, Set<PCQPeptideNode>> getSharedPeptideNodesMap(PCQProteinNode proteinNode1,
			PCQProteinNode proteinNode2, boolean onlySharedByTheseTwoProteins, boolean skipDiscarded) {

		final Map<String, Set<PCQPeptideNode>> map = new THashMap<String, Set<PCQPeptideNode>>();
		if (proteinNode1 == null || proteinNode2 == null) {
			return map;
		}
		final Set<PCQPeptideNode> peptideNodes1 = proteinNode1.getPeptideNodes();
		final Set<PCQPeptideNode> peptideNodes2 = proteinNode2.getPeptideNodes();
		final Set<PCQPeptideNode> totalPeptideNodes = new THashSet<PCQPeptideNode>();
		totalPeptideNodes.addAll(peptideNodes1);
		totalPeptideNodes.addAll(peptideNodes2);

		for (final PCQPeptideNode peptideNode : totalPeptideNodes) {
			if (skipDiscarded && peptideNode.isDiscarded()) {
				continue;
			}
			final Set<PCQProteinNode> proteinNodes = peptideNode.getProteinNodes();
			if (proteinNodes.contains(proteinNode1) && proteinNodes.contains(proteinNode2)) {
				boolean include = false;
				if (onlySharedByTheseTwoProteins) {
					if (proteinNodes.size() == 2) {
						include = true;
					}
				} else {
					include = true;
				}
				if (include) {
					// peptide shared by protein1 and protein2
					final String proteinAccKey = PCQUtils.getProteinNodeAccessionString(proteinNodes);
					if (map.containsKey(proteinAccKey)) {
						map.get(proteinAccKey).add(peptideNode);
					} else {
						final Set<PCQPeptideNode> set = new THashSet<PCQPeptideNode>();
						set.add(peptideNode);
						map.put(proteinAccKey, set);
					}
				}
			}
		}
		return map;
	}

	/**
	 * Gets all the log2Values of the consensus ratios of the individual
	 * {@link QuantifiedPeptideInterface} of the uniques {@link PCQPeptideNode}s
	 * of proteinNode1 with respect to proteinNode2
	 *
	 * @param proteinNode1
	 * @param proteinNode2
	 * @param cond1
	 * @param cond2
	 * @param uniquePepOnly
	 * @param skipDiscarded
	 * @return
	 */
	public static List<Double> getUniqueIndividualPeptideRatioValues(PCQProteinNode proteinNode1,
			PCQProteinNode proteinNode2, QuantCondition cond1, QuantCondition cond2, boolean uniquePepOnly,
			boolean skipDiscarded) {
		final Set<PCQPeptideNode> peptideNodes1 = getUniquePeptideNodes(proteinNode1, proteinNode2, uniquePepOnly,
				skipDiscarded);
		final List<Double> ratios = new ArrayList<Double>();

		for (final PCQPeptideNode peptideNode1 : peptideNodes1) {
			if (skipDiscarded && peptideNode1.isDiscarded()) {
				continue;
			}
			for (final QuantifiedPeptideInterface peptide : peptideNode1.getQuantifiedPeptides()) {
				final QuantRatio consensusRatio = peptide.getConsensusRatio(cond1, cond2);
				if (consensusRatio != null) {
					ratios.add(consensusRatio.getLog2Ratio(cond1, cond2));
				}
			}

		}
		return ratios;
	}

	/**
	 * Gets all the log2Values of the consensus ratios of the individual
	 * {@link QuantifiedPSMInterface} of the uniques {@link PCQPeptideNode}s of
	 * proteinNode1 with respect to proteinNode2
	 *
	 * @param proteinNode1
	 * @param proteinNode2
	 * @param cond1
	 * @param cond2
	 * @param uniquePepOnly
	 * @param skipDiscarded
	 * @return
	 */
	public static List<Double> getUniqueIndividualPSMRatioValues(PCQProteinNode proteinNode1,
			PCQProteinNode proteinNode2, QuantCondition cond1, QuantCondition cond2, boolean uniquePepOnly,
			boolean skipDiscarded) {
		final Set<PCQPeptideNode> peptideNodes1 = getUniquePeptideNodes(proteinNode1, proteinNode2, uniquePepOnly,
				skipDiscarded);
		final List<Double> ratios = new ArrayList<Double>();

		for (final PCQPeptideNode peptideNode1 : peptideNodes1) {
			if (skipDiscarded && peptideNode1.isDiscarded()) {
				continue;
			}
			for (final QuantifiedPSMInterface psm : peptideNode1.getQuantifiedPSMs()) {
				final QuantRatio consensusRatio = psm.getConsensusRatio(cond1, cond2);
				if (consensusRatio != null) {
					ratios.add(consensusRatio.getLog2Ratio(cond1, cond2));
				}
			}

		}
		return ratios;
	}

	/**
	 * Checks to see if the list has both ratios and infinities in it. Returns
	 * true if mixed. Returns false if only one kind.
	 *
	 * @param ratios
	 * @return
	 */
	public static boolean hasRatiosAndINF(List<Double> ratios) {
		int counterINF = 0;
		int counterRatio = 0;
		for (final Double ratio : ratios) {
			if (ratio != null && !ratio.isNaN()) {
				if (Double.isInfinite(ratio)) {
					counterINF++;
				} else {
					counterRatio++;
				}
			} else {
				// means that one of the ratios is null!
			}
		}
		if (counterINF > 0 && counterRatio > 0) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Gets a list of all the infinities from the original list of ratios. If
	 * the ratio is null, it is ignored.
	 *
	 * @param ratios
	 * @return
	 */
	public static List<Double> getINFRatiosValues(List<Double> ratios) {
		final List<Double> INFList = new ArrayList<Double>();
		for (final Double ratio : ratios) {
			if (ratio != null && Double.compare(ratio, Double.NaN) != 0) {
				if (Double.isInfinite(ratio)) {
					INFList.add(ratio);
				}
			}
		}
		return INFList;
	}

	/**
	 * Gets a list of all the ratios (excludes INF) from the original list of
	 * ratios. If the ratio is null, it is ignored.
	 *
	 * @param ratios
	 * @return
	 */
	public static List<Double> getNonINFRatiosValues(List<Double> ratios) {
		final List<Double> RatioList = new ArrayList<Double>();
		for (final Double ratio : ratios) {
			if (ratio != null && Double.compare(ratio, Double.NaN) != 0) {
				if (!Double.isInfinite(ratio)) {
					RatioList.add(ratio);
				}
			}
		}
		return RatioList;
	}

	/**
	 * Checks to see if all the INF passed through are the same
	 *
	 * @param ratios
	 * @return Double, returns INF if all are INF, -INF if all are -INF, and
	 *         null if not the same.
	 */
	public static Double areAllINFValuesSame(List<Double> ratios) {
		int posCount = 0;
		int negCount = 0;
		for (final Double ratio : ratios) {
			if (Double.compare(Double.POSITIVE_INFINITY, ratio) == 0) {
				posCount++;
			}
			if (Double.compare(Double.NEGATIVE_INFINITY, ratio) == 0) {
				negCount++;
			}
		}
		if (posCount > 0 && negCount > 0) {
			return null;
		}
		if (posCount > 0 && negCount == 0) {
			return Double.POSITIVE_INFINITY;
		}

		if (posCount == 0 && negCount > 0) {
			return Double.NEGATIVE_INFINITY;
		}
		return null;
	}

	/**
	 * Checks to see if all the INF passed through are the same
	 *
	 * @param ratios
	 * @return Double, returns INF if all are INF, -INF if all are -INF, and
	 *         null if not the same.
	 */
	public static Double areAllINFSame(List<IonCountRatio> ratios, QuantCondition cond1, QuantCondition cond2) {
		int posCount = 0;
		int negCount = 0;
		for (final IonCountRatio ratio : ratios) {
			if (Double.compare(Double.POSITIVE_INFINITY, ratio.getLog2Ratio(cond1, cond2)) == 0) {
				posCount++;
			}
			if (Double.compare(Double.NEGATIVE_INFINITY, ratio.getLog2Ratio(cond1, cond2)) == 0) {
				negCount++;
			}
		}
		if (posCount > 0 && negCount > 0) {
			return null;
		}
		if (posCount > 0 && negCount == 0) {
			return Double.POSITIVE_INFINITY;
		}

		if (posCount == 0 && negCount > 0) {
			return Double.NEGATIVE_INFINITY;
		}
		return null;
	}

	/**
	 * Gets the number of ions (light or heavy) in the peptide
	 *
	 * @param peptideNode
	 * @return
	 */
	public static int getIonCountFromPeptideNode(PCQPeptideNode peptideNode, QuantCondition condition) {

		final Set<IsobaricQuantifiedPSM> psms = new THashSet<IsobaricQuantifiedPSM>();
		for (final QuantifiedPSMInterface psm : peptideNode.getQuantifiedPSMs()) {
			if (psm instanceof IsobaricQuantifiedPSM) {
				psms.add((IsobaricQuantifiedPSM) psm);
			}
		}
		return getIonCountFromPSMs(psms, condition);

	}

	private static int getIonCountFromPSMs(Collection<IsobaricQuantifiedPSM> psms, QuantCondition condition) {
		int total = 0;
		for (final IsobaricQuantifiedPSM quantifiedPSM : psms) {
			final Set<Ion> ions = quantifiedPSM.getIonsByCondition().get(condition);
			if (ions != null) {
				total = total + ions.size();
			}
		}
		return total;
	}

	public static List<QuantifiedPeptideInterface> getSortedPeptidesBySequence(
			Collection<QuantifiedPeptideInterface> peptides) {
		final List<QuantifiedPeptideInterface> ret = new ArrayList<QuantifiedPeptideInterface>();
		ret.addAll(peptides);
		Collections.sort(ret, new Comparator<QuantifiedPeptideInterface>() {

			@Override
			public int compare(QuantifiedPeptideInterface o1, QuantifiedPeptideInterface o2) {
				return o1.getFullSequence().compareTo(o2.getFullSequence());
			}
		});
		return ret;
	}

	public static List<PCQPeptideNode> getSortedPeptideNodesBySequence(Collection<PCQPeptideNode> peptides) {
		final List<PCQPeptideNode> ret = new ArrayList<PCQPeptideNode>();
		for (final PCQPeptideNode pcqPeptideNode : peptides) {
			if (pcqPeptideNode != null)
				ret.add(pcqPeptideNode);
		}
		Collections.sort(ret, new Comparator<PCQPeptideNode>() {

			@Override
			public int compare(PCQPeptideNode o1, PCQPeptideNode o2) {
				return o1.getFullSequence().compareTo(o2.getFullSequence());
			}
		});
		return ret;
	}

	public static List<PCQProteinNode> getSortedProteinNodesByAcc(Collection<PCQProteinNode> proteinNodesToSort) {
		final List<PCQProteinNode> list = new ArrayList<PCQProteinNode>();
		list.addAll(proteinNodesToSort);
		Collections.sort(list, new Comparator<PCQProteinNode>() {
			@Override
			public int compare(PCQProteinNode o1, PCQProteinNode o2) {
				return o1.getKey().compareTo(o2.getKey());
			}
		});
		return list;
	}

	/**
	 * Get a separated value list of protein accessions after sorting them
	 * alphabetically
	 *
	 * @param proteins
	 * @return
	 */
	public static String getAccessionString(Collection<QuantifiedProteinInterface> proteins) {
		final Set<String> set = new THashSet<String>();
		for (final QuantifiedProteinInterface protein : proteins) {
			set.add(protein.getAccession());
		}
		final List<String> list = new ArrayList<String>();
		list.addAll(set);
		Collections.sort(list);

		final StringBuilder sb = new StringBuilder();
		for (final String acc : list) {
			if (!"".equals(sb.toString()))
				sb.append(PROTEIN_ACC_SEPARATOR);
			sb.append(acc);
		}
		return sb.toString();
	}

	/**
	 * Get a separated value list of protein keys after sorting them
	 * alphabetically
	 *
	 * @param proteins
	 * @return
	 */
	public static String getKeyString(Collection<QuantifiedProteinInterface> proteins) {
		final Set<String> set = new THashSet<String>();
		for (final QuantifiedProteinInterface protein : proteins) {
			set.add(protein.getKey());
		}
		final List<String> list = new ArrayList<String>();
		list.addAll(set);
		Collections.sort(list);

		final StringBuilder sb = new StringBuilder();
		for (final String acc : list) {
			if (!"".equals(sb.toString()))
				sb.append(PROTEIN_ACC_SEPARATOR);
			sb.append(acc);
		}
		return sb.toString();
	}

	/**
	 * Get a CVS list of protein accessions after sorting them alphabetically
	 *
	 * @param proteinNodes
	 * @return
	 */
	public static String getProteinNodeAccessionString(Collection<PCQProteinNode> proteinNodes) {

		final StringBuilder sb = new StringBuilder();
		for (final PCQProteinNode proteinNode : PCQUtils.getSortedProteinNodesByAcc(proteinNodes)) {
			if (!"".equals(sb.toString()))
				sb.append(PROTEIN_ACC_SEPARATOR);
			sb.append(proteinNode.getKey());
		}
		return sb.toString();
	}

	/**
	 * Get a CVS list of peptide sequences after sorting them alphabetically by
	 * the sequence
	 *
	 * @param peptides
	 * @return
	 */
	public static String getPeptidesSequenceString(Collection<QuantifiedPeptideInterface> peptides) {

		final StringBuilder sb = new StringBuilder();
		for (final QuantifiedPeptideInterface peptide : PCQUtils.getSortedPeptidesBySequence(peptides)) {
			if (!"".equals(sb.toString()))
				sb.append("_");
			sb.append(peptide.getSequence());
		}
		return sb.toString();
	}

	/**
	 * Get a CVS list of peptide sequences after sorting them alphabetically by
	 * the sequence
	 *
	 * @param peptideNodes
	 * @return
	 */
	public static String getPeptideNodesSequenceString(Collection<PCQPeptideNode> peptideNodes) {

		final StringBuilder sb = new StringBuilder();
		for (final PCQPeptideNode peptideNode : PCQUtils.getSortedPeptideNodesBySequence(peptideNodes)) {
			if (!"".equals(sb.toString()))
				sb.append("_");
			sb.append(peptideNode.getFullSequence());
		}
		return sb.toString();
	}

	/**
	 * Gets a CSV list of protein descriptions, after sorting them by accession
	 * alphabetically
	 *
	 * @param proteins
	 * @params avoidRedundancy if true, not repeated names are returned
	 * @return
	 */
	public static String getDescriptionStringFromIndividualProteins(Collection<QuantifiedProteinInterface> proteins,
			boolean avoidRedundancy) {
		final StringBuilder sb = new StringBuilder();
		final Set<String> descriptions = new THashSet<String>();
		for (final QuantifiedProteinInterface protein : QuantUtils.getSortedQuantifiedProteinsByAcc(proteins)) {
			final String description = protein.getDescription();
			if (avoidRedundancy && !descriptions.contains(description)) {
				if (!"".equals(sb.toString())) {
					sb.append(PROTEIN_DESCRIPTION_SEPARATOR);
				}
				sb.append(description);
				descriptions.add(description);
			}
		}
		return sb.toString();
	}

	/**
	 * Formats a number, comparing it fo Infinities and printing POS_INF or
	 * NEG_INF, and comparing to Nan, printing DNQ (detected, not quantified)
	 *
	 * @param num
	 * @return
	 */
	public static String format(Double num) {
		if (num != null && !Double.isInfinite(num) && !Double.isNaN(num)) {
			return df.format(num);
		} else if (num != null && Double.isInfinite(num)) {
			if (Double.compare(num, Double.POSITIVE_INFINITY) == 0) {
				return "POS_INF";
			} else if (Double.compare(num, Double.NEGATIVE_INFINITY) == 0) {
				return "NEG_INF";
			}
		} else if (Double.isNaN(num)) {
			return "DNQ"; // means detected, not quantified
		}
		return String.valueOf(num);
	}

	/**
	 * Get peptide nodes licked to the protein node that only are linked to that
	 * one
	 *
	 * @param proteinNode
	 * @return
	 */
	public static Set<PCQPeptideNode> getUniquePeptideNodes(PCQProteinNode proteinNode, boolean skipDiscarded) {
		final Set<PCQPeptideNode> ret = new THashSet<PCQPeptideNode>();
		for (final PCQPeptideNode peptideNode : proteinNode.getPeptideNodes()) {
			if (skipDiscarded && peptideNode.isDiscarded()) {
				continue;
			}
			if (peptideNode.getProteinNodes().size() == 1) {
				ret.add(peptideNode);
			}
		}
		return ret;
	}

	// public static Double getAverageRatio(
	// Set<QuantifiedPeptideInterface> quantifiedPeptides, QuantCondition cond1,
	// QuantCondition cond2) {
	// List<Double> pepRatioList = new ArrayList<Double>();
	// for (QuantifiedPeptideInterface peps : quantifiedPeptides) {
	// Double ratio = peps.getCountRatio(cond1, cond2);
	// pepRatioList.add(ratio);
	// }
	// return ProteinPair.getConsensusRatio(pepRatioList);
	// }
	public static String getGeneNameString(Map<String, Entry> annotatedProteins, ProteinCluster cluster,
			Set<String> validTaxonomies, boolean onlyFirst, boolean skipDiscarded) {
		return getGeneNameString(annotatedProteins, getProteinMap(cluster.getProteinNodes(), skipDiscarded),
				validTaxonomies, onlyFirst);
	}

	public static String getGeneNameString(Map<String, Entry> annotatedProteins,
			Collection<PCQProteinNode> proteinNodes, Set<String> validTaxonomies, boolean onlyFirst,
			boolean skipDiscarded) {
		return getGeneNameString(annotatedProteins, getProteinMap(proteinNodes, skipDiscarded), validTaxonomies,
				onlyFirst);
	}

	public static String getGeneNameString(Map<String, Entry> annotatedProteins, PCQProteinNode proteinNode,
			Set<String> validTaxonomies, boolean onlyFirst, boolean skipDiscarded) {
		return getGeneNameString(annotatedProteins, getProteinMapByProteinAccession(proteinNode, skipDiscarded),
				validTaxonomies, onlyFirst);
	}

	private static Map<String, PCQProteinNode> getProteinMapByProteinAccession(PCQProteinNode protein,
			boolean skipDiscarded) {
		final Set<PCQProteinNode> set = new THashSet<PCQProteinNode>();
		set.add(protein);
		return getProteinMap(set, skipDiscarded);
	}

	private static Map<String, PCQProteinNode> getProteinMap(Collection<PCQProteinNode> proteinNodes,
			boolean skipDiscarded) {
		final Map<String, PCQProteinNode> map = new THashMap<String, PCQProteinNode>();
		for (final PCQProteinNode proteinNode : proteinNodes) {
			if (skipDiscarded && proteinNode.isDiscarded()) {
				continue;
			}
			final String accession = proteinNode.getKey();
			if (accession.contains(PROTEIN_ACC_SEPARATOR)) {
				final String[] split = accession.split(PROTEIN_ACC_SEPARATOR);
				for (final String string : split) {
					if (string.length() == 1)
						continue;
					map.put(string, proteinNode);
				}
			} else {
				map.put(accession, proteinNode);
			}
		}
		return map;
	}

	public static String getSpeciesString(Map<String, Entry> annotatedProteins, Set<PCQProteinNode> proteinNodeSet,
			Set<String> validTaxonomies, boolean skipDiscarded) {

		final Set<String> set = new THashSet<String>();
		for (final PCQProteinNode proteinNode : proteinNodeSet) {
			if (skipDiscarded && proteinNode.isDiscarded()) {
				continue;
			}
			final String rawAcc = proteinNode.getKey();
			final List<String> accs = new ArrayList<String>();
			if (rawAcc.contains(" ")) {
				final String[] split = rawAcc.split(" ");
				for (final String string : split) {
					accs.add(string);
				}
			} else {
				accs.add(rawAcc);
			}

			// int index = 0;
			for (final String acc : accs) {
				if (annotatedProteins != null && annotatedProteins.containsKey(acc)) {
					final String taxon = getTaxonomy(acc, annotatedProteins.get(acc).getOrganism());
					if (taxon != null) {
						boolean valid = false;
						if (validTaxonomies != null && !validTaxonomies.isEmpty()) {
							for (final String skipTaxonomy : validTaxonomies) {
								if (taxon.contains(skipTaxonomy)) {
									valid = true;
								}
							}
						} else {
							valid = true;
						}
						if (!valid) {
							continue;
						}
						set.add(taxon);
					}

				} else {
					// log.warn(acc + " not annotated");
				}
			}
		}
		final List<String> list = new ArrayList<String>();
		list.addAll(set);
		Collections.sort(list);

		final StringBuilder sb = new StringBuilder();

		for (final String string : list) {
			if (!"".equals(sb.toString()))
				sb.append(",");
			sb.append(string);
		}

		return sb.toString();
	}

	public static String getGeneNameString(Map<String, Entry> annotatedProteins,
			Map<String, PCQProteinNode> proteinNodeMap, Set<String> validTaxonomies, boolean onlyFirst) {

		final Set<String> set = new THashSet<String>();
		for (final String rawAcc : proteinNodeMap.keySet()) {

			final List<String> accs = new ArrayList<String>();
			if (rawAcc.contains(" ")) {
				final String[] split = rawAcc.split(" ");
				for (final String string : split) {
					accs.add(string);
				}
			} else {
				accs.add(rawAcc);
			}

			// int index = 0;
			for (final String acc : accs) {
				if (annotatedProteins != null && annotatedProteins.containsKey(acc)) {
					final String taxon = getTaxonomy(acc, annotatedProteins.get(acc).getOrganism());
					if (taxon != null) {
						boolean valid = false;
						if (validTaxonomies != null && !validTaxonomies.isEmpty()) {
							for (final String skipTaxonomy : validTaxonomies) {
								if (taxon.toLowerCase().contains(skipTaxonomy.toLowerCase())) {
									valid = true;
								}
							}
						} else {
							valid = true;
						}
						if (!valid) {
							continue;
						}
					}

					final String geneName = getGeneName(annotatedProteins.get(acc).getGene());
					set.add(geneName);
				} else {
					// log.warn(acc + " not annotated");
				}
			}
		}
		final List<String> list = new ArrayList<String>();
		list.addAll(set);
		Collections.sort(list);
		if (onlyFirst && !list.isEmpty()) {
			return list.get(0);
		}
		final StringBuilder sb = new StringBuilder();
		if (!ProteinClusterQuantParameters.getInstance().isPrintOnlyFirstGene()) {
			for (final String string : list) {
				if (!"".equals(sb.toString()))
					sb.append(",");
				sb.append(string);
			}
		} else {
			if (!list.isEmpty()) {
				sb.append(list.get(0));
			}
		}

		return sb.toString();
	}

	private static String getTaxonomy(String acc, OrganismType organism) {
		if (organism != null && acc != null) {
			if (organism.getName() != null) {
				for (final OrganismNameType organismType : organism.getName()) {
					if (organismType.getType().equals("scientific")) {
						return organismType.getValue();
					}
				}
				if (!organism.getName().isEmpty()) {
					return organism.getName().get(0).getValue();
				}
			}
		}
		return null;
	}

	private static String getGeneName(List<GeneType> gene) {
		for (final GeneType geneType : gene) {
			for (final GeneNameType geneName : geneType.getName()) {
				if (geneName.getType().equals("primary"))
					return geneName.getValue();
			}
		}
		if (gene.isEmpty() || gene.get(0).getName().isEmpty()) {
			return "";
		}
		return gene.get(0).getName().get(0).getValue();
	}

	public static int getIonCount(PCQPeptideNode peptideNode) {
		int ionCount = 0;
		if (peptideNode != null) {
			final Set<QuantifiedPSMInterface> quantifiedPSMs = peptideNode.getQuantifiedPSMs();
			for (final QuantifiedPSMInterface psm : quantifiedPSMs) {
				for (final QuantificationLabel label : QuantificationLabel.values()) {
					ionCount += QuantUtils.getIonCount(psm, label);
				}
			}
		}
		return ionCount;
	}

	public static Map<String, QuantifiedPeptideInterface> getPeptideMapFromClusters(Set<ProteinCluster> clusterSet) {
		final Map<String, QuantifiedPeptideInterface> map = new THashMap<String, QuantifiedPeptideInterface>();
		for (final ProteinCluster proteinCluster : clusterSet) {
			final Set<QuantifiedPeptideInterface> peptideSet = proteinCluster.getPeptideSet();
			for (final QuantifiedPeptideInterface quantifiedPeptide : peptideSet) {
				final String peptideKey = quantifiedPeptide.getKey();
				if (map.containsKey(peptideKey)) {
					throw new IllegalArgumentException(
							"peptide key cannot be repeated. There is an incosistence b etween how Census parser is reading PSMs and creating peptides and how the peptideKey is created here");
				}
				map.put(peptideKey, quantifiedPeptide);
			}
		}
		return map;
	}

	public static List<String> getSortedTaxonomies(Set<QuantifiedProteinInterface> proteinSet) {
		final List<String> ret = new ArrayList<String>();
		for (final QuantifiedProteinInterface protein : proteinSet) {
			final Set<String> taxonomies = protein.getTaxonomies();
			for (final String taxonomy : taxonomies) {
				if (taxonomy != null && !ret.contains(taxonomy)) {
					ret.add(taxonomy);
				}
			}
		}
		Collections.sort(ret);
		return ret;
	}

	public static boolean peptidesShareAllProteins(QuantifiedPeptideInterface peptide1,
			QuantifiedPeptideInterface peptide2) {
		final Set<QuantifiedProteinInterface> proteins1 = peptide1.getQuantifiedProteins();
		final Set<QuantifiedProteinInterface> proteins2 = peptide2.getQuantifiedProteins();

		final Set<String> proteinAccs1 = PCQUtils.getAccessions(proteins1);
		final Set<String> proteinAccs2 = PCQUtils.getAccessions(proteins2);
		if (proteinAccs1.size() == proteinAccs2.size()) {
			for (final String acc1 : proteinAccs1) {
				if (!proteinAccs2.contains(acc1)) {
					return false;
				}
			}
			return true;
		}
		return false;
	}

	private static Set<String> getAccessions(Set<QuantifiedProteinInterface> proteins) {
		final Set<String> accs = new THashSet<String>();
		for (final QuantifiedProteinInterface protein : proteins) {
			accs.add(protein.getAccession());
		}
		return accs;
	}

	public static PCQPeptideNode mergePeptideNodes(PCQPeptideNode peptideNode, PCQPeptideNode peptideNode2) {
		if (peptideNode.hashCode() == peptideNode2.hashCode()) {
			return peptideNode;
		}
		final Set<QuantifiedPeptideInterface> peptideCollection = new THashSet<QuantifiedPeptideInterface>();
		peptideCollection.addAll(peptideNode.getQuantifiedPeptides());
		peptideCollection.addAll(peptideNode2.getQuantifiedPeptides());
		final PCQPeptideNode ret = new PCQPeptideNode(peptideNode.getCluster(), peptideCollection);
		return ret;
	}

	public static void mergeProteinNodes(PCQProteinNode proteinNode, PCQProteinNode proteinNode2) {
		if (proteinNode.hashCode() == proteinNode2.hashCode()) {
			return;
		}
		final Set<QuantifiedProteinInterface> proteins2 = proteinNode2.getQuantifiedProteins();
		for (final QuantifiedProteinInterface quantifiedProteinInterface : proteins2) {
			proteinNode.addProtein(quantifiedProteinInterface);
		}
	}

	public static QuantParser getQuantParser(ProteinClusterQuantParameters params,
			List<Map<QuantCondition, QuantificationLabel>> labelsByConditionsList) throws FileNotFoundException {
		return getQuantParser(params, labelsByConditionsList, params.getQuantInputFileNamesArray());
	}

	public static QuantParser getQuantParser(ProteinClusterQuantParameters params,
			List<Map<QuantCondition, QuantificationLabel>> labelsByConditionsList, final String[] inputFileNamesArray)
			throws FileNotFoundException {
		log.debug("Getting input file parser");
		if (params.getAnalysisInputType() == AnalysisInputType.CENSUS_CHRO) {
			if (params.getMongoDBURI() != null) {
				return PCQUtils.getCensusChroParserUsingMongoDBIndex(params.getMongoDBURI(),
						params.getMongoMassDBName(), params.getMongoSeqDBName(), params.getMongoProtDBName(),
						params.getInputFileFolder(), inputFileNamesArray, labelsByConditionsList,
						params.getNumeratorLabel(), params.getDenominatorLabel(), params.getUniprotReleasesFolder(),
						params.getUniprotVersion(), params.getDecoyRegexp(), params.isIgnoreNotFoundPeptidesInDB(),
						!params.isIgnorePTMs(), params.getPeptideFilterRegexp(), params.getAaQuantified());
			} else {
				return PCQUtils.getCensusChroParser(params.getFastaFile(), params.getInputFileFolder(),
						inputFileNamesArray, labelsByConditionsList, params.getNumeratorLabel(),
						params.getDenominatorLabel(), params.getEnzymeArray(), params.getMissedCleavages(),
						params.isSemiCleavage(), params.getUniprotReleasesFolder(), params.getUniprotVersion(),
						params.getDecoyRegexp(), params.isIgnoreNotFoundPeptidesInDB(), !params.isIgnorePTMs(),
						params.getPeptideFilterRegexp(), params.getAaQuantified(), params.isLookForProteoforms());
			}
		} else if (params.getAnalysisInputType() == AnalysisInputType.CENSUS_OUT) {
			if (params.getMongoDBURI() != null) {
				final CensusOutParser parser = PCQUtils.getCensusOutParserUsingMongoDBIndex(params.getMongoDBURI(),
						params.getMongoMassDBName(), params.getMongoSeqDBName(), params.getMongoProtDBName(),
						params.getInputFileFolder(), inputFileNamesArray, labelsByConditionsList,
						params.getNumeratorLabel(), params.getDenominatorLabel(), params.getUniprotReleasesFolder(),
						params.getUniprotVersion(), params.getDecoyRegexp(), params.isIgnoreNotFoundPeptidesInDB(),
						params.isOnlyOneSpectrumPerChromatographicPeakAndPerSaltStep(), params.isSkipSingletons(),
						!params.isIgnorePTMs(), params.getPeptideFilterRegexp(), params.getAaQuantified());

				return parser;
			} else {
				final CensusOutParser parser = PCQUtils.getCensusOutParser(params.getFastaFile(),
						params.getInputFileFolder(), inputFileNamesArray, labelsByConditionsList,
						params.getNumeratorLabel(), params.getDenominatorLabel(), params.getEnzymeArray(),
						params.getMissedCleavages(), params.isSemiCleavage(), params.getUniprotReleasesFolder(),
						params.getUniprotVersion(), params.getDecoyRegexp(), params.isIgnoreNotFoundPeptidesInDB(),
						params.isOnlyOneSpectrumPerChromatographicPeakAndPerSaltStep(), params.isSkipSingletons(),
						!params.isIgnorePTMs(), params.getPeptideFilterRegexp(), params.getAaQuantified(),
						params.isLookForProteoforms());

				return parser;
			}
		} else if (params.getAnalysisInputType() == AnalysisInputType.SEPARATED_VALUES) {
			if (params.getMongoDBURI() != null) {
				final SeparatedValuesParser parser = PCQUtils.getSeparatedValuesParserUsingMongoDBIndex(
						params.getMongoDBURI(), params.getMongoMassDBName(), params.getMongoSeqDBName(),
						params.getMongoProtDBName(), params.getInputFileFolder(), inputFileNamesArray,
						params.getSeparator(), labelsByConditionsList, params.getNumeratorLabel(),
						params.getDenominatorLabel(), params.getUniprotReleasesFolder(), params.getUniprotVersion(),
						params.getDecoyRegexp(), params.isIgnoreNotFoundPeptidesInDB(), !params.isIgnorePTMs(),
						params.getPeptideFilterRegexp(), params.getAaQuantified());

				return parser;
			} else {
				final SeparatedValuesParser parser = PCQUtils.getSeparatedValuesParser(params.getFastaFile(),
						params.getInputFileFolder(), inputFileNamesArray, params.getSeparator(), labelsByConditionsList,
						params.getNumeratorLabel(), params.getDenominatorLabel(), params.getEnzymeArray(),
						params.getMissedCleavages(), params.isSemiCleavage(), params.getUniprotReleasesFolder(),
						params.getUniprotVersion(), params.getDecoyRegexp(), params.isIgnoreNotFoundPeptidesInDB(),
						!params.isIgnorePTMs(), params.getPeptideFilterRegexp(), params.getAaQuantified(),
						params.isLookForProteoforms());

				return parser;
			}
		}
		return null;

	}

	public static DTASelectParser getDTASelectParser(ProteinClusterQuantParameters params)
			throws FileNotFoundException {
		if (params.getIdentificationInputFileNamesArray() != null
				&& params.getIdentificationInputFileNamesArray().length > 0) {
			return getDTASelectParser(params, params.getIdentificationInputFileNamesArray());
		}
		return null;
	}

	public static DTASelectParser getDTASelectParser(ProteinClusterQuantParameters params,
			final String[] inputFileNamesArray) throws FileNotFoundException {
		log.debug("Getting input file parser");

		if (params.getMongoDBURI() != null) {
			return getDTASelectParserUsingMongoDBIndex(params.getMongoDBURI(), params.getMongoMassDBName(),
					params.getMongoSeqDBName(), params.getMongoProtDBName(), params.getInputFileFolder(),
					inputFileNamesArray, params.getUniprotReleasesFolder(), params.getUniprotVersion(),
					params.getDecoyRegexp(), params.isIgnoreNotFoundPeptidesInDB(), params.getPeptideFilterRegexp());
		} else {
			return getDTASelectParser(params.getFastaFile(), params.getInputFileFolder(), inputFileNamesArray,
					params.getEnzymeArray(), params.getMissedCleavages(), params.isSemiCleavage(),
					params.getUniprotReleasesFolder(), params.getUniprotVersion(), params.getDecoyRegexp(),
					params.isIgnoreNotFoundPeptidesInDB(), params.getPeptideFilterRegexp(),
					params.isLookForProteoforms());
		}
	}

	public static QuantParser getQuantParser(ProteinClusterQuantParameters params,
			Map<QuantCondition, QuantificationLabel> labelsByConditions, String inputFileName)
			throws FileNotFoundException {
		log.debug("Getting input file parser");
		final List<Map<QuantCondition, QuantificationLabel>> labelsByConditionsList = new ArrayList<Map<QuantCondition, QuantificationLabel>>();
		labelsByConditionsList.add(labelsByConditions);
		final String[] inputFileNamesArray = new String[1];
		inputFileNamesArray[0] = inputFileName;
		return getQuantParser(params, labelsByConditionsList, inputFileNamesArray);
	}

	/**
	 * Gets the ratio value from a {@link QuantRatio}. In case of having some
	 * null value, it will return Double.NaN.
	 *
	 * @param ratio
	 * @param cond1
	 * @param cond2
	 * @return
	 */
	public static Double getLog2RatioValue(QuantRatio ratio, QuantCondition cond1, QuantCondition cond2) {
		Double ratioValue = Double.NaN;
		if (ratio != null) {

			final Double log2Ratio = ratio.getLog2Ratio(cond1, cond2);
			if (log2Ratio != null) {
				ratioValue = log2Ratio;
			}

		}
		return ratioValue;
	}

	public static Map<String, Set<QuantifiedProteinInterface>> getIndividualProteinsMap(PCQPeptideNode peptideNode) {
		final Map<String, Set<QuantifiedProteinInterface>> ret = new THashMap<String, Set<QuantifiedProteinInterface>>();
		for (final QuantifiedPeptideInterface peptide : peptideNode.getQuantifiedPeptides()) {
			final Set<QuantifiedProteinInterface> quantifiedProteins = peptide.getQuantifiedProteins();
			for (final QuantifiedProteinInterface quantifiedProteinInterface : quantifiedProteins) {
				final String accession = quantifiedProteinInterface.getAccession();
				if (ret.containsKey(accession)) {
					ret.get(accession).add(quantifiedProteinInterface);
				} else {
					final Set<QuantifiedProteinInterface> set = new THashSet<QuantifiedProteinInterface>();
					set.add(quantifiedProteinInterface);
					ret.put(accession, set);
				}
			}
		}
		return ret;
	}

	/**
	 * Gets the peptide nodes that are unique to a protein node with respect to
	 * another protein node.
	 *
	 * @param proteinNode1
	 * @param proteinNode2
	 * @param uniquePepOnly
	 *            if true, only the peptide nodes that are really unique to
	 *            proteinNode1 are reported. If false, it will return also the
	 *            peptide nodes that are unique to proteinNode1 with respect to
	 *            proteinNode2 but may be shared with a third protein
	 * @param skipDiscarded
	 * @return
	 */
	public static Set<PCQPeptideNode> getUniquePeptideNodes(PCQProteinNode proteinNode1, PCQProteinNode proteinNode2,
			boolean uniquePepOnly, boolean skipDiscarded) {
		final Set<PCQPeptideNode> peptideNodes1 = proteinNode1.getPeptideNodes();
		final Set<PCQPeptideNode> ret = new THashSet<PCQPeptideNode>();

		for (final PCQPeptideNode peptideNode1 : peptideNodes1) {
			if (skipDiscarded && peptideNode1.isDiscarded()) {
				continue;
			}
			// if we require the unique peptide nodes to be really unique (not
			// shared by a third protein)
			if (uniquePepOnly) {
				if (peptideNode1.getProteinNodes().size() == 1) {
					ret.add(peptideNode1);
				} else {
					continue;
				}
				// if we dont care wether the peptide node is shared by a third
				// protein
			} else {
				// if it is shared, not take it
				if (peptideNode1.getProteinNodes().contains(proteinNode2)) {
					continue;
				} else {
					ret.add(peptideNode1);
				}
			}
		}
		return ret;
	}

	public static boolean proteinsShareAllPeptides(Collection<QuantifiedProteinInterface> proteins1,
			Collection<QuantifiedProteinInterface> proteins2) {
		final Set<QuantifiedPeptideInterface> peptides1 = new THashSet<QuantifiedPeptideInterface>();
		for (final QuantifiedProteinInterface protein1 : proteins1) {
			peptides1.addAll(protein1.getQuantifiedPeptides());
		}
		final Set<QuantifiedPeptideInterface> peptides2 = new THashSet<QuantifiedPeptideInterface>();
		for (final QuantifiedProteinInterface protein2 : proteins2) {
			peptides2.addAll(protein2.getQuantifiedPeptides());
		}
		if (peptides1.size() == peptides2.size()) {
			for (final QuantifiedPeptideInterface peptide2 : peptides2) {
				if (!peptides1.contains(peptide2)) {
					return false;
				}
			}
			return true;
		}
		return false;
	}

	public static Set<QuantRatio> getNonInfinityIndividualPSMRatios(PCQPeptideNode peptideNode) {
		if (peptideNode == null) {
			return null;
		}
		final Set<QuantRatio> ret = new THashSet<QuantRatio>();
		final Set<QuantifiedPSMInterface> quantifiedPSMs = peptideNode.getQuantifiedPSMs();
		for (final QuantifiedPSMInterface quantifiedPSMInterface : quantifiedPSMs) {
			ret.addAll(quantifiedPSMInterface.getNonInfinityRatios());
		}
		return ret;
	}

	public static Set<QuantRatio> getNonInfinityIndividualPeptideRatios(PCQPeptideNode peptideNode) {
		if (peptideNode == null) {
			return null;
		}
		final Set<QuantRatio> ret = new THashSet<QuantRatio>();
		final Set<QuantifiedPeptideInterface> quantifiedPeptides = peptideNode.getQuantifiedPeptides();
		for (final QuantifiedPeptideInterface quantifiedPeptide : quantifiedPeptides) {
			ret.addAll(quantifiedPeptide.getNonInfinityRatios());
		}
		return ret;
	}

	public static Double averageOfRatiosTakingIntoAccountInfinitiesAndNans(Collection<QuantRatio> ratios,
			QuantCondition cond1, QuantCondition cond2) {
		final List<Double> ratioValues = new ArrayList<Double>();
		for (final QuantRatio ratio : ratios) {
			if (ratio != null) {
				ratioValues.add(ratio.getLog2Ratio(cond1, cond2));
			}
		}
		if (ratioValues.isEmpty()) {
			return null;
		}
		// check if there are all INFINITIES
		final boolean areInfinities = areAll(Double.POSITIVE_INFINITY, ratioValues)
				|| areAll(Double.NEGATIVE_INFINITY, ratioValues) || areAll(Double.MAX_VALUE, ratioValues)
				|| areAll(-Double.MAX_VALUE, ratioValues);
		if (areInfinities) {
			// return it (we assume is only one sign of the infinities here
			return ratioValues.iterator().next();
		} else {
			// if they are all Nan,return nan
			if (areAll(Double.NaN, ratioValues)) {
				return Double.NaN;
			}
			// return an average of the non infinities
			final List<Double> nonInfinityNonNanValues = new ArrayList<Double>();
			for (final Double ratioValue : ratioValues) {
				if (!ratioValue.isInfinite() && !ratioValue.isNaN()) {
					nonInfinityNonNanValues.add(ratioValue);
				}
			}
			// report the average
			return Maths.mean(nonInfinityNonNanValues.toArray(new Double[0]));
		}
	}

	public static Double stdevOfRatiosTakingIntoAccountInfinitiesAndNans(Collection<QuantRatio> ratios,
			QuantCondition cond1, QuantCondition cond2) {
		final List<Double> ratioValues = new ArrayList<Double>();
		for (final QuantRatio ratio : ratios) {
			if (ratio != null) {
				ratioValues.add(ratio.getLog2Ratio(cond1, cond2));
			}
		}
		if (ratioValues.isEmpty()) {
			return null;
		}
		// check if there are all INFINITIES
		final boolean areInfinities = areAll(Double.POSITIVE_INFINITY, ratioValues)
				|| areAll(Double.NEGATIVE_INFINITY, ratioValues) || areAll(Double.MAX_VALUE, ratioValues)
				|| areAll(-Double.MAX_VALUE, ratioValues);
		if (areInfinities) {
			// return it (we assume is only one sign of the infinities here
			return ratioValues.iterator().next();
		} else {
			// if they are all Nan,return nan
			if (areAll(Double.NaN, ratioValues)) {
				return Double.NaN;
			}
			// return an average of the non infinities
			final List<Double> nonInfinityNonNanValues = new ArrayList<Double>();
			for (final Double ratioValue : ratioValues) {
				if (!ratioValue.isInfinite() && !ratioValue.isNaN()) {
					nonInfinityNonNanValues.add(ratioValue);
				}
			}
			// report the stdev
			return Maths.stddev(nonInfinityNonNanValues.toArray(new Double[0]));
		}
	}

	private static boolean areAll(double value, Collection<Double> values) {
		if (values.isEmpty()) {
			return false;
		}
		for (final Double double1 : values) {
			if (!double1.equals(value)) {
				return false;
			}
		}
		return true;
	}

	public static QuantRatio getRepresentativeRatioForPeptideNode(PCQPeptideNode peptideNode, QuantCondition cond1,
			QuantCondition cond2, String replicateName, boolean skipDiscarded) {
		final Set<PCQPeptideNode> set = new THashSet<PCQPeptideNode>();
		set.add(peptideNode);
		return getRepresentativeRatioForPeptideNodes(set, cond1, cond2, replicateName, skipDiscarded);
	}

	/**
	 * Gets the representative log2 ratio value for the classifications of the
	 * protein pairs, which is:<br>
	 * - in case of being isobaric isotopologues (analysisInputType=
	 * {@link AnalysisInputType}=CENSUS_CHRO) and SanXot is enabled, the average
	 * of the peptideNode Ri ratios coming from SanXot<br>
	 * - in case of being isobaric isotopologues and SanXot is not enabled,
	 * depending on the parameter isobaricRatioType, it will be the average of
	 * the peptideNode average Rc or Ri ratios of the individual peptides in the
	 * node <br>
	 * - in case of other quantification techniques, if SanXot is enabled, it is
	 * the average of the SanXot ratios of the peptide nodes coming from SanXot.
	 * <br>
	 * - in case of other quantification techniques and SanXot is not enabled,
	 * it is the average of the PSM ratios.<br>
	 *
	 *
	 * @param peptideNodes
	 * @param cond1
	 * @param cond2
	 * @param skipDiscarded
	 * @return
	 */
	public static QuantRatio getRepresentativeRatioForPeptideNodes(Collection<PCQPeptideNode> peptideNodes,
			QuantCondition cond1, QuantCondition cond2, String replicateName, boolean skipDiscarded) {
		if (peptideNodes.isEmpty()) {
			return CensusRatio.getNaNRatio(cond1, cond2, AggregationLevel.PEPTIDE_NODE, "RATIO");
		}
		final List<Integer> quantifiedSitePositionInPeptideList = new ArrayList<Integer>();
		final List<QuantRatio> toAverage = new ArrayList<QuantRatio>();
		String avgRatioDescription = "";
		// SANXOT
		// it doesn't matter from what is coming from.
		// we use an average of the sanxot ratios
		final ProteinClusterQuantParameters params = ProteinClusterQuantParameters.getInstance();
		if (params.isPerformRatioIntegration()) {
			avgRatioDescription = peptideNodes.size() < 2 ? "Integrated ratios" : "Avg of integrated ratios";
			// SANXOT ENABLED or use of Ri ratios
			// get the average of the sanxot ratios
			for (final PCQPeptideNode peptideNode : peptideNodes) {
				if (skipDiscarded && peptideNode.isDiscarded()) {
					continue;
				}
				final QuantRatio consensusRatio = peptideNode.getSanXotRatio(cond1, cond2, replicateName);
				if (consensusRatio != null) {
					toAverage.add(consensusRatio);
					if (consensusRatio.getQuantifiedSitePositionInPeptide() != null) {
						quantifiedSitePositionInPeptideList.add(consensusRatio.getQuantifiedSitePositionInPeptide());
					}
				}
			}
		} else {
			// only for collapseBySites=TRUE
			final Set<Pair<IsobaricQuantifiedPeptide, PositionInPeptide>> isobaricpeptidesAndPositionsInPeptides = new THashSet<Pair<IsobaricQuantifiedPeptide, PositionInPeptide>>();
			final Set<Pair<QuantifiedPeptideInterface, PositionInPeptide>> peptidesAndPositionsInPeptides = new THashSet<Pair<QuantifiedPeptideInterface, PositionInPeptide>>();
			//
			if (params.isCollapseBySites()) {
				for (final PCQPeptideNode peptideNode : peptideNodes) {
					if (skipDiscarded && peptideNode.isDiscarded()) {
						continue;
					}
					for (final QuantifiedPeptideInterface peptide : peptideNode.getQuantifiedPeptides()) {
						if (peptide instanceof IsobaricQuantifiedPeptide) {
							final Pair<IsobaricQuantifiedPeptide, PositionInPeptide> pair = new Pair<IsobaricQuantifiedPeptide, PositionInPeptide>(
									(IsobaricQuantifiedPeptide) peptide, peptideNode.getPositionInPeptide(peptide));
							isobaricpeptidesAndPositionsInPeptides.add(pair);
						} else {
							final Pair<QuantifiedPeptideInterface, PositionInPeptide> pair = new Pair<QuantifiedPeptideInterface, PositionInPeptide>(
									peptide, peptideNode.getPositionInPeptide(peptide));
							peptidesAndPositionsInPeptides.add(pair);
						}
					}
				}
			}
			if (params.getAnalysisInputType() == AnalysisInputType.CENSUS_CHRO) {

				if (params.getIsobaricRatioType() == IsobaricRatioType.Ri) {

					avgRatioDescription = peptideNodes.size() < 2 ? "Ri ratio" : "Avg of Ri ratios";
					for (final PCQPeptideNode peptideNode : peptideNodes) {
						if (skipDiscarded && peptideNode.isDiscarded()) {
							continue;
						}
						if (params.isCollapseBySites()) {
							final QuantRatio isobaricRatiosForSiteSpecificPeptideNode = getAverageIsobaricRatioForSiteSpecificPeptideNode(
									peptideNode, cond1, cond2);
							if (isobaricRatiosForSiteSpecificPeptideNode != null) {
								toAverage.add(isobaricRatiosForSiteSpecificPeptideNode);
							}
						} else {
							final QuantRatio averageRatio = QuantUtils.getAverageRatio(QuantUtils.getNonInfinityRatios(
									QuantUtils.getRatiosByName(peptideNode.getRatios(), getRatioNameByAnalysisType())),
									AggregationLevel.PEPTIDE_NODE);
							if (averageRatio != null) {
								toAverage.add(averageRatio);
							}
						}

					}
				} else if (params.getIsobaricRatioType() == IsobaricRatioType.Rc) {
					// use of Rc ratios get the average of the normalized Rc
					// ratios for each peptide node
					avgRatioDescription = peptideNodes.size() < 2 ? "Normalized Rc ratio"
							: "Avg of normalized Rc ratios";
					for (final PCQPeptideNode peptideNode : peptideNodes) {
						if (skipDiscarded && peptideNode.isDiscarded()) {
							continue;
						}
						if (params.isCollapseBySites()) {
							final IonCountRatio normalizedIonCountRatioForPeptidesForQuantifiedSites = QuantUtils
									.getNormalizedIonCountRatioForPeptidesForQuantifiedSites(
											isobaricpeptidesAndPositionsInPeptides, cond1, cond2, replicateName,
											params.getAaQuantified());
							if (normalizedIonCountRatioForPeptidesForQuantifiedSites != null) {
								toAverage.add(normalizedIonCountRatioForPeptidesForQuantifiedSites);
							}
						} else {
							final IonCountRatio normalizedIonCountRatioForPeptideNode = getNormalizedIonCountRatioForPeptideNode(
									peptideNode, cond1, cond2, replicateName);
							if (normalizedIonCountRatioForPeptideNode != null) {
								toAverage.add(normalizedIonCountRatioForPeptideNode);
							}
						}
					}
				} else {
					throw new IllegalArgumentException("Option not possible");
				}

			} else {
				// no census chro
				// SILAC or others,

				// get the average of the individual psm ratios
				avgRatioDescription = "Avg of individual PSM ratios";
				for (final PCQPeptideNode peptideNode : peptideNodes) {
					if (skipDiscarded && peptideNode.isDiscarded()) {
						continue;
					}
					if (params.isCollapseBySites()) {
						final QuantRatio ratioForSiteSpecificPeptideNode = getAverageRatioForSiteSpecificPeptideNode(
								peptideNode, cond1, cond2);
						if (ratioForSiteSpecificPeptideNode != null) {
							toAverage.add(ratioForSiteSpecificPeptideNode);
						}
					} else {
						final Set<QuantifiedPSMInterface> quantifiedPSMs = peptideNode.getQuantifiedPSMs();
						for (final QuantifiedPSMInterface psm : quantifiedPSMs) {
							if (replicateName != null && !psm.getFileNames().contains(replicateName)) {
								continue;
							}
							if (psm instanceof QuantifiedPSM) {
								final QuantRatio validRatio = QuantUtils.getRatioByName(psm,
										getRatioNameByAnalysisType());
								if (validRatio != null) {
									toAverage.add(validRatio);
									if (validRatio.getQuantifiedSitePositionInPeptide() != null) {
										quantifiedSitePositionInPeptideList
												.add(validRatio.getQuantifiedSitePositionInPeptide());
									}
								}
							} else {
								throw new IllegalArgumentException(
										"In case of SILAC,  it has to be a QuantifiedPSMFromCensusOut");
							}
						}
					}
				}
			}

		}
		// if there is only one, return it in order to not loose the extended
		// class and description of the ratio
		if (toAverage.size() == 1) {
			return toAverage.get(0);
		}
		final Double finalValue = PCQUtils.averageOfRatiosTakingIntoAccountInfinitiesAndNans(toAverage, cond1, cond2);
		if (finalValue != null) {
			final CensusRatio censusRatio = new CensusRatio(finalValue, true, cond1, cond2,
					AggregationLevel.PEPTIDE_NODE, avgRatioDescription);
			if (quantifiedSitePositionInPeptideList.size() == 1) {
				censusRatio.setQuantifiedSitePositionInPeptide(quantifiedSitePositionInPeptideList.get(0));
			}
			final Double stdev = PCQUtils.stdevOfRatiosTakingIntoAccountInfinitiesAndNans(toAverage, cond1, cond2);
			if (stdev != null) {
				censusRatio.setRatioScore(new RatioScore(String.valueOf(stdev), "STDEV", "Standard deviation of ratios",
						"Standard deviation of the ratios averaged"));
			}

			return censusRatio;
		}
		return CensusRatio.getNaNRatio(cond1, cond2, AggregationLevel.PEPTIDE_NODE, "RATIO");

	}

	/**
	 * The peptide node MUST have been created as a site specific peptide node.
	 * <br>
	 * It retrieves an average of all the non-infinity isobaric ratios per
	 * peptide that are valid for the site that is quantified.
	 * 
	 * @param peptideNode
	 * @param cond1
	 * @param cond2
	 * @return
	 */
	private static QuantRatio getAverageRatioForSiteSpecificPeptideNode(PCQPeptideNode peptideNode,
			QuantCondition cond1, QuantCondition cond2) {
		final List<QuantRatio> toAverage = new ArrayList<QuantRatio>();
		final List<Pair<QuantifiedPeptideInterface, PositionInPeptide>> peptidesWithPositionsInPeptide = peptideNode
				.getPeptidesWithPositionsInPeptide();
		for (final Pair<QuantifiedPeptideInterface, PositionInPeptide> pair : peptidesWithPositionsInPeptide) {
			final QuantifiedPeptideInterface peptide = pair.getFirstelement();
			final int positionInPeptide = pair.getSecondElement().getPosition();
			final List<QuantRatio> ratios = QuantUtils.getRatiosByName(peptide, getRatioNameByAnalysisType());
			for (final QuantRatio quantRatio : ratios) {
				final Integer quantifiedSitePositionInPeptide = quantRatio.getQuantifiedSitePositionInPeptide();
				if (quantifiedSitePositionInPeptide != null
						&& quantifiedSitePositionInPeptide.equals(positionInPeptide)) {
					toAverage.add(quantRatio);
				}
			}
		}
		final Double finalValue = PCQUtils.averageOfRatiosTakingIntoAccountInfinitiesAndNans(toAverage, cond1, cond2);
		if (finalValue != null) {
			final CensusRatio censusRatio = new CensusRatio(finalValue, true, cond1, cond2,
					AggregationLevel.PEPTIDE_NODE, "Avg ratios for site");
			final Double stdev = PCQUtils.stdevOfRatiosTakingIntoAccountInfinitiesAndNans(toAverage, cond1, cond2);
			if (stdev != null) {
				censusRatio.setRatioScore(new RatioScore(String.valueOf(stdev), "STDEV", "Standard deviation of ratios",
						"Standard deviation of the ratios averaged"));
			}
			return censusRatio;
		}
		return CensusRatio.getNaNRatio(cond1, cond2, AggregationLevel.PEPTIDE_NODE, "RATIO");
	}

	/**
	 * The peptide node MUST have been created as a site specific peptide node.
	 * <br>
	 * It retrieves an average of all the non-infinity isobaric ratios per
	 * peptide that are valid for the site that is quantified.
	 * 
	 * @param peptideNode
	 * @param cond1
	 * @param cond2
	 * @return
	 */
	private static QuantRatio getAverageIsobaricRatioForSiteSpecificPeptideNode(PCQPeptideNode peptideNode,
			QuantCondition cond1, QuantCondition cond2) {
		final List<QuantRatio> toAverage = new ArrayList<QuantRatio>();

		final List<Pair<QuantifiedPeptideInterface, PositionInPeptide>> peptidesWithPositionsInPeptide = peptideNode
				.getPeptidesWithPositionsInPeptide();
		for (final Pair<QuantifiedPeptideInterface, PositionInPeptide> pair : peptidesWithPositionsInPeptide) {
			final QuantifiedPeptideInterface peptide = pair.getFirstelement();
			if (peptide instanceof IsobaricQuantifiedPeptide) {
				final int positionInPeptide = pair.getSecondElement().getPosition();
				final Set<IsoRatio> isoRatios = QuantUtils
						.getIsobaricRatiosForSiteFromPeptide((IsobaricQuantifiedPeptide) peptide, positionInPeptide);
				toAverage.addAll(isoRatios);
			}
		}
		final Double finalValue = PCQUtils.averageOfRatiosTakingIntoAccountInfinitiesAndNans(toAverage, cond1, cond2);
		if (finalValue != null) {
			final CensusRatio censusRatio = new CensusRatio(finalValue, true, cond1, cond2,
					AggregationLevel.PEPTIDE_NODE, "Avg Ri ratios for site");
			final Double stdev = PCQUtils.stdevOfRatiosTakingIntoAccountInfinitiesAndNans(toAverage, cond1, cond2);
			if (stdev != null) {
				censusRatio.setRatioScore(new RatioScore(String.valueOf(stdev), "STDEV", "Standard deviation of ratios",
						"Standard deviation of the ratios averaged"));
			}
			return censusRatio;
		}
		return CensusRatio.getNaNRatio(cond1, cond2, AggregationLevel.PEPTIDE_NODE, "RATIO");
	}

	/**
	 * Gets the individual log2 ratio values for the classification1
	 * (statistical outlier test) of the protein pairs, which is:<br>
	 * - in case of being isobaric isotopologues (analysisInputType=
	 * {@link AnalysisInputType}=CENSUS_CHRO) and SanXot is enabled, the
	 * individual Ri ratios of the PSMs<br>
	 * - in case of being isobaric isotopolofues and SanXot is not enabled, the
	 * individual Rc ratios of the individual peptides in the node<br>
	 * - in case of other quantification techniques, if SanXot is enabled, it is
	 * the individual psm ratios.<br>
	 * - in case of other quantification techniques and SanXot is not enabled,
	 * it is also the individual psm ratios.<br>
	 *
	 *
	 * @param peptideNodes
	 * @param cond1
	 * @param cond2
	 * @param skipDiscarded
	 * @return
	 **/
	public static List<Double> getIndividualRepresentativeLog2ValuesForEachPeptideForProteinPairAnalysis(
			Collection<PCQPeptideNode> peptideNodes, QuantCondition cond1, QuantCondition cond2,
			boolean skipDiscarded) {

		final List<Double> toAverage = new ArrayList<Double>();
		// ISOTOPOLOGUES
		if (ProteinClusterQuantParameters.getInstance().getAnalysisInputType() == AnalysisInputType.CENSUS_CHRO) {

			if (ProteinClusterQuantParameters.getInstance().isPerformRatioIntegration()) {
				// SANXOT ENABLED
				// get the individual Ri of the PSMs
				for (final PCQPeptideNode peptideNode : peptideNodes) {
					final Set<QuantifiedPSMInterface> quantifiedPSMs = peptideNode.getQuantifiedPSMs();
					for (final QuantifiedPSMInterface psm : quantifiedPSMs) {
						if (psm instanceof IsobaricQuantifiedPSM) {
							final Set<IsoRatio> nonInfinityIsoRatios = ((IsobaricQuantifiedPSM) psm)
									.getNonInfinityIsoRatios();
							for (final IsoRatio isoRatio : nonInfinityIsoRatios) {
								toAverage.add(isoRatio.getLog2Ratio(cond1, cond2));
							}
						}
					}
				}
			} else {
				// SANXOT DISABLED
				// get the individual Rc ratios of the individual peptides of
				// the node
				for (final PCQPeptideNode peptideNode : peptideNodes) {
					final Set<QuantifiedPeptideInterface> peptides = peptideNode.getQuantifiedPeptides();
					for (final QuantifiedPeptideInterface peptide : peptides) {
						if (peptide instanceof IsobaricQuantifiedPeptide) {
							toAverage.add(QuantUtils
									.getIonCountRatioForPeptide((IsobaricQuantifiedPeptide) peptide, cond1, cond2)
									.getLog2Ratio(cond1, cond2));
						}
					}
				}
			}

		} else {
			// SILAC or others,
			// get the individual psm ratios
			for (final PCQPeptideNode peptideNode : peptideNodes) {
				final Set<QuantifiedPSMInterface> quantifiedPSMs = peptideNode.getQuantifiedPSMs();
				for (final QuantifiedPSMInterface psm : quantifiedPSMs) {
					if (psm instanceof QuantifiedPSM) {
						final QuantRatio validRatio = QuantUtils.getRepresentativeRatio(psm);
						if (validRatio != null) {
							toAverage.add(validRatio.getLog2Ratio(cond1, cond2));
						}
					} else {
						throw new IllegalArgumentException(
								"In case of SILAC, it has to be a QuantifiedPSMFromCensusOut");
					}
				}
			}

		}
		return toAverage;
	}

	/**
	 * 
	 * @param peptideNode
	 * @param cond1
	 * @param cond2
	 * @param replicateName
	 * @return null if not isobaric peptides are in the node. a
	 *         {@link IonCountRatio} with Nan value if no ion counts are
	 *         present. an {@link IonCountRatio} otherwise.
	 */
	public static IonCountRatio getNormalizedIonCountRatioForPeptideNode(PCQPeptideNode peptideNode,
			QuantCondition cond1, QuantCondition cond2, String replicateName) {
		final List<IsobaricQuantifiedPeptide> isobaricPeptides = new ArrayList<IsobaricQuantifiedPeptide>();
		// return null if the peptide node doesn't contains any isobaric PSM
		boolean found = false;
		for (final QuantifiedPeptideInterface peptide : peptideNode.getQuantifiedPeptides()) {
			if (peptide instanceof IsobaricQuantifiedPeptide) {
				found = true;
			}
		}
		if (!found) {
			// the peptide node is not formed with isobaric peptides
			return null;
		}

		// when collapseBySites is TRUE we also need:
		final List<Pair<QuantifiedPeptideInterface, PositionInPeptide>> peptidesAndPositionsInPeptides = peptideNode
				.getPeptidesWithPositionsInPeptide();
		final List<Pair<IsobaricQuantifiedPeptide, PositionInPeptide>> isobaricPeptidesAndPositionsInPeptides = new ArrayList<Pair<IsobaricQuantifiedPeptide, PositionInPeptide>>();
		for (final Pair<QuantifiedPeptideInterface, PositionInPeptide> pair : peptidesAndPositionsInPeptides) {
			if (pair.getFirstelement() instanceof IsobaricQuantifiedPeptide) {
				isobaricPeptides.add((IsobaricQuantifiedPeptide) pair.getFirstelement());
				isobaricPeptidesAndPositionsInPeptides.add(new Pair<IsobaricQuantifiedPeptide, PositionInPeptide>(
						(IsobaricQuantifiedPeptide) pair.getFirstelement(), pair.getSecondElement()));
			}
		}
		if (peptidesAndPositionsInPeptides.isEmpty()) {
			return IonCountRatio.NAN_RATIO;
		}
		if (ProteinClusterQuantParameters.getInstance().isCollapseBySites()) {

			return QuantUtils.getNormalizedIonCountRatioForPeptidesForQuantifiedSites(
					isobaricPeptidesAndPositionsInPeptides, cond1, cond2, replicateName,
					ProteinClusterQuantParameters.getInstance().getAaQuantified());
		} else {
			return QuantUtils.getNormalizedIonCountRatioForPeptides(isobaricPeptides, cond1, cond2, replicateName);
		}
	}

	public static String escapeInfinity(Double ratio) {
		if (ratio == null) {
			return null;
		}
		if (Double.isInfinite(ratio)) {
			return "'" + String.valueOf(ratio);
		} else if (Double.compare(Double.MAX_VALUE, ratio) == 0) {
			return "'+INF";
		} else if (Double.compare(-Double.MAX_VALUE, ratio) == 0) {
			return "'-INF";
		} else {
			return String.valueOf(ratio);
		}
	}

	/**
	 * Returns the proteinPTMKey of a protein having a peptide with a PTM.<br>
	 * If a protein P12345 has the peptide ABCDE(+80)FGHI starting at position
	 * 20, the protein key would be P12345_+80(24).
	 *
	 * @param protein
	 * @param peptide
	 * @return
	 */
	public static String getProteinPTMKey(String accession, QuantifiedPeptideInterface... peptides) {
		final DecimalFormat df = new DecimalFormat("+#.##");
		final List<StringPosition> ptmPositionsInProtein = new ArrayList<StringPosition>();
		for (final QuantifiedPeptideInterface peptide : peptides) {
			ptmPositionsInProtein.addAll(QuantUtils.getPTMPositionsInProtein(accession, peptide,
					ProteinClusterQuantParameters.getInstance().getUniprotVersion(),
					ProteinClusterQuantParameters.getInstance().getUniprotReleasesFolder()));
		}
		final StringBuilder sb = new StringBuilder();
		sb.append(accession);
		if (!ptmPositionsInProtein.isEmpty()) {
			sb.append("_");
		}
		// sort by position
		Collections.sort(ptmPositionsInProtein, new Comparator<StringPosition>() {
			@Override
			public int compare(StringPosition o1, StringPosition o2) {
				return Integer.compare(o1.position, o2.position);
			}
		});
		// to not repeat:
		final TIntHashSet positions = new TIntHashSet();
		boolean first = true;
		for (final StringPosition ptmPositionInProtein : ptmPositionsInProtein) {
			if (positions.contains(ptmPositionInProtein.position)) {
				continue;
			}
			if (!first) {
				sb.append(",");
			}
			positions.add(ptmPositionInProtein.position);
			String string = ptmPositionInProtein.string;
			try {

				string = df.format(Double.valueOf(string));

			} catch (final NumberFormatException e) {

			}
			sb.append(ptmPositionInProtein.position).append("(").append(string).append(")");
			first = false;
		}
		return sb.toString();
	}

	public static Set<QuantifiedProteinInterface> getProteinsFromPeptides(
			Collection<QuantifiedPeptideInterface> peptides) {
		final Set<QuantifiedProteinInterface> ret = new THashSet<QuantifiedProteinInterface>();

		for (final QuantifiedPeptideInterface peptide : peptides) {
			ret.addAll(peptide.getQuantifiedProteins());
		}

		return ret;
	}

	public static String getSpeciesString(Set<String> taxonomies) {
		final StringBuilder sb = new StringBuilder();
		final List<String> list = new ArrayList<String>();
		list.addAll(taxonomies);
		Collections.sort(list);
		for (final String taxonomy : list) {
			if (!"".equals(sb.toString())) {
				sb.append(",");
			}
			sb.append(taxonomy);
		}
		return sb.toString();

	}

	public static String getPeptideNodeKeyByProteinSite(char[] aaQuantified, PCQPeptideNode pcqPeptideNode) {
		// TODO Auto-generated method stub
		return null;
	}

	public static String getTaxonomyString(PCQProteinNode proteinNode) {
		final StringBuilder sb = new StringBuilder();

		final Set<String> taxonomies = proteinNode.getTaxonomies();
		final List<String> list = new ArrayList<String>();
		for (final String tax : taxonomies) {
			if (tax != null) {
				list.add(tax);
			}
		}
		Collections.sort(list);
		for (final String tax : list) {
			if (!"".equals(sb.toString())) {
				sb.append(",");
			}
			sb.append(tax);
		}

		return sb.toString();
	}

	public static boolean areEquals(Collection<PositionInProtein> col1, Collection<PositionInProtein> col2) {
		if (col1.size() == col2.size()) {
			for (final PositionInProtein positionInProtein1 : col1) {
				if (!col2.contains(positionInProtein1)) {
					return false;
				}
			}
			for (final PositionInProtein positionInProtein2 : col2) {
				if (!col1.contains(positionInProtein2)) {
					return false;
				}
			}
			return true;
		}
		return false;
	}

	public static List<PositionInProtein> getUniqueToFirst(Collection<PositionInProtein> col1,
			Collection<PositionInProtein> col2) {
		final List<PositionInProtein> ret = new ArrayList<PositionInProtein>();

		for (final PositionInProtein positionInProtein1 : col1) {
			if (!col2.contains(positionInProtein1)) {
				ret.add(positionInProtein1);
			}
		}

		return ret;
	}

	/**
	 * Get a string such as: P12345#12#P23456#456 from, in this example, two
	 * proteins with position 12 and 456 respectively
	 * 
	 * @param keys1
	 * @return
	 */
	public static String getPositionsInProteinsKey(List<PositionInProtein> keys1) {
		Collections.sort(keys1, new Comparator<PositionInProtein>() {

			@Override
			public int compare(PositionInProtein o1, PositionInProtein o2) {
				final int compareTo = o1.getKey().compareTo(o2.getKey());
				if (compareTo == 0) {
					return Integer.compare(o1.getPosition(), o2.getPosition());
				}
				return compareTo;
			}
		});
		final StringBuilder sb = new StringBuilder();
		for (final PositionInProtein positionInProtein : keys1) {
			if (!"".equals(sb.toString())) {

				sb.append(KEY_SEPARATOR);
			}
			sb.append(positionInProtein);
		}
		return sb.toString();
	}

	public static boolean shareAny(List<PositionInProtein> col1, List<PositionInProtein> col2) {
		for (final PositionInProtein positionInProtein1 : col1) {
			if (col2.contains(positionInProtein1)) {
				return true;
			}
		}
		for (final PositionInProtein positionInProtein : col2) {
			if (col1.contains(positionInProtein)) {
				return true;
			}
		}
		return false;
	}

	public static List<PositionInProtein> getInCommon(List<PositionInProtein> col1, List<PositionInProtein> col2) {
		final List<PositionInProtein> ret = new ArrayList<PositionInProtein>();

		for (final PositionInProtein positionInProtein1 : col1) {
			if (col2.contains(positionInProtein1)) {
				ret.add(positionInProtein1);
			}
		}

		return ret;
	}

	public static boolean containsAny(String sequence, char[] quantifiedAAs) {
		for (final char c : quantifiedAAs) {
			if (sequence.contains(String.valueOf(c))) {
				return true;
			}
		}
		return false;
	}

	public static int howManyContains(String sequence, char[] quantifiedAAs) {
		int ret = 0;
		for (final char c : quantifiedAAs) {
			ret += StringUtils.allPositionsOf(sequence, c).size();
		}
		return ret;
	}

	/**
	 * this method returns the name of the ratio that is used depending on the
	 * analysis type.<br>
	 * This function would be the one to be called to get the ratios from PSMs
	 * and peptides, by calling to
	 * QuantUtils.getRatiosByName(getRatioNameByAnalysisType())
	 * 
	 * @return
	 */
	public static String getRatioNameByAnalysisType() {
		final ProteinClusterQuantParameters params = ProteinClusterQuantParameters.getInstance();
		final AnalysisInputType analysisInputType = params.getAnalysisInputType();
		if (analysisInputType == null) {
			return null;
		}
		switch (analysisInputType) {
		case CENSUS_CHRO:
			if (params.getIsobaricRatioType() == IsobaricRatioType.Rc) {
				return CensusChroParser.ISOBARIC_COUNT_RATIO;
			} else if (params.getIsobaricRatioType() == IsobaricRatioType.Ri) {
				return CensusChroParser.ISOBARIC_INTENSITY_RATIO;
			}
			break;
		case CENSUS_OUT:
			return CensusOutParser.AREA_RATIO;

		case SEPARATED_VALUES:
			return SeparatedValuesParser.RATIO;

		default:
			break;
		}
		throw new IllegalArgumentException("No ratio name defined for " + analysisInputType + " analysis");
	}
}
