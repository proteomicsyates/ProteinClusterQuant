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
import java.util.HashMap;
import java.util.HashSet;
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
import edu.scripps.yates.census.read.model.Ion;
import edu.scripps.yates.census.read.model.IonSerie.IonSerieType;
import edu.scripps.yates.census.read.model.IsobaricQuantifiedPSM;
import edu.scripps.yates.census.read.model.IsobaricQuantifiedPeptide;
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
import edu.scripps.yates.pcq.ProteinClusterQuantParameters;
import edu.scripps.yates.pcq.model.IonCountRatio;
import edu.scripps.yates.pcq.model.PCQPeptideNode;
import edu.scripps.yates.pcq.model.PCQProteinNode;
import edu.scripps.yates.pcq.model.ProteinCluster;
import edu.scripps.yates.utilities.alignment.nwalign.NWAlign;
import edu.scripps.yates.utilities.alignment.nwalign.NWResult;
import edu.scripps.yates.utilities.model.enums.AggregationLevel;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;

public class PCQUtils {

	public static final String PROTEIN_DESCRIPTION_SEPARATOR = "####";
	public static DecimalFormat df = new DecimalFormat("#.#");
	private static Map<String, DBIndexInterface> indexByFastaIndexKey = new HashMap<String, DBIndexInterface>();
	private static final Map<String, UniprotProteinLocalRetriever> uplrMap = new HashMap<String, UniprotProteinLocalRetriever>();
	private final static Logger log = Logger.getLogger(PCQUtils.class);
	public static final String PROTEIN_ACC_SEPARATOR = " ";
	public static final double factor = 1.2;
	private static Map<String, QuantParser> quantParsersByFileNamesKey = new HashMap<String, QuantParser>();

	public static Map<String, Set<NWResult>> alignPeptides(List<QuantifiedPeptideInterface> peptideList,
			QuantCondition cond1, QuantCondition cond2, FileWriter Output2) throws IOException {
		Set<NWResult> gAS = new HashSet<NWResult>();
		Map<String, Set<NWResult>> gAM = new HashMap<String, Set<NWResult>>();

		try {

			for (int i = 0; i < peptideList.size(); i++) {
				QuantifiedPeptideInterface pep1 = peptideList.get(i);
				if (i % 25 == 0) {
					log.info("Aligning Peptides " + i + "/" + peptideList.size());
				}

				for (int j = i + 1; j < peptideList.size(); j++) {
					QuantifiedPeptideInterface pep2 = peptideList.get(j);
					NWResult pepResults = NWAlign.needlemanWunsch(pep1.getSequence(), pep2.getSequence(), -11, -1);

					if ((pepResults.getFinalAlignmentScore() >= ProteinClusterQuantParameters.getInstance()
							.getFinalAlignmentScore()
							&& pepResults.getSequenceIdentity() >= ProteinClusterQuantParameters.getInstance()
									.getSequenceIdentity())
							&& pepResults.getMaxConsecutiveIdenticalAlignment() >= ProteinClusterQuantParameters
									.getInstance().getMinConsecutiveIdenticalAlignment()) {
						// pepResults call to function to put result in map
						// pass map, result, call on both pep1
						gAM = putResultInMap(gAM, pepResults, gAS, pep1);
						gAM = putResultInMap(gAM, pepResults, gAS, pep2);

						// Set<QuantifiedProteinInterface> protSet =
						// pep1.getQuantifiedProteins();
						//
						// if ((pepResults.getFinalAlignmentScore() >= 30 &&
						// pepResults.getSequenceIdentity() >= 0.8)
						// && pepResults.getMaxConsecutiveIdenticalAlignment()
						// >= 6)
						// {
						// for (QuantifiedProteinInterface quantifiedProtein :
						// protSet) {
						// print p1 and pep1
						/*
						 * System.out.print(quantifiedProtein.getAccession() +
						 * "\t"); System.out.print(pep1.getSequence() + "\t");
						 * System.out.print("1" + "\t"); System.out.print("0" +
						 * "\t"); log.info("1");
						 */

						// print pep1 and pep2
						Output2.append(pep1.getSequence() + "\t");
						Output2.append(pep2.getSequence() + "\t");
						Double ratio1 = null;
						if (pep1 instanceof IsobaricQuantifiedPeptide) {
							ratio1 = ((IsobaricQuantifiedPeptide) pep1).getCountRatio(cond1, cond2);
						}
						if (ratio1 == Double.NEGATIVE_INFINITY) {
							Output2.append("NEG_INF" + "\t");
						} else if (ratio1 == Double.POSITIVE_INFINITY) {
							Output2.append("POS_INF" + "\t");
						} else {
							Output2.append(ratio1 + "\t");
						}
						Double ratio2 = null;
						if (pep2 instanceof IsobaricQuantifiedPeptide) {
							ratio2 = ((IsobaricQuantifiedPeptide) pep2).getCountRatio(cond1, cond2);
						}
						if (ratio2 == Double.NEGATIVE_INFINITY) {
							Output2.append("NEG_INF" + "\t");
						} else if (ratio2 == Double.POSITIVE_INFINITY) {
							Output2.append("POS_INF" + "\t");
						} else {
							Output2.append(ratio2 + "\t");
						}

						Output2.append(pepResults.getSequenceIdentity() + "\t");
						Output2.append("1" + "\n");
						Output2.flush();

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
			return gAM;
		} finally

		{
			if (Output2 != null) {
				Output2.close();
			}
		}

	}

	private static Map<String, Set<NWResult>> putResultInMap(Map<String, Set<NWResult>> gAM, NWResult result,
			Set<NWResult> gAS, QuantifiedPeptideInterface pep) {
		String seq = pep.getSequence();

		// if gAM has the sequence, use that sequence to get the result
		if ((gAM.containsKey(seq))) {
			gAM.get(seq).add(result);
		}
		// if gAM does not have the sequence, make a new set of results, add the
		// results to this set, and then put the sequence and the set in gAM
		else {
			Set<NWResult> set = new HashSet<NWResult>();
			set.add(result);
			gAM.put(seq, set);
		}
		// return to alignPeptides
		return gAM;
	}

	// merges the clusters if there is a similar peptide pair between the two
	// clusters
	public static ProteinCluster mergeClusters(ProteinCluster cluster, ProteinCluster cluster2) {
		for (QuantifiedPeptideInterface peptide : cluster2.getPeptideSet()) {
			cluster.addIndividualQuantifiedPeptide(peptide);
		}

		for (QuantifiedProteinInterface protein : cluster2.getProteinSet()) {
			cluster.addIndividualQuantifiedProtein(protein);
		}
		return cluster;
	}

	public static CensusChroParser getCensusChroParser(File fastaFile, File inputFileFolder, String[] fileNames,
			Map<QuantCondition, QuantificationLabel> labelsByConditions, char[] enzymeArray, int missedCleavages,
			File uniprotReleasesFolder, String uniprotVersion, String decoyRegexp, boolean ignoreNotFoundPeptidesInDB)
			throws FileNotFoundException {
		List<Map<QuantCondition, QuantificationLabel>> list = new ArrayList<Map<QuantCondition, QuantificationLabel>>();
		for (int i = 0; i < fileNames.length; i++) {
			list.add(labelsByConditions);
		}
		return getCensusChroParser(fastaFile, inputFileFolder, fileNames, list, enzymeArray, missedCleavages,
				uniprotReleasesFolder, uniprotVersion, decoyRegexp, ignoreNotFoundPeptidesInDB);
	}

	public static CensusChroParser getCensusChroParser(File fastaFile, File inputFilefolder, String[] fileNames,
			List<Map<QuantCondition, QuantificationLabel>> labelsByConditions, char[] enzymeArray, int missedCleavages,
			File uniprotReleasesFolder, String uniprotVersion, String decoyRegexp, boolean ignoreNotFoundPeptidesInDB)
			throws FileNotFoundException {
		// Set parser (6 files) to peptides
		List<RemoteSSHFileReference> xmlFiles = new ArrayList<RemoteSSHFileReference>();
		for (String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.add(new RemoteSSHFileReference(inputXmlFile));
		}
		String fileNamesKey = getFileNamesKey(xmlFiles);
		if (quantParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return (CensusChroParser) quantParsersByFileNamesKey.get(fileNamesKey);
		}

		CensusChroParser parser = new CensusChroParser(xmlFiles, labelsByConditions, QuantificationLabel.LIGHT,
				QuantificationLabel.HEAVY);
		try {
			parser.addIonExclusion(IonSerieType.B, 1);
			parser.addIonExclusion(IonSerieType.Y, 1);
			parser.setDecoyPattern(decoyRegexp);
			parser.setDistinguishModifiedPeptides(false);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
			final DBIndexInterface fastaDBIndex = getFastaDBIndex(fastaFile, enzymeArray, missedCleavages);
			parser.setDbIndex(fastaDBIndex);
			parser.enableProteinMergingBySecondaryAccessions(
					getUniprotProteinLocalRetrieverByFolder(uniprotReleasesFolder), uniprotVersion);

			return parser;
		} finally {
			addQuantParserToStaticMap(fileNamesKey, parser);
		}
	}

	private static String getFileNamesKey(List<RemoteSSHFileReference> xmlFiles) {
		StringBuilder sb = new StringBuilder();
		for (RemoteSSHFileReference remoteSSHFileReference : xmlFiles) {
			sb.append(remoteSSHFileReference.getOutputFile().getAbsolutePath());
		}
		return sb.toString();
	}

	private static CensusOutParser getCensusOutParser(File fastaFile, File inputFilefolder, String[] fileNames,
			List<Map<QuantCondition, QuantificationLabel>> labelsByConditions, char[] enzymeArray, int missedCleavages,
			File uniprotReleasesFolder, String uniprotVersion, String decoyRegexp, boolean ignoreNotFoundPeptidesInDB,
			boolean onlyOneSpectrumPerChromatographicPeakAndPerSaltStep, boolean skipSingletons)
			throws FileNotFoundException {
		// Set parser (6 files) to peptides
		List<RemoteSSHFileReference> xmlFiles = new ArrayList<RemoteSSHFileReference>();

		for (String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.add(new RemoteSSHFileReference(inputXmlFile));
		}
		String fileNamesKey = getFileNamesKey(xmlFiles);
		if (quantParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return (CensusOutParser) quantParsersByFileNamesKey.get(fileNamesKey);
		}
		CensusOutParser parser = new CensusOutParser(xmlFiles, labelsByConditions, QuantificationLabel.LIGHT,
				QuantificationLabel.HEAVY);
		try {
			parser.setDecoyPattern(decoyRegexp);
			parser.setDistinguishModifiedPeptides(false);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
			parser.setOnlyOneSpectrumPerChromatographicPeakAndPerSaltStep(
					onlyOneSpectrumPerChromatographicPeakAndPerSaltStep);
			parser.setSkipSingletons(skipSingletons);
			final DBIndexInterface fastaDBIndex = getFastaDBIndex(fastaFile, enzymeArray, missedCleavages);
			parser.setDbIndex(fastaDBIndex);
			parser.enableProteinMergingBySecondaryAccessions(
					getUniprotProteinLocalRetrieverByFolder(uniprotReleasesFolder), uniprotVersion);
			return parser;
		} finally {
			addQuantParserToStaticMap(fileNamesKey, parser);

		}
	}

	private static SeparatedValuesParser getSeparatedValuesParser(File fastaFile, File inputFilefolder,
			String[] fileNames, String separator, List<Map<QuantCondition, QuantificationLabel>> labelsByConditions,
			char[] enzymeArray, int missedCleavages, File uniprotReleasesFolder, String uniprotVersion,
			String decoyRegexp, boolean ignoreNotFoundPeptidesInDB) throws FileNotFoundException {
		// Set parser (6 files) to peptides
		List<RemoteSSHFileReference> xmlFiles = new ArrayList<RemoteSSHFileReference>();

		for (String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.add(new RemoteSSHFileReference(inputXmlFile));
		}
		String fileNamesKey = getFileNamesKey(xmlFiles);
		if (quantParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return (SeparatedValuesParser) quantParsersByFileNamesKey.get(fileNamesKey);
		}
		SeparatedValuesParser parser = new SeparatedValuesParser(xmlFiles, separator, labelsByConditions,
				QuantificationLabel.LIGHT, QuantificationLabel.HEAVY);
		try {
			parser.setDecoyPattern(decoyRegexp);
			parser.setDistinguishModifiedPeptides(false);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
			final DBIndexInterface fastaDBIndex = getFastaDBIndex(fastaFile, enzymeArray, missedCleavages);
			parser.setDbIndex(fastaDBIndex);
			parser.enableProteinMergingBySecondaryAccessions(
					getUniprotProteinLocalRetrieverByFolder(uniprotReleasesFolder), uniprotVersion);
			return parser;
		} finally {
			addQuantParserToStaticMap(fileNamesKey, parser);
		}
	}

	private static void addQuantParserToStaticMap(String fileNamesKey, QuantParser parser) {
		quantParsersByFileNamesKey.put(fileNamesKey, parser);
		log.info(quantParsersByFileNamesKey.size() + " parsers stored.");
	}

	private static DBIndexInterface getFastaDBIndex(File fastaFile, char[] enzymeArray, int missedCleavages) {
		if (fastaFile != null) {

			DBIndexSearchParams defaultDBIndexParams = DBIndexInterface.getDefaultDBIndexParams(fastaFile);
			String fastaIndexKey = IndexUtil.createFullIndexFileName(defaultDBIndexParams);
			if (indexByFastaIndexKey.containsKey(fastaIndexKey)) {
				return indexByFastaIndexKey.get(fastaIndexKey);
			}
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeArr(enzymeArray, missedCleavages, false);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeOffset(0);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeNocutResidues("");
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setH2OPlusProtonAdded(true);
			DBIndexInterface dbIndex = new DBIndexInterface(defaultDBIndexParams);
			indexByFastaIndexKey.put(fastaIndexKey, dbIndex);
			return dbIndex;
		}
		return null;
	}

	private static CensusChroParser getCensusChroParserUsingMongoDBIndex(String mongoDBURI, String mongoMassDBName,
			String mongoSeqDBName, String mongoProtDBName, File inputFilefolder, String[] fileNames,
			List<Map<QuantCondition, QuantificationLabel>> labelsByConditions, File uniprotReleasesFolder,
			String uniprotVersion, String decoyRegexp, boolean ignoreNotFoundPeptidesInDB)
			throws FileNotFoundException {
		// Set parser (6 files) to peptides
		List<RemoteSSHFileReference> xmlFiles = new ArrayList<RemoteSSHFileReference>();

		for (String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.add(new RemoteSSHFileReference(inputXmlFile));
		}
		String fileNamesKey = getFileNamesKey(xmlFiles);
		if (quantParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return (CensusChroParser) quantParsersByFileNamesKey.get(fileNamesKey);
		}
		CensusChroParser parser = new CensusChroParser(xmlFiles, labelsByConditions, QuantificationLabel.LIGHT,
				QuantificationLabel.HEAVY);
		try {
			parser.addIonExclusion(IonSerieType.B, 1);
			parser.addIonExclusion(IonSerieType.Y, 1);
			parser.setDecoyPattern(decoyRegexp);
			parser.setDistinguishModifiedPeptides(false);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
			final DBIndexInterface mongoDBIndex = getMongoDBIndex(mongoDBURI, mongoMassDBName, mongoSeqDBName,
					mongoProtDBName);
			parser.setDbIndex(mongoDBIndex);
			parser.enableProteinMergingBySecondaryAccessions(
					getUniprotProteinLocalRetrieverByFolder(uniprotReleasesFolder), uniprotVersion);

			return parser;
		} finally {
			addQuantParserToStaticMap(fileNamesKey, parser);

		}
	}

	public static UniprotProteinLocalRetriever getUniprotProteinLocalRetrieverByFolder(File uniprotReleasesFolder) {
		if (!uplrMap.containsKey(uniprotReleasesFolder.getAbsolutePath())) {
			uplrMap.put(uniprotReleasesFolder.getAbsolutePath(),
					new UniprotProteinLocalRetriever(uniprotReleasesFolder, true));
		}
		return uplrMap.get(uniprotReleasesFolder.getAbsolutePath());
	}

	private static CensusOutParser getCensusOutParserUsingMongoDBIndex(String mongoDBURI, String mongoMassDBName,
			String mongoSeqDBName, String mongoProtDBName, File inputFilefolder, String[] fileNames,
			List<Map<QuantCondition, QuantificationLabel>> labelsByConditions, File uniprotReleasesFolder,
			String uniprotVersion, String decoyRegexp, boolean ignoreNotFoundPeptidesInDB,
			boolean onlyOneSpectrumPerChromatographicPeakAndPerSaltStep, boolean skipSingletons)
			throws FileNotFoundException {
		// Set parser (6 files) to peptides
		List<RemoteSSHFileReference> xmlFiles = new ArrayList<RemoteSSHFileReference>();

		for (String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.add(new RemoteSSHFileReference(inputXmlFile));
		}
		String fileNamesKey = getFileNamesKey(xmlFiles);
		if (quantParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return (CensusOutParser) quantParsersByFileNamesKey.get(fileNamesKey);
		}
		CensusOutParser parser = new CensusOutParser(xmlFiles, labelsByConditions, QuantificationLabel.LIGHT,
				QuantificationLabel.HEAVY);
		try {
			parser.setDecoyPattern(decoyRegexp);
			parser.setDistinguishModifiedPeptides(false);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
			parser.setOnlyOneSpectrumPerChromatographicPeakAndPerSaltStep(
					onlyOneSpectrumPerChromatographicPeakAndPerSaltStep);
			parser.setSkipSingletons(skipSingletons);
			DBIndexInterface dbIndex = getMongoDBIndex(mongoDBURI, mongoMassDBName, mongoSeqDBName, mongoProtDBName);
			parser.setDbIndex(dbIndex);
			parser.enableProteinMergingBySecondaryAccessions(
					getUniprotProteinLocalRetrieverByFolder(uniprotReleasesFolder), uniprotVersion);
			return parser;
		} finally {
			addQuantParserToStaticMap(fileNamesKey, parser);

		}
	}

	private static SeparatedValuesParser getSeparatedValuesParserUsingMongoDBIndex(String mongoDBURI,
			String mongoMassDBName, String mongoSeqDBName, String mongoProtDBName, File inputFilefolder,
			String[] fileNames, String separator, List<Map<QuantCondition, QuantificationLabel>> labelsByConditions,
			File uniprotReleasesFolder, String uniprotVersion, String decoyRegexp, boolean ignoreNotFoundPeptidesInDB)
			throws FileNotFoundException {
		// Set parser (6 files) to peptides
		List<RemoteSSHFileReference> xmlFiles = new ArrayList<RemoteSSHFileReference>();

		for (String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.add(new RemoteSSHFileReference(inputXmlFile));
		}
		String fileNamesKey = getFileNamesKey(xmlFiles);
		if (quantParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return (SeparatedValuesParser) quantParsersByFileNamesKey.get(fileNamesKey);
		}
		SeparatedValuesParser parser = new SeparatedValuesParser(xmlFiles, separator, labelsByConditions,
				QuantificationLabel.LIGHT, QuantificationLabel.HEAVY);
		try {
			parser.setDecoyPattern(decoyRegexp);
			parser.setDistinguishModifiedPeptides(false);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
			DBIndexInterface dbIndex = getMongoDBIndex(mongoDBURI, mongoMassDBName, mongoSeqDBName, mongoProtDBName);
			parser.setDbIndex(dbIndex);
			parser.enableProteinMergingBySecondaryAccessions(
					getUniprotProteinLocalRetrieverByFolder(uniprotReleasesFolder), uniprotVersion);

			return parser;
		} finally {
			addQuantParserToStaticMap(fileNamesKey, parser);
		}
	}

	private static DBIndexInterface getMongoDBIndex(String mongoDBURI, String mongoMassDBName, String mongoSeqDBName,
			String mongoProtDBName) {

		if (mongoDBURI != null) {
			log.info("Using mongoDB at: " + mongoDBURI + " massDBName:" + mongoMassDBName + " seqDBName:"
					+ mongoSeqDBName + " protDBName:" + mongoProtDBName);

			DBIndexSearchParamsImpl params = new DBIndexSearchParamsImpl(mongoDBURI, mongoMassDBName, mongoSeqDBName,
					mongoProtDBName);
			DBIndexInterface dbIndex = new DBIndexInterface(params);
			return dbIndex;
		}
		return null;
	}

	public static boolean shareAtLeastOnePeptideNode(PCQProteinNode proteinNode1, PCQProteinNode proteinNode2) {

		Set<PCQPeptideNode> peptideNodes1 = proteinNode1.getPeptideNodes();
		for (PCQPeptideNode peptideNode : peptideNodes1) {
			if (peptideNode.isDiscarded()) {
				continue;
			}
			if (proteinNode2.getPeptideNodes().contains(peptideNode)) {
				return true;
			}
		}
		return false;
	}

	public static boolean shareAtLeastOnePeptideBySimilarity(PCQProteinNode proteinNode1, PCQProteinNode proteinNode2,
			Map<String, Set<NWResult>> gAM) {
		Set<PCQPeptideNode> peptideNodes1 = proteinNode1.getPeptideNodes();

		for (PCQPeptideNode peptideNode1 : peptideNodes1) {
			if (peptideNode1.isDiscarded()) {
				continue;
			}
			String peptideSeq1 = peptideNode1.getSequence();
			if (gAM != null && gAM.containsKey(peptideSeq1)) {
				Set<NWResult> alignments = gAM.get(peptideSeq1);

				for (NWResult nwResult : alignments) {
					String peptideSeq2 = nwResult.getSeq1();
					if (peptideSeq1.equals(peptideSeq2)) {
						peptideSeq2 = nwResult.getSeq2();
					}

					if (proteinNodeContainsPeptideSequence(proteinNode2, peptideSeq2)) {
						return true;
					}
				}
			}
		}
		return false;
	}

	private static Boolean proteinNodeContainsPeptideSequence(PCQProteinNode proteinNode, String peptideSeq) {
		Set<QuantifiedPeptideInterface> peptides = proteinNode.getQuantifiedPeptides();
		for (QuantifiedPeptideInterface peptide : peptides) {
			if (peptide.getSequence().equals(peptideSeq)) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Get the set of shared {@link PCQPeptideNode} of the
	 * {@link PCQProteinNode} and get the consensus ratio, calling to
	 * {@link PCQUtils}.
	 * {@link #getConsensusRatio(peptideNodes, condition1, condition2)}
	 *
	 * @param proteinNode1
	 * @param proteinNode2
	 * @param cond1
	 * @param cond2
	 * @return
	 */
	public static QuantRatio getSharedConsensusPeptideRatio(PCQProteinNode proteinNode1, PCQProteinNode proteinNode2,
			QuantCondition cond1, QuantCondition cond2) {
		final Set<PCQPeptideNode> sharedPeptides = getSharedPeptideNodeSet(proteinNode1, proteinNode2, false);
		return getConsensusRatio(sharedPeptides, cond1, cond2, null);
	}

	/**
	 * Gets a consensus {@link QuantRatio} from a set of {@link PCQPeptideNode}.
	 * If the sets contains only one {@link PCQPeptideNode}, the consensus
	 * {@link QuantRatio} will be the one from that peptide. If there is more
	 * than one {@link PCQPeptideNode}, then an average {@link QuantRatio} will
	 * be returned.
	 *
	 * @param peptideNodes
	 * @param cond1
	 * @param cond2
	 * @param replicateName
	 *            if not null, get the consensus ratio calling to
	 *            {@link #getConsensusRatio(Collection, QuantCondition, QuantCondition, replicateName)}
	 * @return
	 */
	public static QuantRatio getConsensusRatio(Collection<PCQPeptideNode> peptideNodes, QuantCondition cond1,
			QuantCondition cond2, String replicateName) {
		// if there is only one peptide, get the consensus ratio from it
		if (peptideNodes.size() == 1) {
			final PCQPeptideNode peptideNode = peptideNodes.iterator().next();
			if (!peptideNode.isDiscarded()) {
				return peptideNode.getConsensusRatio(cond1, cond2, replicateName);
			}
			return null;
		}

		// otherwise calculate an average over the non infinity ratios
		final Set<QuantRatio> consensusRatios = getConsensusRatios(peptideNodes, cond1, cond2, replicateName);
		Set<QuantRatio> nonInfinityRatios = QuantUtils.getNonInfinityRatios(consensusRatios);
		return QuantUtils.getAverageRatio(nonInfinityRatios, AggregationLevel.PEPTIDE);

	}

	/**
	 * Gets a {@link Set} of individual {@link QuantRatio} as consensus ratios
	 * for each of the {@link PCQPeptideNode} in the {@link Collection}, calling
	 * to getConsensusRatio() in each {@link PCQPeptideNode}
	 *
	 * @param peptideNodes
	 * @param cond1
	 * @param cond2
	 * @param replicateName
	 *            if not null, get only the consensus {@link QuantRatio} for
	 *            that replicate
	 * @return
	 */
	private static Set<QuantRatio> getConsensusRatios(Collection<PCQPeptideNode> peptideNodes, QuantCondition cond1,
			QuantCondition cond2, String replicateName) {
		Set<QuantRatio> ret = new HashSet<QuantRatio>();
		for (PCQPeptideNode peptideNode : peptideNodes) {
			if (peptideNode.isDiscarded()) {
				continue;
			}
			if (replicateName != null) {
				final QuantRatio consensusRatio = peptideNode.getConsensusRatio(cond1, cond2, replicateName);
				if (consensusRatio != null) {
					ret.add(consensusRatio);
				}
			} else {
				final QuantRatio consensusRatio = peptideNode.getConsensusRatio(cond1, cond2);
				if (consensusRatio != null) {
					ret.add(consensusRatio);
				}
			}
		}
		return ret;
	}

	public static Set<PCQPeptideNode> getSharedPeptideNodeSet(PCQProteinNode proteinNode1, PCQProteinNode proteinNode2,
			boolean onlysharedByThisToProteins) {
		Set<PCQPeptideNode> peptideNodes1 = proteinNode1.getPeptideNodes();
		Set<PCQPeptideNode> peptideNodes2 = proteinNode2.getPeptideNodes();
		Set<PCQPeptideNode> ret = new HashSet<PCQPeptideNode>();
		for (PCQPeptideNode peptideNode1 : peptideNodes1) {
			if (peptideNode1.isDiscarded()) {
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
			PCQProteinNode proteinNode2, QuantCondition cond1, QuantCondition cond2) {
		Set<QuantifiedPeptideInterface> peptides1 = proteinNode1.getQuantifiedPeptides();
		Set<QuantifiedPeptideInterface> peptides2 = proteinNode2.getQuantifiedPeptides();
		List<Double> ratioValues = new ArrayList<Double>();

		for (QuantifiedPeptideInterface peptide1 : peptides1) {
			if (peptides2.contains(peptide1)) {
				final QuantRatio consensusRatio = peptide1.getConsensusRatio(cond1, cond2);
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
			PCQProteinNode proteinNode2, boolean onlySharedByTheseTwoProteins) {

		Set<PCQPeptideNode> peptidesNodes1 = proteinNode1.getPeptideNodes();
		Set<PCQPeptideNode> peptidesNodes2 = proteinNode2.getPeptideNodes();
		Set<PCQPeptideNode> totalPeptideNodes = new HashSet<PCQPeptideNode>();
		totalPeptideNodes.addAll(peptidesNodes1);
		totalPeptideNodes.addAll(peptidesNodes2);
		Map<String, Set<PCQPeptideNode>> map = new HashMap<String, Set<PCQPeptideNode>>();
		for (PCQPeptideNode peptideNode : totalPeptideNodes) {
			if (peptideNode.isDiscarded()) {
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
					for (PCQProteinNode proteinNode : proteinNodes) {
						if (proteinNode.isDiscarded()) {
							continue;
						}
						String proteinAccKey = proteinNode.getAccession();
						if (map.containsKey(proteinAccKey)) {
							map.get(proteinAccKey).add(peptideNode);
						} else {
							Set<PCQPeptideNode> set = new HashSet<PCQPeptideNode>();
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

		Map<String, Set<PCQPeptideNode>> map = new HashMap<String, Set<PCQPeptideNode>>();
		if (proteinNode1 == null || proteinNode2 == null) {
			return map;
		}
		Set<PCQPeptideNode> peptideNodes1 = proteinNode1.getPeptideNodes();
		Set<PCQPeptideNode> peptideNodes2 = proteinNode2.getPeptideNodes();
		Set<PCQPeptideNode> totalPeptideNodes = new HashSet<PCQPeptideNode>();
		totalPeptideNodes.addAll(peptideNodes1);
		totalPeptideNodes.addAll(peptideNodes2);

		for (PCQPeptideNode peptideNode : totalPeptideNodes) {
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
					String proteinAccKey = PCQUtils.getProteinNodeAccessionString(proteinNodes);
					if (map.containsKey(proteinAccKey)) {
						map.get(proteinAccKey).add(peptideNode);
					} else {
						Set<PCQPeptideNode> set = new HashSet<PCQPeptideNode>();
						set.add(peptideNode);
						map.put(proteinAccKey, set);
					}
				}
			}
		}
		return map;
	}

	/**
	 * Gets the consensus {@link QuantRatio} for the peptides unique to protein1
	 * in respect to protein2
	 *
	 * @param proteinNode1
	 * @param proteinNode2
	 * @param cond1
	 * @param cond2
	 * @param uniquePepOnly
	 *            if true, any peptide unique to protein1 in respect to protein2
	 *            that belongs to another protein will not be included.
	 * @return
	 */
	public static QuantRatio getUniqueConsensusPeptideRatio(PCQProteinNode proteinNode1, PCQProteinNode proteinNode2,
			QuantCondition cond1, QuantCondition cond2, boolean uniquePepOnly) {

		final Set<PCQPeptideNode> uniquePeptideNodes = getUniquePeptideNodes(proteinNode1, proteinNode2, uniquePepOnly,
				true);
		if (uniquePeptideNodes.isEmpty()) {
			return null;
		}
		// in case of isobaric ratios, calculate a consensus isobaric ratio
		if (uniquePeptideNodes.iterator().next().getQuantifiedPeptides().iterator()
				.next() instanceof IsobaricQuantifiedPeptide) {
			Set<IsobaricQuantifiedPeptide> isobaricQuantPeptides = new HashSet<IsobaricQuantifiedPeptide>();
			for (PCQPeptideNode peptideNode : uniquePeptideNodes) {
				if (peptideNode.isDiscarded()) {
					continue;
				}
				for (QuantifiedPeptideInterface peptide : peptideNode.getQuantifiedPeptides()) {
					if (peptide.isDiscarded()) {
						continue;
					}
					isobaricQuantPeptides.add((IsobaricQuantifiedPeptide) peptide);
				}

			}
			return getConsensusIonCountRatio(isobaricQuantPeptides, cond1, cond2);
		} else {
			return getConsensusRatio(uniquePeptideNodes, cond1, cond2, null);
		}

	}

	/**
	 * Gets a consensus {@link IonCountRatio} from a set of
	 * {@link IsobaricQuantifiedPeptide} where the ions from each
	 * {@link QuantCondition} are pulled together and normalized by the number
	 * of {@link QuantifiedPSMInterface} per {@link QuantifiedPeptideInterface}
	 *
	 * @param isobaricQuantifiedPeptides
	 * @param cond1
	 * @param cond2
	 * @return
	 */
	public static IonCountRatio getConsensusIonCountRatio(Set<IsobaricQuantifiedPeptide> isobaricQuantifiedPeptides,
			QuantCondition cond1, QuantCondition cond2) {

		IonCountRatio ratio = new IonCountRatio(AggregationLevel.PEPTIDE);
		for (IsobaricQuantifiedPeptide isoPeptide : isobaricQuantifiedPeptides) {
			final int numPSMs = isoPeptide.getQuantifiedPSMs().size();
			// get number of ions in one condition, and normalize by the
			// number of PSMs
			int peakCount1 = 0;
			if (isoPeptide.getIonsByCondition().containsKey(cond1)) {
				peakCount1 = isoPeptide.getIonsByCondition().get(cond1).size();
			}
			double normalizedPeakCount1 = peakCount1 * 1.0 / numPSMs;
			int peakCount2 = 0;
			if (isoPeptide.getIonsByCondition().containsKey(cond2)) {
				peakCount2 = isoPeptide.getIonsByCondition().get(cond2).size();
			}
			double normalizedPeakCount2 = peakCount2 * 1.0 / numPSMs;
			ratio.addIonCount(cond1, normalizedPeakCount1);
			ratio.addIonCount(cond2, normalizedPeakCount2);

		}
		return ratio;
	}

	/**
	 * explainwhatis doing
	 *
	 * @param proteinNode1
	 * @param proteinNode2
	 * @param cond1
	 * @param cond2
	 * @return
	 */
	public static List<Double> getUniqueIndividualPeptideRatioValues(PCQProteinNode proteinNode1,
			PCQProteinNode proteinNode2, QuantCondition cond1, QuantCondition cond2, boolean uniquePepOnly) {
		Set<PCQPeptideNode> peptideNodes1 = getUniquePeptideNodes(proteinNode1, proteinNode2, uniquePepOnly, true);
		List<Double> ratios = new ArrayList<Double>();

		for (PCQPeptideNode peptideNode1 : peptideNodes1) {
			if (peptideNode1.isDiscarded()) {
				continue;
			}
			for (QuantifiedPeptideInterface peptide : peptideNode1.getQuantifiedPeptides()) {
				final QuantRatio consensusRatio = peptide.getConsensusRatio(cond1, cond2);
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
		for (Double ratio : ratios) {
			if (ratio != null && Double.compare(ratio, Double.NaN) != 0) {
				if (Double.isInfinite(ratio)) {
					counterINF++;
				} else {
					counterRatio++;
				}
			} else {
				// means that one of the threeCompare is null!
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
		List<Double> INFList = new ArrayList<Double>();
		for (Double ratio : ratios) {
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
		List<Double> RatioList = new ArrayList<Double>();
		for (Double ratio : ratios) {
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
		for (Double ratio : ratios) {
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
		for (IonCountRatio ratio : ratios) {
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

		Set<IsobaricQuantifiedPSM> psms = new HashSet<IsobaricQuantifiedPSM>();
		for (QuantifiedPSMInterface psm : peptideNode.getQuantifiedPSMs()) {
			if (psm instanceof IsobaricQuantifiedPSM) {
				psms.add((IsobaricQuantifiedPSM) psm);
			}
		}
		return getIonCountFromPSMs(psms, condition);

	}

	private static int getIonCountFromPSMs(Collection<IsobaricQuantifiedPSM> psms, QuantCondition condition) {
		int total = 0;
		for (IsobaricQuantifiedPSM quantifiedPSM : psms) {
			Set<Ion> ions = quantifiedPSM.getIonsByCondition().get(condition);
			if (ions != null) {
				total = total + ions.size();
			}
		}
		return total;
	}

	public static List<QuantifiedPeptideInterface> getSortedPeptidesBySequence(
			Collection<QuantifiedPeptideInterface> peptides) {
		List<QuantifiedPeptideInterface> ret = new ArrayList<QuantifiedPeptideInterface>();
		ret.addAll(peptides);
		Collections.sort(ret, new Comparator<QuantifiedPeptideInterface>() {

			@Override
			public int compare(QuantifiedPeptideInterface o1, QuantifiedPeptideInterface o2) {
				return o1.getSequence().compareTo(o2.getSequence());
			}
		});
		return ret;
	}

	public static List<PCQPeptideNode> getSortedPeptideNodesBySequence(Collection<PCQPeptideNode> peptides) {
		List<PCQPeptideNode> ret = new ArrayList<PCQPeptideNode>();
		ret.addAll(peptides);
		Collections.sort(ret, new Comparator<PCQPeptideNode>() {

			@Override
			public int compare(PCQPeptideNode o1, PCQPeptideNode o2) {
				return o1.getSequence().compareTo(o2.getSequence());
			}
		});
		return ret;
	}

	public static List<PCQProteinNode> getSortedProteinNodesByAcc(Collection<PCQProteinNode> proteinNodesToSort) {
		List<PCQProteinNode> list = new ArrayList<PCQProteinNode>();
		list.addAll(proteinNodesToSort);
		Collections.sort(list, new Comparator<PCQProteinNode>() {
			@Override
			public int compare(PCQProteinNode o1, PCQProteinNode o2) {
				return o1.getAccession().compareTo(o2.getAccession());
			}
		});
		return list;
	}

	/**
	 * Get a CVS list of protein accessions after sorting them alphabetically
	 *
	 * @param proteins
	 * @return
	 */
	public static String getAccessionString(Collection<QuantifiedProteinInterface> proteins) {
		Set<String> set = new HashSet<String>();
		for (QuantifiedProteinInterface protein : proteins) {
			set.add(protein.getAccession());
		}
		List<String> list = new ArrayList<String>();
		list.addAll(set);
		Collections.sort(list);

		StringBuilder sb = new StringBuilder();
		for (String acc : list) {
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

		StringBuilder sb = new StringBuilder();
		for (PCQProteinNode proteinNode : PCQUtils.getSortedProteinNodesByAcc(proteinNodes)) {
			if (!"".equals(sb.toString()))
				sb.append(PROTEIN_ACC_SEPARATOR);
			sb.append(proteinNode.getAccession());
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

		StringBuilder sb = new StringBuilder();
		for (QuantifiedPeptideInterface peptide : PCQUtils.getSortedPeptidesBySequence(peptides)) {
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

		StringBuilder sb = new StringBuilder();
		for (PCQPeptideNode peptideNode : PCQUtils.getSortedPeptideNodesBySequence(peptideNodes)) {
			if (!"".equals(sb.toString()))
				sb.append("_");
			sb.append(peptideNode.getSequence());
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
		StringBuilder sb = new StringBuilder();
		Set<String> descriptions = new HashSet<String>();
		for (QuantifiedProteinInterface protein : QuantUtils.getSortedQuantifiedProteinsByAcc(proteins)) {
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
	public static Set<PCQPeptideNode> getUniquePeptideNodes(PCQProteinNode proteinNode) {
		Set<PCQPeptideNode> ret = new HashSet<PCQPeptideNode>();
		for (PCQPeptideNode peptideNode : proteinNode.getPeptideNodes()) {
			if (peptideNode.isDiscarded()) {
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
			Set<String> validTaxonomies, boolean onlyFirst) {
		return getGeneNameString(annotatedProteins, getProteinMap(cluster.getProteinNodes()), validTaxonomies,
				onlyFirst);
	}

	public static String getGeneNameString(Map<String, Entry> annotatedProteins,
			Collection<PCQProteinNode> proteinNodes, Set<String> validTaxonomies, boolean onlyFirst) {
		return getGeneNameString(annotatedProteins, getProteinMap(proteinNodes), validTaxonomies, onlyFirst);
	}

	public static String getGeneNameString(Map<String, Entry> annotatedProteins, PCQProteinNode proteinNode,
			Set<String> validTaxonomies, boolean onlyFirst) {
		return getGeneNameString(annotatedProteins, getProteinMap(proteinNode), validTaxonomies, onlyFirst);
	}

	private static Map<String, PCQProteinNode> getProteinMap(PCQProteinNode protein) {
		Set<PCQProteinNode> set = new HashSet<PCQProteinNode>();
		set.add(protein);
		return getProteinMap(set);
	}

	private static Map<String, PCQProteinNode> getProteinMap(Collection<PCQProteinNode> proteins) {
		Map<String, PCQProteinNode> map = new HashMap<String, PCQProteinNode>();
		for (PCQProteinNode proteinNode : proteins) {
			if (proteinNode.isDiscarded()) {
				continue;
			}
			final String accession = proteinNode.getAccession();
			if (accession.contains(PROTEIN_ACC_SEPARATOR)) {
				final String[] split = accession.split(PROTEIN_ACC_SEPARATOR);
				for (String string : split) {
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
			Set<String> validTaxonomies) {

		Set<String> set = new HashSet<String>();
		for (PCQProteinNode proteinNode : proteinNodeSet) {
			if (proteinNode.isDiscarded()) {
				continue;
			}
			String rawAcc = proteinNode.getAccession();
			List<String> accs = new ArrayList<String>();
			if (rawAcc.contains(" ")) {
				final String[] split = rawAcc.split(" ");
				for (String string : split) {
					accs.add(string);
				}
			} else {
				accs.add(rawAcc);
			}

			// int index = 0;
			for (String acc : accs) {
				if (annotatedProteins != null && annotatedProteins.containsKey(acc)) {
					String taxon = getTaxonomy(acc, annotatedProteins.get(acc).getOrganism());
					if (taxon != null) {
						boolean valid = false;
						if (validTaxonomies != null && !validTaxonomies.isEmpty()) {
							for (String skipTaxonomy : validTaxonomies) {
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
		List<String> list = new ArrayList<String>();
		list.addAll(set);
		Collections.sort(list);

		StringBuilder sb = new StringBuilder();

		for (String string : list) {
			if (!"".equals(sb.toString()))
				sb.append(",");
			sb.append(string);
		}

		return sb.toString();
	}

	public static String getGeneNameString(Map<String, Entry> annotatedProteins,
			Map<String, PCQProteinNode> proteinNodeMap, Set<String> validTaxonomies, boolean onlyFirst) {

		Set<String> set = new HashSet<String>();
		for (String rawAcc : proteinNodeMap.keySet()) {

			List<String> accs = new ArrayList<String>();
			if (rawAcc.contains(" ")) {
				final String[] split = rawAcc.split(" ");
				for (String string : split) {
					accs.add(string);
				}
			} else {
				accs.add(rawAcc);
			}

			// int index = 0;
			for (String acc : accs) {
				if (annotatedProteins != null && annotatedProteins.containsKey(acc)) {
					String taxon = getTaxonomy(acc, annotatedProteins.get(acc).getOrganism());
					if (taxon != null) {
						boolean valid = false;
						if (validTaxonomies != null && !validTaxonomies.isEmpty()) {
							for (String skipTaxonomy : validTaxonomies) {
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
					}

					final String geneName = getGeneName(annotatedProteins.get(acc).getGene());
					set.add(geneName);
				} else {
					// log.warn(acc + " not annotated");
				}
			}
		}
		List<String> list = new ArrayList<String>();
		list.addAll(set);
		Collections.sort(list);
		if (onlyFirst && !list.isEmpty()) {
			return list.get(0);
		}
		StringBuilder sb = new StringBuilder();
		if (!ProteinClusterQuantParameters.getInstance().isPrintOnlyFirstGene()) {
			for (String string : list) {
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
				for (OrganismNameType organismType : organism.getName()) {
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
		for (GeneType geneType : gene) {
			for (GeneNameType geneName : geneType.getName()) {
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
			for (QuantifiedPSMInterface psm : quantifiedPSMs) {
				for (QuantificationLabel label : QuantificationLabel.values()) {
					ionCount += QuantUtils.getIonCount(psm, label);
				}
			}
		}
		return ionCount;
	}

	public static Map<String, QuantifiedPeptideInterface> getPeptideMapFromClusters(Set<ProteinCluster> clusterSet) {
		Map<String, QuantifiedPeptideInterface> map = new HashMap<String, QuantifiedPeptideInterface>();
		for (ProteinCluster proteinCluster : clusterSet) {
			final Set<QuantifiedPeptideInterface> peptideSet = proteinCluster.getPeptideSet();
			for (QuantifiedPeptideInterface quantifiedPeptide : peptideSet) {
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
		List<String> ret = new ArrayList<String>();
		for (QuantifiedProteinInterface protein : proteinSet) {
			final Set<String> taxonomies = protein.getTaxonomies();
			for (String taxonomy : taxonomies) {
				if (!ret.contains(taxonomy)) {
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

		Set<String> proteinAccs1 = PCQUtils.getAccessions(proteins1);
		Set<String> proteinAccs2 = PCQUtils.getAccessions(proteins2);
		if (proteinAccs1.size() == proteinAccs2.size()) {
			for (String acc1 : proteinAccs1) {
				if (!proteinAccs2.contains(acc1)) {
					return false;
				}
			}
			return true;
		}
		return false;
	}

	private static Set<String> getAccessions(Set<QuantifiedProteinInterface> proteins) {
		Set<String> accs = new HashSet<String>();
		for (QuantifiedProteinInterface protein : proteins) {
			accs.add(protein.getAccession());
		}
		return accs;
	}

	public static PCQPeptideNode mergePeptideNodes(PCQPeptideNode peptideNode, PCQPeptideNode peptideNode2) {
		if (peptideNode.hashCode() == peptideNode2.hashCode()) {
			return peptideNode;
		}
		Set<QuantifiedPeptideInterface> peptideCollection = new HashSet<QuantifiedPeptideInterface>();
		peptideCollection.addAll(peptideNode.getQuantifiedPeptides());
		peptideCollection.addAll(peptideNode2.getQuantifiedPeptides());
		PCQPeptideNode ret = new PCQPeptideNode(peptideNode.getCluster(), peptideCollection);
		return ret;
	}

	public static void mergeProteinNodes(PCQProteinNode proteinNode, PCQProteinNode proteinNode2) {
		if (proteinNode.hashCode() == proteinNode2.hashCode()) {
			return;
		}
		final Set<QuantifiedProteinInterface> proteins2 = proteinNode2.getQuantifiedProteins();
		for (QuantifiedProteinInterface quantifiedProteinInterface : proteins2) {
			proteinNode.addProtein(quantifiedProteinInterface);
		}
	}

	public static QuantParser getQuantParser(ProteinClusterQuantParameters params,
			List<Map<QuantCondition, QuantificationLabel>> labelsByConditionsList, final String[] inputFileNamesArray)
			throws FileNotFoundException {
		log.debug("Getting input file parser");
		if (params.getInputType() == AnalysisInputType.CENSUS_CHRO) {
			if (params.getMongoDBURI() != null) {
				return PCQUtils.getCensusChroParserUsingMongoDBIndex(params.getMongoDBURI(),
						params.getMongoMassDBName(), params.getMongoSeqDBName(), params.getMongoProtDBName(),
						params.getInputFileFolder(), inputFileNamesArray, labelsByConditionsList,
						params.getUniprotReleasesFolder(), params.getUniprotVersion(), params.getDecoyRegexp(),
						params.isIgnoreNotFoundPeptidesInDB());
			} else {
				return PCQUtils.getCensusChroParser(params.getFastaFile(), params.getInputFileFolder(),
						inputFileNamesArray, labelsByConditionsList, params.getEnzymeArray(),
						params.getMissedCleavages(), params.getUniprotReleasesFolder(), params.getUniprotVersion(),
						params.getDecoyRegexp(), params.isIgnoreNotFoundPeptidesInDB());
			}
		} else if (params.getInputType() == AnalysisInputType.CENSUS_OUT) {
			if (params.getMongoDBURI() != null) {
				final CensusOutParser parser = PCQUtils.getCensusOutParserUsingMongoDBIndex(params.getMongoDBURI(),
						params.getMongoMassDBName(), params.getMongoSeqDBName(), params.getMongoProtDBName(),
						params.getInputFileFolder(), inputFileNamesArray, labelsByConditionsList,
						params.getUniprotReleasesFolder(), params.getUniprotVersion(), params.getDecoyRegexp(),
						params.isIgnoreNotFoundPeptidesInDB(),
						params.isOnlyOneSpectrumPerChromatographicPeakAndPerSaltStep(), params.isSkipSingletons());

				return parser;
			} else {
				final CensusOutParser parser = PCQUtils.getCensusOutParser(params.getFastaFile(),
						params.getInputFileFolder(), inputFileNamesArray, labelsByConditionsList,
						params.getEnzymeArray(), params.getMissedCleavages(), params.getUniprotReleasesFolder(),
						params.getUniprotVersion(), params.getDecoyRegexp(), params.isIgnoreNotFoundPeptidesInDB(),
						params.isOnlyOneSpectrumPerChromatographicPeakAndPerSaltStep(), params.isSkipSingletons());

				return parser;
			}
		} else if (params.getInputType() == AnalysisInputType.SEPARATED_VALUES) {
			if (params.getMongoDBURI() != null) {
				final SeparatedValuesParser parser = PCQUtils.getSeparatedValuesParserUsingMongoDBIndex(
						params.getMongoDBURI(), params.getMongoMassDBName(), params.getMongoSeqDBName(),
						params.getMongoProtDBName(), params.getInputFileFolder(), inputFileNamesArray,
						params.getSeparator(), labelsByConditionsList, params.getUniprotReleasesFolder(),
						params.getUniprotVersion(), params.getDecoyRegexp(), params.isIgnoreNotFoundPeptidesInDB());

				return parser;
			} else {
				final SeparatedValuesParser parser = PCQUtils.getSeparatedValuesParser(params.getFastaFile(),
						params.getInputFileFolder(), inputFileNamesArray, params.getSeparator(), labelsByConditionsList,
						params.getEnzymeArray(), params.getMissedCleavages(), params.getUniprotReleasesFolder(),
						params.getUniprotVersion(), params.getDecoyRegexp(), params.isIgnoreNotFoundPeptidesInDB());

				return parser;
			}
		}
		throw new IllegalArgumentException("inputType is not recognized");
	}

	public static QuantParser getQuantParser(ProteinClusterQuantParameters params,
			Map<QuantCondition, QuantificationLabel> labelsByConditions, String inputFileName)
			throws FileNotFoundException {
		log.debug("Getting input file parser");
		List<Map<QuantCondition, QuantificationLabel>> labelsByConditionsList = new ArrayList<Map<QuantCondition, QuantificationLabel>>();
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
	public static Double getRatioValue(QuantRatio ratio, QuantCondition cond1, QuantCondition cond2) {
		Double ratioValue = Double.NaN;
		if (ratio != null) {

			Double log2Ratio = ratio.getLog2Ratio(cond1, cond2);
			if (log2Ratio != null) {
				ratioValue = log2Ratio;
			}

		}
		return ratioValue;
	}

	public static Map<String, Set<QuantifiedProteinInterface>> getIndividualProteinsMap(PCQPeptideNode peptideNode) {
		Map<String, Set<QuantifiedProteinInterface>> ret = new HashMap<String, Set<QuantifiedProteinInterface>>();
		for (QuantifiedPeptideInterface peptide : peptideNode.getQuantifiedPeptides()) {
			final Set<QuantifiedProteinInterface> quantifiedProteins = peptide.getQuantifiedProteins();
			for (QuantifiedProteinInterface quantifiedProteinInterface : quantifiedProteins) {
				final String accession = quantifiedProteinInterface.getAccession();
				if (ret.containsKey(accession)) {
					ret.get(accession).add(quantifiedProteinInterface);
				} else {
					Set<QuantifiedProteinInterface> set = new HashSet<QuantifiedProteinInterface>();
					set.add(quantifiedProteinInterface);
					ret.put(accession, set);
				}
			}
		}
		return ret;
	}

	public static Set<PCQPeptideNode> getUniquePeptideNodes(PCQProteinNode proteinNode1, PCQProteinNode proteinNode2,
			boolean uniquePepOnly, boolean skipDiscarded) {
		Set<PCQPeptideNode> peptideNodes1 = proteinNode1.getPeptideNodes();
		Set<PCQPeptideNode> ret = new HashSet<PCQPeptideNode>();

		for (PCQPeptideNode peptideNode1 : peptideNodes1) {
			if (skipDiscarded && peptideNode1.isDiscarded()) {
				continue;
			}
			if (uniquePepOnly) {
				if (peptideNode1.getProteinNodes().size() == 1) {
					ret.add(peptideNode1);
				} else {
					continue;
				}
			} else {
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
		Set<QuantifiedPeptideInterface> peptides1 = new HashSet<QuantifiedPeptideInterface>();
		for (QuantifiedProteinInterface protein1 : proteins1) {
			peptides1.addAll(protein1.getQuantifiedPeptides());
		}
		Set<QuantifiedPeptideInterface> peptides2 = new HashSet<QuantifiedPeptideInterface>();
		for (QuantifiedProteinInterface protein2 : proteins2) {
			peptides2.addAll(protein2.getQuantifiedPeptides());
		}
		if (peptides1.size() == peptides2.size()) {
			for (QuantifiedPeptideInterface peptide2 : peptides2) {
				if (!peptides1.contains(peptide2)) {
					return false;
				}
			}
			return true;
		}
		return false;
	}

}
