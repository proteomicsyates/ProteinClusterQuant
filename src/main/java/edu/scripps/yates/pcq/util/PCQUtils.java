
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
import edu.scripps.yates.census.read.model.CensusRatio;
import edu.scripps.yates.census.read.model.Ion;
import edu.scripps.yates.census.read.model.IonCountRatio;
import edu.scripps.yates.census.read.model.IonSerie.IonSerieType;
import edu.scripps.yates.census.read.model.IsoRatio;
import edu.scripps.yates.census.read.model.IsobaricQuantifiedPSM;
import edu.scripps.yates.census.read.model.IsobaricQuantifiedPeptide;
import edu.scripps.yates.census.read.model.QuantifiedPSM;
import edu.scripps.yates.census.read.model.QuantifiedPeptide;
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
import edu.scripps.yates.pcq.ProteinClusterQuantParameters;
import edu.scripps.yates.pcq.model.PCQPeptideNode;
import edu.scripps.yates.pcq.model.PCQProteinNode;
import edu.scripps.yates.pcq.model.ProteinCluster;
import edu.scripps.yates.pcq.xgmml.util.AlignedPeptides;
import edu.scripps.yates.pcq.xgmml.util.AlignmentSet;
import edu.scripps.yates.utilities.alignment.nwalign.NWAlign;
import edu.scripps.yates.utilities.alignment.nwalign.NWResult;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.model.enums.AggregationLevel;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;
import edu.scripps.yates.utilities.util.StringPosition;

public class PCQUtils {

	public static final String PROTEIN_DESCRIPTION_SEPARATOR = "####";
	public static DecimalFormat df = new DecimalFormat("#.#");
	private static Map<String, DBIndexInterface> indexByFastaIndexKey = new HashMap<String, DBIndexInterface>();
	private static final Map<String, UniprotProteinLocalRetriever> uplrMap = new HashMap<String, UniprotProteinLocalRetriever>();
	private final static Logger log = Logger.getLogger(PCQUtils.class);
	public static final String PROTEIN_ACC_SEPARATOR = " ";
	public static final double factor = 1.2;
	public static final String FDR_CONFIDENCE_SCORE_NAME = "FDR";
	private static Map<String, QuantParser> quantParsersByFileNamesKey = new HashMap<String, QuantParser>();
	private static Map<String, DTASelectParser> dtaSelectParsersByFileNamesKey = new HashMap<String, DTASelectParser>();

	public static AlignmentSet alignPeptides(List<QuantifiedPeptideInterface> peptideList, QuantCondition cond1,
			QuantCondition cond2, FileWriter Output2) throws IOException {
		AlignmentSet ret = new AlignmentSet();

		try {
			log.info("Aligning " + peptideList.size() + " peptides between them");
			for (int i = 0; i < peptideList.size(); i++) {
				QuantifiedPeptideInterface pep1 = peptideList.get(i);
				if (i % 25 == 0) {
					log.info("Aligning Peptides " + i + "/" + peptideList.size());
				}

				for (int j = i + 1; j < peptideList.size(); j++) {
					QuantifiedPeptideInterface pep2 = peptideList.get(j);
					NWResult alignment = NWAlign.needlemanWunsch(pep1.getSequence(), pep2.getSequence(), -11, -1);

					if ((alignment.getFinalAlignmentScore() >= ProteinClusterQuantParameters.getInstance()
							.getFinalAlignmentScore()
							&& alignment.getSequenceIdentity() >= ProteinClusterQuantParameters.getInstance()
									.getSequenceIdentity())
							&& alignment.getMaxConsecutiveIdenticalAlignment() >= ProteinClusterQuantParameters
									.getInstance().getMinConsecutiveIdenticalAlignment()) {
						// store aligment
						ret.addAlignment(new AlignedPeptides(alignment, pep1, pep2));
						// print pep1 and pep2
						Output2.append(pep1.getSequence() + "\t");
						Output2.append(pep2.getSequence() + "\t");
						Double ratio1 = null;
						if (pep1 instanceof IsobaricQuantifiedPeptide) {
							ratio1 = ((IsobaricQuantifiedPeptide) pep1).getIonCountRatio(cond1, cond2)
									.getLog2Ratio(cond1, cond2);
						} else if (pep1 instanceof QuantifiedPeptide) {
							ratio1 = ((QuantifiedPeptide) pep1).getConsensusRatio(cond1, cond2).getLog2Ratio(cond1,
									cond2);
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
							ratio2 = ((IsobaricQuantifiedPeptide) pep2).getIonCountRatio(cond1, cond2)
									.getLog2Ratio(cond1, cond2);
						} else if (pep2 instanceof QuantifiedPeptide) {
							ratio2 = ((QuantifiedPeptide) pep2).getConsensusRatio(cond1, cond2).getLog2Ratio(cond1,
									cond2);
						}
						if (ratio2 == Double.NEGATIVE_INFINITY) {
							Output2.append("NEG_INF" + "\t");
						} else if (ratio2 == Double.POSITIVE_INFINITY) {
							Output2.append("POS_INF" + "\t");
						} else {
							Output2.append(ratio2 + "\t");
						}

						Output2.append(alignment.getSequenceIdentity() + "\t");
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
			return ret;
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
			boolean semiCleavage, File uniprotReleasesFolder, String uniprotVersion, String decoyRegexp,
			boolean ignoreNotFoundPeptidesInDB, boolean distinguishModifiedPeptides, String peptideFilterRegexp)
			throws FileNotFoundException {
		List<Map<QuantCondition, QuantificationLabel>> list = new ArrayList<Map<QuantCondition, QuantificationLabel>>();
		for (int i = 0; i < fileNames.length; i++) {
			list.add(labelsByConditions);
		}
		return getCensusChroParser(fastaFile, inputFileFolder, fileNames, list, enzymeArray, missedCleavages,
				semiCleavage, uniprotReleasesFolder, uniprotVersion, decoyRegexp, ignoreNotFoundPeptidesInDB,
				distinguishModifiedPeptides, peptideFilterRegexp);
	}

	public static CensusChroParser getCensusChroParser(File fastaFile, File inputFilefolder, String[] fileNames,
			List<Map<QuantCondition, QuantificationLabel>> labelsByConditions, char[] enzymeArray, int missedCleavages,
			boolean semiCleavage, File uniprotReleasesFolder, String uniprotVersion, String decoyRegexp,
			boolean ignoreNotFoundPeptidesInDB, boolean distinguishModifiedPeptides, String peptideFilterRegexp)
			throws FileNotFoundException {
		List<RemoteSSHFileReference> xmlFiles = new ArrayList<RemoteSSHFileReference>();
		// Set parser (6 files) to peptides
		for (String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.add(new RemoteSSHFileReference(inputXmlFile));
		}
		String fileNamesKey = getFileNamesKey(fileNames);
		if (quantParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return (CensusChroParser) quantParsersByFileNamesKey.get(fileNamesKey);
		}

		CensusChroParser parser = new CensusChroParser(xmlFiles, labelsByConditions, QuantificationLabel.LIGHT,
				QuantificationLabel.HEAVY);
		try {
			parser.setRetrieveFastaIsoforms(false);
			parser.addIonExclusion(IonSerieType.B, 1);
			parser.addIonExclusion(IonSerieType.Y, 1);
			parser.setDecoyPattern(decoyRegexp);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
			final DBIndexInterface fastaDBIndex = getFastaDBIndex(fastaFile, enzymeArray, missedCleavages, semiCleavage,
					peptideFilterRegexp);
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
			String uniprotVersion, String decoyRegexp, boolean ignoreNotFoundPeptidesInDB, String peptideFilterRegexp)
			throws FileNotFoundException {
		// Set parser (6 files) to peptides
		Map<String, RemoteSSHFileReference> xmlFiles = new HashMap<String, RemoteSSHFileReference>();
		for (String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.put(fileName, new RemoteSSHFileReference(inputXmlFile));
		}
		String fileNamesKey = getFileNamesKey(fileNames);
		if (dtaSelectParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return dtaSelectParsersByFileNamesKey.get(fileNamesKey);
		}

		DTASelectParser parser = new DTASelectParser(xmlFiles);
		try {

			parser.setDecoyPattern(decoyRegexp);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
			final DBIndexInterface fastaDBIndex = getFastaDBIndex(fastaFile, enzymeArray, missedCleavages, semiCleavage,
					peptideFilterRegexp);
			parser.setDbIndex(fastaDBIndex);
			parser.enableProteinMergingBySecondaryAccessions(
					getUniprotProteinLocalRetrieverByFolder(uniprotReleasesFolder), uniprotVersion);

			return parser;
		} finally {
			addDTASelectParserToStaticMap(fileNamesKey, parser);
		}
	}

	private static String getFileNamesKey(String[] xmlFiles) {
		StringBuilder sb = new StringBuilder();
		for (String key : xmlFiles) {
			sb.append(key);
		}
		return sb.toString();
	}

	private static CensusOutParser getCensusOutParser(File fastaFile, File inputFilefolder, String[] fileNames,
			List<Map<QuantCondition, QuantificationLabel>> labelsByConditions, char[] enzymeArray, int missedCleavages,
			boolean semiCleavage, File uniprotReleasesFolder, String uniprotVersion, String decoyRegexp,
			boolean ignoreNotFoundPeptidesInDB, boolean onlyOneSpectrumPerChromatographicPeakAndPerSaltStep,
			boolean skipSingletons, boolean distinguishModifiedPeptides, String peptideFilterRegexp)
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
		String fileNamesKey = getFileNamesKey(fileNames);
		if (quantParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return (CensusOutParser) quantParsersByFileNamesKey.get(fileNamesKey);
		}
		CensusOutParser parser = new CensusOutParser(xmlFiles, labelsByConditions, QuantificationLabel.LIGHT,
				QuantificationLabel.HEAVY);
		try {
			parser.setRetrieveFastaIsoforms(false);
			parser.setDecoyPattern(decoyRegexp);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
			parser.setOnlyOneSpectrumPerChromatographicPeakAndPerSaltStep(
					onlyOneSpectrumPerChromatographicPeakAndPerSaltStep);
			parser.setSkipSingletons(skipSingletons);
			final DBIndexInterface fastaDBIndex = getFastaDBIndex(fastaFile, enzymeArray, missedCleavages, semiCleavage,
					peptideFilterRegexp);
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
			char[] enzymeArray, int missedCleavages, boolean semiCleavage, File uniprotReleasesFolder,
			String uniprotVersion, String decoyRegexp, boolean ignoreNotFoundPeptidesInDB,
			boolean distinguishModifiedPeptides, String peptideFilterRegexp) throws FileNotFoundException {
		// Set parser (6 files) to peptides
		List<RemoteSSHFileReference> xmlFiles = new ArrayList<RemoteSSHFileReference>();

		for (String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.add(new RemoteSSHFileReference(inputXmlFile));
		}
		String fileNamesKey = getFileNamesKey(fileNames);
		if (quantParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return (SeparatedValuesParser) quantParsersByFileNamesKey.get(fileNamesKey);
		}
		SeparatedValuesParser parser = new SeparatedValuesParser(xmlFiles, separator, labelsByConditions,
				QuantificationLabel.LIGHT, QuantificationLabel.HEAVY);
		try {
			parser.setRetrieveFastaIsoforms(false);
			parser.setDecoyPattern(decoyRegexp);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
			final DBIndexInterface fastaDBIndex = getFastaDBIndex(fastaFile, enzymeArray, missedCleavages, semiCleavage,
					peptideFilterRegexp);
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

	private static void addDTASelectParserToStaticMap(String fileNamesKey, DTASelectParser parser) {
		dtaSelectParsersByFileNamesKey.put(fileNamesKey, parser);
		log.info(dtaSelectParsersByFileNamesKey.size() + " parsers stored.");
	}

	private static DBIndexInterface getFastaDBIndex(File fastaFile, char[] enzymeArray, int missedCleavages,
			boolean semicleavage, String peptideFilterRegexp) {
		if (fastaFile != null) {

			DBIndexSearchParams defaultDBIndexParams = DBIndexInterface.getDefaultDBIndexParams(fastaFile);
			String fastaIndexKey = IndexUtil.createFullIndexFileName(defaultDBIndexParams);
			if (indexByFastaIndexKey.containsKey(fastaIndexKey)) {
				return indexByFastaIndexKey.get(fastaIndexKey);
			}
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeArr(enzymeArray, missedCleavages, semicleavage);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeOffset(0);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeNocutResidues("");
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setH2OPlusProtonAdded(true);
			if (peptideFilterRegexp != null) {
				((DBIndexSearchParamsImpl) defaultDBIndexParams)
						.setPeptideFilter(new PeptideFilterBySequence(peptideFilterRegexp));
			}
			DBIndexInterface dbIndex = new DBIndexInterface(defaultDBIndexParams);
			indexByFastaIndexKey.put(fastaIndexKey, dbIndex);
			return dbIndex;
		}
		return null;
	}

	private static CensusChroParser getCensusChroParserUsingMongoDBIndex(String mongoDBURI, String mongoMassDBName,
			String mongoSeqDBName, String mongoProtDBName, File inputFilefolder, String[] fileNames,
			List<Map<QuantCondition, QuantificationLabel>> labelsByConditions, File uniprotReleasesFolder,
			String uniprotVersion, String decoyRegexp, boolean ignoreNotFoundPeptidesInDB,
			boolean distinguishModifiedPeptides, String peptideFilterRegexp) throws FileNotFoundException {
		// Set parser (6 files) to peptides
		List<RemoteSSHFileReference> xmlFiles = new ArrayList<RemoteSSHFileReference>();

		for (String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.add(new RemoteSSHFileReference(inputXmlFile));
		}
		String fileNamesKey = getFileNamesKey(fileNames);
		if (quantParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return (CensusChroParser) quantParsersByFileNamesKey.get(fileNamesKey);
		}
		CensusChroParser parser = new CensusChroParser(xmlFiles, labelsByConditions, QuantificationLabel.LIGHT,
				QuantificationLabel.HEAVY);
		try {
			parser.setRetrieveFastaIsoforms(false);
			parser.addIonExclusion(IonSerieType.B, 1);
			parser.addIonExclusion(IonSerieType.Y, 1);
			parser.setDecoyPattern(decoyRegexp);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
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
			boolean onlyOneSpectrumPerChromatographicPeakAndPerSaltStep, boolean skipSingletons,
			boolean distinguishModifiedPeptides, String peptideFilterRegexp) throws FileNotFoundException {
		// Set parser (6 files) to peptides
		List<RemoteSSHFileReference> xmlFiles = new ArrayList<RemoteSSHFileReference>();

		for (String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.add(new RemoteSSHFileReference(inputXmlFile));
		}
		String fileNamesKey = getFileNamesKey(fileNames);
		if (quantParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return (CensusOutParser) quantParsersByFileNamesKey.get(fileNamesKey);
		}
		CensusOutParser parser = new CensusOutParser(xmlFiles, labelsByConditions, QuantificationLabel.LIGHT,
				QuantificationLabel.HEAVY);
		try {
			parser.setRetrieveFastaIsoforms(false);
			parser.setDecoyPattern(decoyRegexp);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
			parser.setOnlyOneSpectrumPerChromatographicPeakAndPerSaltStep(
					onlyOneSpectrumPerChromatographicPeakAndPerSaltStep);
			parser.setSkipSingletons(skipSingletons);
			DBIndexInterface dbIndex = getMongoDBIndex(mongoDBURI, mongoMassDBName, mongoSeqDBName, mongoProtDBName,
					peptideFilterRegexp);
			parser.setDbIndex(dbIndex);
			parser.enableProteinMergingBySecondaryAccessions(
					getUniprotProteinLocalRetrieverByFolder(uniprotReleasesFolder), uniprotVersion);
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
		Map<String, RemoteSSHFileReference> xmlFiles = new HashMap<String, RemoteSSHFileReference>();

		for (String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.put(fileName, new RemoteSSHFileReference(inputXmlFile));
		}
		String fileNamesKey = getFileNamesKey(fileNames);
		if (dtaSelectParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return dtaSelectParsersByFileNamesKey.get(fileNamesKey);
		}
		DTASelectParser parser = new DTASelectParser(xmlFiles);
		try {
			parser.setDecoyPattern(decoyRegexp);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
			DBIndexInterface dbIndex = getMongoDBIndex(mongoDBURI, mongoMassDBName, mongoSeqDBName, mongoProtDBName,
					peptideFilterRegexp);
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
			File uniprotReleasesFolder, String uniprotVersion, String decoyRegexp, boolean ignoreNotFoundPeptidesInDB,
			boolean distinguishModifiedPeptides, String peptideFilterRegexp) throws FileNotFoundException {
		// Set parser (6 files) to peptides
		List<RemoteSSHFileReference> xmlFiles = new ArrayList<RemoteSSHFileReference>();

		for (String fileName : fileNames) {
			final File inputXmlFile = new File(inputFilefolder.getAbsolutePath() + File.separator + fileName);
			if (!inputXmlFile.exists()) {
				throw new FileNotFoundException(inputXmlFile.getAbsolutePath() + " doesn't exist");
			}
			xmlFiles.add(new RemoteSSHFileReference(inputXmlFile));
		}
		String fileNamesKey = getFileNamesKey(fileNames);
		if (quantParsersByFileNamesKey.containsKey(fileNamesKey)) {
			return (SeparatedValuesParser) quantParsersByFileNamesKey.get(fileNamesKey);
		}
		SeparatedValuesParser parser = new SeparatedValuesParser(xmlFiles, separator, labelsByConditions,
				QuantificationLabel.LIGHT, QuantificationLabel.HEAVY);
		try {
			parser.setRetrieveFastaIsoforms(false);
			parser.setDecoyPattern(decoyRegexp);
			parser.setIgnoreNotFoundPeptidesInDB(ignoreNotFoundPeptidesInDB);
			DBIndexInterface dbIndex = getMongoDBIndex(mongoDBURI, mongoMassDBName, mongoSeqDBName, mongoProtDBName,
					peptideFilterRegexp);
			parser.setDbIndex(dbIndex);
			parser.enableProteinMergingBySecondaryAccessions(
					getUniprotProteinLocalRetrieverByFolder(uniprotReleasesFolder), uniprotVersion);

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
			DBIndexSearchParamsImpl params = new DBIndexSearchParamsImpl(mongoDBURI, mongoMassDBName, mongoSeqDBName,
					mongoProtDBName, peptideFilter);
			DBIndexInterface dbIndex = new DBIndexInterface(params);
			return dbIndex;
		}
		return null;
	}

	public static boolean shareAtLeastOnePeptideNode(PCQProteinNode proteinNode1, PCQProteinNode proteinNode2,
			boolean skipDiscardedPeptideNodes) {

		Set<PCQPeptideNode> peptideNodes1 = proteinNode1.getPeptideNodes();
		for (PCQPeptideNode peptideNode : peptideNodes1) {
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
		Set<PCQPeptideNode> peptideNodes1 = proteinNode1.getPeptideNodes();

		for (PCQPeptideNode peptideNode1 : peptideNodes1) {
			if (skipDiscardedPeptideNodes && peptideNode1.isDiscarded()) {
				continue;
			}
			QuantifiedPeptideInterface peptide1 = peptideNode1.getItemsInNode().iterator().next();
			if (peptideAlignments != null) {
				Set<AlignedPeptides> alignments = peptideAlignments.getAlignmentsForPeptide(peptide1);

				for (AlignedPeptides alignment : alignments) {
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
		Set<QuantifiedPeptideInterface> peptides = proteinNode.getQuantifiedPeptides();
		for (QuantifiedPeptideInterface peptide : peptides) {
			if (peptide.getSequence().equals(peptideSeq)) {
				return true;
			}
		}
		return false;
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
			QuantCondition cond2, String replicateName, boolean skipDiscarded) {
		if (peptideNodes.isEmpty()) {
			return null;
		}
		// if there is only one peptide, get the consensus ratio from it
		if (peptideNodes.size() == 1) {
			final PCQPeptideNode peptideNode = peptideNodes.iterator().next();
			if (skipDiscarded && !peptideNode.isDiscarded()) {
				return getRepresentativeRatioForPeptideNode(peptideNode, cond1, cond2, true, replicateName);
			}
			return null;
		}

		// in case of isobaric ratios, calculate a consensus ion count ratio
		if (peptideNodes.iterator().next().getQuantifiedPeptides().iterator()
				.next() instanceof IsobaricQuantifiedPeptide) {
			Set<IsobaricQuantifiedPeptide> isobaricQuantPeptides = new HashSet<IsobaricQuantifiedPeptide>();
			for (PCQPeptideNode peptideNode : peptideNodes) {
				if (skipDiscarded && peptideNode.isDiscarded()) {
					continue;
				}
				for (QuantifiedPeptideInterface peptide : peptideNode.getQuantifiedPeptides()) {
					isobaricQuantPeptides.add((IsobaricQuantifiedPeptide) peptide);
				}

			}
			return QuantUtils.getNormalizedIonCountRatioForPeptides(isobaricQuantPeptides, cond1, cond2, replicateName);
		} else {
			// otherwise calculate an average over the non infinity ratios
			final Set<QuantRatio> consensusRatios = getConsensusRatios(peptideNodes, cond1, cond2, replicateName,
					skipDiscarded);
			Set<QuantRatio> nonInfinityRatios = QuantUtils.getNonInfinityRatios(consensusRatios);
			return QuantUtils.getAverageRatio(nonInfinityRatios, AggregationLevel.PEPTIDE);
		}
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
			QuantCondition cond2, String replicateName, boolean skipDiscarded) {

		Set<QuantRatio> ret = new HashSet<QuantRatio>();
		if (peptideNodes.isEmpty()) {
			return ret;
		}
		for (PCQPeptideNode peptideNode : peptideNodes) {
			if (skipDiscarded && peptideNode.isDiscarded()) {
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
			boolean onlysharedByThisToProteins, boolean skipDiscarded) {
		Set<PCQPeptideNode> peptideNodes1 = proteinNode1.getPeptideNodes();
		Set<PCQPeptideNode> peptideNodes2 = proteinNode2.getPeptideNodes();
		Set<PCQPeptideNode> ret = new HashSet<PCQPeptideNode>();
		for (PCQPeptideNode peptideNode1 : peptideNodes1) {
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

		List<Double> ratioValues = new ArrayList<Double>();

		for (PCQPeptideNode peptideNode : sharedPeptideNodes) {
			if (skipDiscarded && peptideNode.isDiscarded()) {
				continue;
			}
			for (QuantifiedPeptideInterface peptide : peptideNode.getQuantifiedPeptides()) {
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

		List<Double> ratioValues = new ArrayList<Double>();

		for (PCQPeptideNode peptideNode : sharedPeptideNodes) {
			if (skipDiscarded && peptideNode.isDiscarded()) {
				continue;
			}
			for (QuantifiedPSMInterface psm : peptideNode.getQuantifiedPSMs()) {
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

		Set<PCQPeptideNode> peptidesNodes1 = proteinNode1.getPeptideNodes();
		Set<PCQPeptideNode> peptidesNodes2 = proteinNode2.getPeptideNodes();
		Set<PCQPeptideNode> totalPeptideNodes = new HashSet<PCQPeptideNode>();
		totalPeptideNodes.addAll(peptidesNodes1);
		totalPeptideNodes.addAll(peptidesNodes2);
		Map<String, Set<PCQPeptideNode>> map = new HashMap<String, Set<PCQPeptideNode>>();
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
					for (PCQProteinNode proteinNode : proteinNodes) {
						if (skipDiscarded && proteinNode.isDiscarded()) {
							continue;
						}
						String proteinAccKey = proteinNode.getKey();
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
		Set<PCQPeptideNode> peptideNodes1 = getUniquePeptideNodes(proteinNode1, proteinNode2, uniquePepOnly,
				skipDiscarded);
		List<Double> ratios = new ArrayList<Double>();

		for (PCQPeptideNode peptideNode1 : peptideNodes1) {
			if (skipDiscarded && peptideNode1.isDiscarded()) {
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
		Set<PCQPeptideNode> peptideNodes1 = getUniquePeptideNodes(proteinNode1, proteinNode2, uniquePepOnly,
				skipDiscarded);
		List<Double> ratios = new ArrayList<Double>();

		for (PCQPeptideNode peptideNode1 : peptideNodes1) {
			if (skipDiscarded && peptideNode1.isDiscarded()) {
				continue;
			}
			for (QuantifiedPSMInterface psm : peptideNode1.getQuantifiedPSMs()) {
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
		for (Double ratio : ratios) {
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
				return o1.getFullSequence().compareTo(o2.getFullSequence());
			}
		});
		return ret;
	}

	public static List<PCQPeptideNode> getSortedPeptideNodesBySequence(Collection<PCQPeptideNode> peptides) {
		List<PCQPeptideNode> ret = new ArrayList<PCQPeptideNode>();
		for (PCQPeptideNode pcqPeptideNode : peptides) {
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
		List<PCQProteinNode> list = new ArrayList<PCQProteinNode>();
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
	 * Get a separated value list of protein keys after sorting them
	 * alphabetically
	 *
	 * @param proteins
	 * @return
	 */
	public static String getKeyString(Collection<QuantifiedProteinInterface> proteins) {
		Set<String> set = new HashSet<String>();
		for (QuantifiedProteinInterface protein : proteins) {
			set.add(protein.getKey());
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
	public static Set<PCQPeptideNode> getUniquePeptideNodes(PCQProteinNode proteinNode, boolean skipDiscarded) {
		Set<PCQPeptideNode> ret = new HashSet<PCQPeptideNode>();
		for (PCQPeptideNode peptideNode : proteinNode.getPeptideNodes()) {
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
		Set<PCQProteinNode> set = new HashSet<PCQProteinNode>();
		set.add(protein);
		return getProteinMap(set, skipDiscarded);
	}

	private static Map<String, PCQProteinNode> getProteinMap(Collection<PCQProteinNode> proteinNodes,
			boolean skipDiscarded) {
		Map<String, PCQProteinNode> map = new HashMap<String, PCQProteinNode>();
		for (PCQProteinNode proteinNode : proteinNodes) {
			if (skipDiscarded && proteinNode.isDiscarded()) {
				continue;
			}
			final String accession = proteinNode.getKey();
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
			Set<String> validTaxonomies, boolean skipDiscarded) {

		Set<String> set = new HashSet<String>();
		for (PCQProteinNode proteinNode : proteinNodeSet) {
			if (skipDiscarded && proteinNode.isDiscarded()) {
				continue;
			}
			String rawAcc = proteinNode.getKey();
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
			List<Map<QuantCondition, QuantificationLabel>> labelsByConditionsList) throws FileNotFoundException {
		return getQuantParser(params, labelsByConditionsList, params.getQuantInputFileNamesArray());
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
						params.isIgnoreNotFoundPeptidesInDB(), !params.isIgnorePTMs(), params.getPeptideFilterRegexp());
			} else {
				return PCQUtils.getCensusChroParser(params.getFastaFile(), params.getInputFileFolder(),
						inputFileNamesArray, labelsByConditionsList, params.getEnzymeArray(),
						params.getMissedCleavages(), params.isSemiCleavage(), params.getUniprotReleasesFolder(),
						params.getUniprotVersion(), params.getDecoyRegexp(), params.isIgnoreNotFoundPeptidesInDB(),
						!params.isIgnorePTMs(), params.getPeptideFilterRegexp());
			}
		} else if (params.getInputType() == AnalysisInputType.CENSUS_OUT) {
			if (params.getMongoDBURI() != null) {
				final CensusOutParser parser = PCQUtils.getCensusOutParserUsingMongoDBIndex(params.getMongoDBURI(),
						params.getMongoMassDBName(), params.getMongoSeqDBName(), params.getMongoProtDBName(),
						params.getInputFileFolder(), inputFileNamesArray, labelsByConditionsList,
						params.getUniprotReleasesFolder(), params.getUniprotVersion(), params.getDecoyRegexp(),
						params.isIgnoreNotFoundPeptidesInDB(),
						params.isOnlyOneSpectrumPerChromatographicPeakAndPerSaltStep(), params.isSkipSingletons(),
						!params.isIgnorePTMs(), params.getPeptideFilterRegexp());

				return parser;
			} else {
				final CensusOutParser parser = PCQUtils.getCensusOutParser(params.getFastaFile(),
						params.getInputFileFolder(), inputFileNamesArray, labelsByConditionsList,
						params.getEnzymeArray(), params.getMissedCleavages(), params.isSemiCleavage(),
						params.getUniprotReleasesFolder(), params.getUniprotVersion(), params.getDecoyRegexp(),
						params.isIgnoreNotFoundPeptidesInDB(),
						params.isOnlyOneSpectrumPerChromatographicPeakAndPerSaltStep(), params.isSkipSingletons(),
						!params.isIgnorePTMs(), params.getPeptideFilterRegexp());

				return parser;
			}
		} else if (params.getInputType() == AnalysisInputType.SEPARATED_VALUES) {
			if (params.getMongoDBURI() != null) {
				final SeparatedValuesParser parser = PCQUtils.getSeparatedValuesParserUsingMongoDBIndex(
						params.getMongoDBURI(), params.getMongoMassDBName(), params.getMongoSeqDBName(),
						params.getMongoProtDBName(), params.getInputFileFolder(), inputFileNamesArray,
						params.getSeparator(), labelsByConditionsList, params.getUniprotReleasesFolder(),
						params.getUniprotVersion(), params.getDecoyRegexp(), params.isIgnoreNotFoundPeptidesInDB(),
						!params.isIgnorePTMs(), params.getPeptideFilterRegexp());

				return parser;
			} else {
				final SeparatedValuesParser parser = PCQUtils.getSeparatedValuesParser(params.getFastaFile(),
						params.getInputFileFolder(), inputFileNamesArray, params.getSeparator(), labelsByConditionsList,
						params.getEnzymeArray(), params.getMissedCleavages(), params.isSemiCleavage(),
						params.getUniprotReleasesFolder(), params.getUniprotVersion(), params.getDecoyRegexp(),
						params.isIgnoreNotFoundPeptidesInDB(), !params.isIgnorePTMs(), params.getPeptideFilterRegexp());

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
					params.isIgnoreNotFoundPeptidesInDB(), params.getPeptideFilterRegexp());
		}
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
	public static Double getLog2RatioValue(QuantRatio ratio, QuantCondition cond1, QuantCondition cond2) {
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
		Set<PCQPeptideNode> peptideNodes1 = proteinNode1.getPeptideNodes();
		Set<PCQPeptideNode> ret = new HashSet<PCQPeptideNode>();

		for (PCQPeptideNode peptideNode1 : peptideNodes1) {
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

	public static Set<QuantRatio> getNonInfinityIndividualPSMRatios(PCQPeptideNode peptideNode) {
		if (peptideNode == null) {
			return null;
		}
		Set<QuantRatio> ret = new HashSet<QuantRatio>();
		final Set<QuantifiedPSMInterface> quantifiedPSMs = peptideNode.getQuantifiedPSMs();
		for (QuantifiedPSMInterface quantifiedPSMInterface : quantifiedPSMs) {
			ret.addAll(quantifiedPSMInterface.getNonInfinityRatios());
		}
		return ret;
	}

	public static Set<QuantRatio> getNonInfinityIndividualPeptideRatios(PCQPeptideNode peptideNode) {
		if (peptideNode == null) {
			return null;
		}
		Set<QuantRatio> ret = new HashSet<QuantRatio>();
		final Set<QuantifiedPeptideInterface> quantifiedPeptides = peptideNode.getQuantifiedPeptides();
		for (QuantifiedPeptideInterface quantifiedPeptide : quantifiedPeptides) {
			ret.addAll(quantifiedPeptide.getNonInfinityRatios());
		}
		return ret;
	}

	public static Double averageOfRatiosTakingIntoAccountInfinitiesAndNans(Collection<QuantRatio> ratios,
			QuantCondition cond1, QuantCondition cond2) {
		List<Double> ratioValues = new ArrayList<Double>();
		for (QuantRatio ratio : ratios) {
			if (ratio != null) {
				ratioValues.add(ratio.getLog2Ratio(cond1, cond2));
			}
		}
		if (ratioValues.isEmpty()) {
			return null;
		}
		// check if there are all INFINITIES
		boolean areInfinities = areAll(Double.POSITIVE_INFINITY, ratioValues)
				|| areAll(Double.NEGATIVE_INFINITY, ratioValues);
		if (areInfinities) {
			// return it (we assume is only one sign of the infinities here
			return ratioValues.iterator().next();
		} else {
			// if they are all Nan,return nan
			if (areAll(Double.NaN, ratioValues)) {
				return Double.NaN;
			}
			// return an average of the non infinities
			List<Double> nonInfinityNonNanValues = new ArrayList<Double>();
			for (Double ratioValue : ratioValues) {
				if (!ratioValue.isInfinite() && !ratioValue.isNaN()) {
					nonInfinityNonNanValues.add(ratioValue);
				}
			}
			// report the average
			return Maths.mean(nonInfinityNonNanValues.toArray(new Double[0]));
		}
	}

	private static boolean areAll(double value, Collection<Double> values) {
		if (values.isEmpty()) {
			return false;
		}
		for (Double double1 : values) {
			if (!double1.equals(value)) {
				return false;
			}
		}
		return true;
	}

	public static QuantRatio getRepresentativeRatioForPeptideNode(PCQPeptideNode peptideNode, QuantCondition cond1,
			QuantCondition cond2, boolean skipDiscarded, String replicateName) {
		Set<PCQPeptideNode> set = new HashSet<PCQPeptideNode>();
		set.add(peptideNode);
		return getRepresentativeRatioForPeptideNodes(set, cond1, cond2, skipDiscarded, replicateName);
	}

	/**
	 * Gets the representative log2 ratio value for the classifications of the
	 * protein pairs, which is:<br>
	 * - in case of being isobaric isotopologues (analysisInputType=
	 * {@link AnalysisInputType}=CENSUS_CHRO) and SanXot is enabled, the average
	 * of the peptideNode Ri ratios comming from SanXot<br>
	 * - in case of being isobaric isotopolofues and SanXot is not enabled, the
	 * average of the peptideNode average Rc ratios of the individual peptides
	 * in the node<br>
	 * - in case of other quantification techniques, if SanXot is enabled, it is
	 * the average of the consensus ratios of the peptide nodes comming from
	 * SanXot.<br>
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
			QuantCondition cond1, QuantCondition cond2, boolean skipDiscarded, String replicateName) {
		if (peptideNodes.isEmpty()) {
			return new CensusRatio(Double.NaN, false, cond1, cond2, AggregationLevel.PEPTIDE_NODE, null);
		}
		List<QuantRatio> toAverage = new ArrayList<QuantRatio>();
		// ISOTOPOLOGUES
		if (ProteinClusterQuantParameters.getInstance().getInputType() == AnalysisInputType.CENSUS_CHRO) {

			if (ProteinClusterQuantParameters.getInstance().isPerformRatioIntegration()) {
				// SANXOT ENABLED
				// get the average of the sanxot Ri
				for (PCQPeptideNode peptideNode : peptideNodes) {
					if (skipDiscarded && peptideNode.isDiscarded()) {
						continue;
					}
					final QuantRatio consensusRatio = peptideNode.getConsensusRatio(cond1, cond2, replicateName);
					if (consensusRatio != null) {
						toAverage.add(consensusRatio);
					}
				}
			} else {
				// SANXOT DISABLED
				// get the average of the normalized Rc ratios for each peptide
				// node
				for (PCQPeptideNode peptideNode : peptideNodes) {
					if (skipDiscarded && peptideNode.isDiscarded()) {
						continue;
					}
					final IonCountRatio normalizedIonCountRatioForPeptideNode = getNormalizedIonCountRatioForPeptideNode(
							peptideNode, cond1, cond2, replicateName);
					if (normalizedIonCountRatioForPeptideNode != null) {
						toAverage.add(normalizedIonCountRatioForPeptideNode);
					}
				}
			}

		} else {
			// SILAC or others,
			if (ProteinClusterQuantParameters.getInstance().isPerformRatioIntegration()) {
				// SANXOT ENABLED
				// get the average of the sanxot ratios
				for (PCQPeptideNode peptideNode : peptideNodes) {
					if (skipDiscarded && peptideNode.isDiscarded()) {
						continue;
					}
					final QuantRatio consensusRatio = peptideNode.getConsensusRatio(cond1, cond2, replicateName);
					if (consensusRatio != null) {
						toAverage.add(consensusRatio);
					}
				}
			} else {
				// SANXOT DISABLED
				// get the average of the individual psm ratios
				for (PCQPeptideNode peptideNode : peptideNodes) {
					if (skipDiscarded && peptideNode.isDiscarded()) {
						continue;
					}
					final Set<QuantifiedPSMInterface> quantifiedPSMs = peptideNode.getQuantifiedPSMs();
					for (QuantifiedPSMInterface psm : quantifiedPSMs) {
						if (replicateName != null && !psm.getFileNames().contains(replicateName)) {
							continue;
						}
						if (psm instanceof QuantifiedPSM) {
							QuantRatio validRatio = QuantUtils.getRatioValidForAnalysis((QuantifiedPSM) psm);
							if (validRatio != null) {
								toAverage.add(validRatio);
							}
						} else {
							throw new IllegalArgumentException(
									"In case of SILAC,  it has to be a QuantifiedPSMFromCensusOut");
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
			return new CensusRatio(finalValue, true, cond1, cond2, AggregationLevel.PEPTIDE_NODE, null);
		}
		return CensusRatio.NAN_RATIO;

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
	public static List<Double> getIndividualRepresentativeLog2ValuesForEachPeptide(
			Collection<PCQPeptideNode> peptideNodes, QuantCondition cond1, QuantCondition cond2,
			boolean skipDiscarded) {

		List<Double> toAverage = new ArrayList<Double>();
		// ISOTOPOLOGUES
		if (ProteinClusterQuantParameters.getInstance().getInputType() == AnalysisInputType.CENSUS_CHRO) {

			if (ProteinClusterQuantParameters.getInstance().isPerformRatioIntegration()) {
				// SANXOT ENABLED
				// get the individual Ri of the PSMs
				for (PCQPeptideNode peptideNode : peptideNodes) {
					final Set<QuantifiedPSMInterface> quantifiedPSMs = peptideNode.getQuantifiedPSMs();
					for (QuantifiedPSMInterface psm : quantifiedPSMs) {
						if (psm instanceof IsobaricQuantifiedPSM) {
							final Set<IsoRatio> nonInfinityIsoRatios = ((IsobaricQuantifiedPSM) psm)
									.getNonInfinityIsoRatios();
							for (IsoRatio isoRatio : nonInfinityIsoRatios) {
								toAverage.add(isoRatio.getLog2Ratio(cond1, cond2));
							}
						}
					}
				}
			} else {
				// SANXOT DISABLED
				// get the individual Rc ratios of the individual peptides of
				// the node
				for (PCQPeptideNode peptideNode : peptideNodes) {
					final Set<QuantifiedPeptideInterface> peptides = peptideNode.getQuantifiedPeptides();
					for (QuantifiedPeptideInterface peptide : peptides) {
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
			for (PCQPeptideNode peptideNode : peptideNodes) {
				final Set<QuantifiedPSMInterface> quantifiedPSMs = peptideNode.getQuantifiedPSMs();
				for (QuantifiedPSMInterface psm : quantifiedPSMs) {
					if (psm instanceof QuantifiedPSM) {
						QuantRatio validRatio = QuantUtils.getRatioValidForAnalysis((QuantifiedPSM) psm);
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

	public static IonCountRatio getNormalizedIonCountRatioForPeptideNode(PCQPeptideNode peptideNode,
			QuantCondition cond1, QuantCondition cond2, String replicateName) {
		Set<IsobaricQuantifiedPeptide> isobaricQuantifiedPeptides = new HashSet<IsobaricQuantifiedPeptide>();
		final Set<QuantifiedPeptideInterface> peptides = peptideNode.getQuantifiedPeptides();
		for (QuantifiedPeptideInterface quantifiedPeptideInterface : peptides) {
			if (quantifiedPeptideInterface instanceof IsobaricQuantifiedPeptide) {
				isobaricQuantifiedPeptides.add((IsobaricQuantifiedPeptide) quantifiedPeptideInterface);
			}
		}
		if (isobaricQuantifiedPeptides.isEmpty()) {
			return IonCountRatio.NAN_RATIO;
		}
		return QuantUtils.getNormalizedIonCountRatioForPeptides(isobaricQuantifiedPeptides, cond1, cond2,
				replicateName);
	}

	public static String escapeInfinity(Double ratio) {
		if (ratio == null) {
			return null;
		}
		if (Double.isInfinite(ratio)) {
			return "'" + String.valueOf(ratio);
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
		DecimalFormat df = new DecimalFormat("+#.##");
		List<StringPosition> ptmPositionsInProtein = new ArrayList<StringPosition>();
		for (QuantifiedPeptideInterface peptide : peptides) {
			ptmPositionsInProtein.addAll(QuantUtils.getPTMPositionsInProtein(accession, peptide,
					ProteinClusterQuantParameters.getInstance().getUniprotVersion(),
					ProteinClusterQuantParameters.getInstance().getUniprotReleasesFolder()));
		}
		StringBuilder sb = new StringBuilder();
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
		Set<Integer> positions = new HashSet<Integer>();
		boolean first = true;
		for (StringPosition ptmPositionInProtein : ptmPositionsInProtein) {
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

			} catch (NumberFormatException e) {

			}
			sb.append(ptmPositionInProtein.position).append("(").append(string).append(")");
			first = false;
		}
		return sb.toString();
	}

	public static Set<QuantifiedProteinInterface> getProteinsFromPeptides(
			Collection<QuantifiedPeptideInterface> peptides) {
		Set<QuantifiedProteinInterface> ret = new HashSet<QuantifiedProteinInterface>();

		for (QuantifiedPeptideInterface peptide : peptides) {
			ret.addAll(peptide.getQuantifiedProteins());
		}

		return ret;
	}

	public static String getSpeciesString(Set<String> taxonomies) {
		StringBuilder sb = new StringBuilder();
		List<String> list = new ArrayList<String>();
		list.addAll(taxonomies);
		Collections.sort(list);
		for (String taxonomy : list) {
			if (!"".equals(sb.toString())) {
				sb.append(",");
			}
			sb.append(taxonomy);
		}
		return sb.toString();

	}
}
