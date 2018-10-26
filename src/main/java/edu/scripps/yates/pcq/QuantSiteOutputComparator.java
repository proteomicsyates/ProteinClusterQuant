package edu.scripps.yates.pcq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.pcq.compare.model.QuantifiedSite;
import edu.scripps.yates.pcq.compare.model.QuantifiedSiteSet;
import edu.scripps.yates.utilities.appversion.AppVersion;
import edu.scripps.yates.utilities.maths.PValueCorrectionType;
import edu.scripps.yates.utilities.sequence.PositionInPeptide;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import smile.netlib.NLMatrix;

public class QuantSiteOutputComparator {
	private static Logger log = Logger.getLogger(QuantSiteOutputComparator.class);
	private static Options options;
	private final static String currentFolder = System.getProperty("user.dir");
	private final List<File> inputFiles = new ArrayList<File>();
	private final double rInf;
	private final String outputFileName;
	private final PValueCorrectionType pValueCorrectionMethod = PValueCorrectionType.BY;

	public QuantSiteOutputComparator(File inputFile1, File inputFile2, double rInf, String outputFileName) {
		inputFiles.add(inputFile1);
		inputFiles.add(inputFile2);
		this.rInf = rInf;
		this.outputFileName = outputFileName;
	}

	public static void main(String[] args) {
		final AppVersion version = ProteinClusterQuant.getVersion();
		System.out.println("Running Quant Site comparator version " + version.toString());
		setupCommandLineOptions();
		final CommandLineParser parser = new BasicParser();
		QuantSiteOutputComparator quantSiteComparator = null;
		try {
			final CommandLine cmd = parser.parse(options, args);
			if (cmd.getOptionValue("RInf") == null) {
				throw new Exception("Provide input parameter file with 'RInf' option");
			}
			if (cmd.getOptionValue("f1") == null) {
				throw new Exception("Provide input parameter file with 'f1' option");
			}
			if (cmd.getOptionValue("f2") == null) {
				throw new Exception("Provide input parameter file with 'f2' option");
			}
			if (cmd.getOptionValue("out") == null) {
				throw new Exception("Provide input parameter file with 'out' option");
			}
			File inputFile1 = new File(cmd.getOptionValue("f1"));
			File inputFile2 = new File(cmd.getOptionValue("f2"));
			if (!inputFile1.exists()) {
				inputFile1 = new File(currentFolder + File.separator + cmd.getOptionValue("f1"));
				if (!inputFile1.exists()) {
					throw new Exception("Input file 1 not found");
				}
			}
			if (!inputFile2.exists()) {
				inputFile2 = new File(currentFolder + File.separator + cmd.getOptionValue("f2"));
				if (!inputFile2.exists()) {
					throw new Exception("Input file 2 not found");
				}
			}
			log.info("File 1: " + inputFile1.getAbsolutePath());
			log.info("File 2: " + inputFile2.getAbsolutePath());
			Double rInf = null;
			try {
				rInf = Double.valueOf(cmd.getOptionValue("RInf"));
			} catch (final NumberFormatException e) {
				throw new Exception("Option 'RInf' must be numerical");
			}
			final String outputFileName = cmd.getOptionValue("out");
			quantSiteComparator = new QuantSiteOutputComparator(inputFile1, inputFile2, rInf, outputFileName);

		} catch (final Exception e) {
			e.printStackTrace();
			errorInParameters();
		}
		if (quantSiteComparator != null) {
			try {
				quantSiteComparator.run();
			} catch (final IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
		} else {
			System.exit(-1);
		}
		log.info("Program exited correctly.");
		System.exit(0);

	}

	private void run() throws IOException {
		QuantifiedSiteSet quantSites = null;
		for (final File file : inputFiles) {
			if (quantSites == null) {
				quantSites = readFile(file);
			} else {
				final QuantifiedSiteSet quantSites2 = readFile(file);
				quantSites = mergeQuantifiedSiteSets(quantSites, quantSites2);
			}

		}
		writeTripletsOutput(quantSites);
		// writePairWisePValueMatrixes(quantSites);
	}

	private void writePairWisePValueMatrixes(QuantifiedSiteSet quantSites) throws IOException {

		final Map<String, NLMatrix> matrixMap = new THashMap<String, NLMatrix>();
		for (final QuantifiedSite quantSite : quantSites.getSortedByRatios()) {
			final NLMatrix matrix = new NLMatrix(quantSite.getNumExperiments(), quantSite.getNumExperiments());
			matrixMap.put(quantSite.getNodeKey(), matrix);
			for (int sampleIndex1 = 0; sampleIndex1 < quantSite.getNumExperiments(); sampleIndex1++) {
				for (int sampleIndex2 = sampleIndex1 + 1; sampleIndex2 < quantSite
						.getNumExperiments(); sampleIndex2++) {
					final double pValue = performTTest(quantSite, sampleIndex1, sampleIndex2);
					matrix.set(sampleIndex1, sampleIndex2, pValue);
				}
			}

		}

		log.info("Writting output file with " + matrixMap.size() + "matrixes ...");
		final File outputFile = new File(
				currentFolder + File.separator + FilenameUtils.getBaseName(outputFileName) + "_PVALUES_matrixes.txt");
		final FileWriter fw = new FileWriter(outputFile);
		// header
		fw.write(QuantifiedSite.NODE_KEY + "\t");
		for (int index = 0; index < quantSites.getNumExperiments(); index++) {
			fw.write(QuantifiedSite.LOG2RATIO + "_" + (index + 1) + "\t");
			fw.write(QuantifiedSite.RATIOSCOREVALUE + "_" + (index + 1) + "\t");
			// fw.write(QuantifiedSite.NUMPSMS + "_" + (index + 1) + "\t");
			fw.write(QuantifiedSite.NUMMEASUREMENTS + "_" + (index + 1) + "\t");
		}
		fw.write(QuantifiedSite.SEQUENCE + "\t");
		fw.write(QuantifiedSite.POSITIONSINPEPTIDE + "\t");
		fw.write(QuantifiedSite.PROTEINS + "\t");
		fw.write(QuantifiedSite.GENES + "\t");
		fw.write("p-value" + "\t");
		fw.write("corrected p-value (" + pValueCorrectionMethod + ")" + "\t");
		fw.write("\n");
		for (final QuantifiedSite quantSite : quantSites) {
			final QuantifiedSite quantifiedSite = quantSites.getQuantifiedSitesByKey().get(quantSite.getNodeKey());
			fw.write(quantifiedSite.getNodeKey() + "\t");
			for (int index = 0; index < quantSites.getNumExperiments(); index++) {
				fw.write(replaceInfinite(quantifiedSite.getLog2Ratio(index)) + "\t");
				fw.write(quantifiedSite.getRatioScoreValue(index) + "\t");
				fw.write(quantifiedSite.getNumMeasurements(index) + "\t");
			}

			fw.write(quantifiedSite.getSequence() + "\t");
			fw.write(getSeparatedValueStringFromChars(quantifiedSite.getPositionsInPeptide(), "-") + "\t");
			fw.write(quantifiedSite.getProteins() + "\t");
			fw.write(quantifiedSite.getGenes() + "\t");

			fw.write("\n");
		}
		fw.close();
		log.info("Output file written at: '" + outputFile.getAbsolutePath() + "'");
	}

	private void writeTripletsOutput(QuantifiedSiteSet mergedQuantSites) throws IOException {

		// get it sorted
		final List<QuantifiedSite> sortedQuantifiedSites = mergedQuantSites.getSortedByRatios();

		log.info("Writting output file...");
		final File outputFile = new File(currentFolder + File.separator + outputFileName);
		final FileWriter fw = new FileWriter(outputFile);
		// header
		fw.write(QuantifiedSite.NODE_KEY + "\t");
		for (int index = 0; index < mergedQuantSites.getNumExperiments(); index++) {
			fw.write(QuantifiedSite.LOG2RATIO + "_" + (index + 1) + "\t");
			fw.write(QuantifiedSite.RATIOSCOREVALUE + "_" + (index + 1) + "\t");
			// fw.write(QuantifiedSite.NUMPSMS + "_" + (index + 1) + "\t");
			fw.write(QuantifiedSite.NUMMEASUREMENTS + "_" + (index + 1) + "\t");
		}
		fw.write(QuantifiedSite.SEQUENCE + "\t");
		fw.write(QuantifiedSite.POSITIONSINPEPTIDE + "\t");
		fw.write(QuantifiedSite.PROTEINS + "\t");
		fw.write(QuantifiedSite.GENES + "\t");
		fw.write("p-value" + "\t");
		fw.write("corrected p-value (" + pValueCorrectionMethod + ")" + "\t");
		fw.write("\n");
		for (final QuantifiedSite quantSite : sortedQuantifiedSites) {
			final QuantifiedSite quantifiedSite = mergedQuantSites.getQuantifiedSitesByKey()
					.get(quantSite.getNodeKey());
			fw.write(quantifiedSite.getNodeKey() + "\t");
			for (int index = 0; index < mergedQuantSites.getNumExperiments(); index++) {
				fw.write(replaceInfinite(quantifiedSite.getLog2Ratio(index)) + "\t");
				fw.write(quantifiedSite.getRatioScoreValue(index) + "\t");
				fw.write(quantifiedSite.getNumMeasurements(index) + "\t");
			}

			fw.write(quantifiedSite.getSequence() + "\t");
			fw.write(getSeparatedValueStringFromChars(quantifiedSite.getPositionsInPeptide(), "-") + "\t");
			fw.write(quantifiedSite.getProteins() + "\t");
			fw.write(quantifiedSite.getGenes() + "\t");

			fw.write("\n");
		}
		fw.close();
		log.info("Output file written at: '" + outputFile.getAbsolutePath() + "'");

	}

	private double performTTest(QuantifiedSite quantSite, int sampleIndex1, int sampleIndex2) {
		final double mean1 = quantSite.getLog2Ratio(sampleIndex1);
		final double mean2 = quantSite.getLog2Ratio(sampleIndex2);
		final double stdev1 = quantSite.getRatioScoreValue(sampleIndex1);
		final double var1 = Math.pow(stdev1, 2);
		final double stdev2 = quantSite.getRatioScoreValue(sampleIndex2);
		final double var2 = Math.pow(stdev2, 2);
		final int n1 = quantSite.getNumMeasurements(sampleIndex1);
		final int n2 = quantSite.getNumMeasurements(sampleIndex2);

		final edu.scripps.yates.utilities.maths.TTest pValue = edu.scripps.yates.utilities.maths.TTest.test(mean1, var1,
				n1, mean2, var2, n2, true);
		return pValue.pvalue;

	}

	public static String getSeparatedValueStringFromChars(Collection<PositionInPeptide> collection, String separator) {
		final StringBuilder sb = new StringBuilder();
		for (final Object c : collection) {
			if (!"".equals(sb.toString())) {
				sb.append(separator);
			}
			sb.append(c.toString());
		}
		return sb.toString();
	}

	private String replaceInfinite(Double log2Ratio) {
		if (log2Ratio == null) {
			return "";
		}
		if (Double.isInfinite(log2Ratio)) {
			if (Double.POSITIVE_INFINITY == log2Ratio) {
				return String.valueOf(rInf);
			} else if (Double.NEGATIVE_INFINITY == log2Ratio) {
				return String.valueOf(-rInf);
			}
		}
		return String.valueOf(log2Ratio);
	}

	private QuantifiedSiteSet mergeQuantifiedSiteSets(QuantifiedSiteSet quantSites1, QuantifiedSiteSet quantSites2) {
		log.info("Merging two quantified site sets");
		final QuantifiedSiteSet ret = new QuantifiedSiteSet();
		final Set<String> keySet1 = quantSites1.getQuantifiedSitesByKey().keySet();
		for (final String keySet : keySet1) {
			if (quantSites2.getQuantifiedSitesByKey().contains(keySet)) {
				final QuantifiedSite quantSite1 = quantSites1.getQuantifiedSitesByKey().get(keySet);
				final QuantifiedSite quantSite2 = quantSites2.getQuantifiedSitesByKey().get(keySet);
				if (Double.isNaN(quantSite1.getLog2Ratio(0))) {
					if (Double.isNaN(quantSite2.getLog2Ratio(0))) {
						continue;
					}
				}
				// use the first for the result, adding the values of the second
				quantSite1.addLog2Ratio(quantSite2.getLog2Ratio(0));
				quantSite1.addNumMeasurements(quantSite2.getNumMeasurements(0));
				quantSite1.addNumPeptides(quantSite2.getNumPeptides(0));
				quantSite1.addNumPSMs(quantSite2.getNumPSMs(0));
				quantSite1.addRatioScoreValue(quantSite2.getRatioScoreValue(0));
				quantSite1.getPositionsInPeptide().addAll(quantSite2.getPositionsInPeptide());
				ret.add(quantSite1);
			}
		}
		log.info("Result of merging is " + ret.getQuantifiedSitesByKey().size() + " quantified sites from "
				+ quantSites1.getQuantifiedSitesByKey().size() + " + " + quantSites2.getQuantifiedSitesByKey().size());
		log.info("Now having " + ret.getNumExperiments() + " experiments");
		return ret;
	}

	private QuantifiedSiteSet readFile(File inputFile) throws IOException {
		int numLine = 1;
		String line = null;
		try {
			log.info("Reading input file: '" + inputFile.getAbsolutePath() + "'");
			final BufferedReader br = new BufferedReader(new FileReader(inputFile));
			final QuantifiedSiteSet ret = new QuantifiedSiteSet();

			final TObjectIntHashMap<String> indexesByHeaders = new TObjectIntHashMap<String>();
			while ((line = br.readLine()) != null) {
				try {
					final String[] split = line.split("\t");
					if (numLine == 1) {
						for (int i = 0; i < split.length; i++) {
							indexesByHeaders.put(split[i], i);
						}
					} else {
						final QuantifiedSite quantSite = new QuantifiedSite(split, indexesByHeaders);
						ret.add(quantSite);
					}
				} finally {
					numLine++;
				}
			}
			br.close();
			log.info(ret.getQuantifiedSitesByKey().size() + " quantified sites read from input file '"
					+ inputFile.getAbsolutePath() + "'");
			return ret;
		} catch (final Exception e) {
			log.error("LINE  " + numLine + ": " + line);
			log.error("Error reading at line " + numLine + " of file '" + inputFile.getAbsolutePath() + "'");
			throw e;
		}
	}

	private static void setupCommandLineOptions() {
		// create Options object
		options = new Options();
		options.addOption("RInf", true, "[MANDATORY] -RInf replaces +/- Infinity with a user defined (+/-)value");
		options.addOption("f1", true,
				"[MANDATORY] Full path to output file for peptideNodeTable_perSite, of a PCQ run to compare");
		options.addOption("f2", true,
				"[MANDATORY] Full path to output file for peptideNodeTable_perSite, of a PCQ run to compare");
		options.addOption("out", true, "[MANDATORY] Output file name that will be created in the current folder");

	}

	private static void errorInParameters() {
		// automatically generate the help statement
		final HelpFormatter formatter = new HelpFormatter();

		formatter.printHelp(
				"`java -jar QuantSiteoutputComparator.jar -RInf [n] -f1 [nnn_peptideNodeTable_perSite.txt_file1] -f2 [mmm_peptideNodeTable_perSite.txt_file2] -out [output_file_name]",
				"\n\n\n", options, "\n\nContact Salvador Martinez-Bartolome at salvador@scripps.edu for more help");

		System.exit(0);
	}

}
