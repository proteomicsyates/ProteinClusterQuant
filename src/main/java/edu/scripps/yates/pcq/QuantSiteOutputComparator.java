package edu.scripps.yates.pcq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Set;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.log4j.Logger;

import edu.scripps.yates.pcq.compare.model.QuantifiedSite;
import edu.scripps.yates.pcq.compare.model.QuantifiedSiteSet;
import edu.scripps.yates.utilities.sequence.PositionInPeptide;
import gnu.trove.map.hash.TObjectIntHashMap;

public class QuantSiteOutputComparator {
	private static Logger log = Logger.getLogger(QuantSiteOutputComparator.class);
	private static Options options;
	private final static String currentFolder = System.getProperty("user.dir");
	private final File inputFile1;
	private final File inputFile2;
	private final double rInf;
	private final String outputFileName;

	public QuantSiteOutputComparator(File inputFile1, File inputFile2, double rInf, String outputFileName) {
		this.inputFile1 = inputFile1;
		this.inputFile2 = inputFile2;
		this.rInf = rInf;
		this.outputFileName = outputFileName;
	}

	public static void main(String[] args) {
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
		final QuantifiedSiteSet quantSites1 = readFile(inputFile1);
		final QuantifiedSiteSet quantSites2 = readFile(inputFile2);
		writeOutput(quantSites1, quantSites2);
	}

	private void writeOutput(QuantifiedSiteSet quantSites1, QuantifiedSiteSet quantSites2) throws IOException {

		// construct a list of the merged quant sites
		final QuantifiedSiteSet mergedQuantSites = mergeQuantifiedSiteSets(quantSites1, quantSites2);

		// get it sorted
		final List<QuantifiedSite> sortedQuantifiedSites = mergedQuantSites.getSortedByRatios();

		log.info("Writting output file...");
		final File outputFile = new File(currentFolder + File.separator + outputFileName);
		final FileWriter fw = new FileWriter(outputFile);
		// header
		fw.write(QuantifiedSite.NODE_KEY + "\t");
		fw.write(QuantifiedSite.LOG2RATIO + "\t");
		fw.write(QuantifiedSite.RATIOSCOREVALUE + "\t");
		fw.write(QuantifiedSite.NUMPSMS + "\t");
		fw.write(QuantifiedSite.NUMMEASUREMENTS + "\t");

		fw.write(QuantifiedSite.LOG2RATIO + "_2" + "\t");
		fw.write(QuantifiedSite.RATIOSCOREVALUE + "_2" + "\t");
		fw.write(QuantifiedSite.NUMPSMS + "_2" + "\t");
		fw.write(QuantifiedSite.NUMMEASUREMENTS + "_2" + "\t");
		fw.write(QuantifiedSite.SEQUENCE + "\t");
		fw.write(QuantifiedSite.POSITIONSINPEPTIDE + "\t");
		fw.write(QuantifiedSite.PROTEINS + "\t");
		fw.write(QuantifiedSite.GENES + "\t");
		fw.write("\n");
		for (final QuantifiedSite quantifiedSite : sortedQuantifiedSites) {
			fw.write(quantifiedSite.getNodeKey() + "\t");
			fw.write(replaceInfinite(quantifiedSite.getLog2Ratio()) + "\t");
			fw.write(quantifiedSite.getRatioScoreValue() + "\t");
			fw.write(quantifiedSite.getNumPSMs() + "\t");
			fw.write(quantifiedSite.getNumMeasurements() + "\t");

			fw.write(replaceInfinite(quantifiedSite.getLog2Ratio2()) + "\t");
			fw.write(quantifiedSite.getRatioScoreValue2() + "\t");
			fw.write(quantifiedSite.getNumPSMs2() + "\t");
			fw.write(quantifiedSite.getNumMeasurements2() + "\t");
			fw.write(quantifiedSite.getSequence() + "\t");

			fw.write(getSeparatedValueStringFromChars(quantifiedSite.getPositionsInPeptide(), "-") + "\t");
			fw.write(quantifiedSite.getProteins() + "\t");
			fw.write(quantifiedSite.getGenes() + "\t");
			fw.write("\n");
		}
		fw.close();
		log.info("Output file written at: '" + outputFile.getAbsolutePath() + "'");

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
				if (Double.isNaN(quantSite1.getLog2Ratio())) {
					if (Double.isNaN(quantSite2.getLog2Ratio())) {
						continue;
					}
				}
				// use the first for the result, adding the values of the second
				quantSite1.setLog2Ratio2(quantSite2.getLog2Ratio());
				quantSite1.setNumMeasurements2(quantSite2.getNumMeasurements());
				quantSite1.setNumPeptides2(quantSite2.getNumPeptides());
				quantSite1.setNumPSMs2(quantSite2.getNumPSMs());
				quantSite1.setRatioScoreValue2(quantSite2.getRatioScoreValue());
				quantSite1.getPositionsInPeptide().addAll(quantSite2.getPositionsInPeptide());
				ret.add(quantSite1);
			}
		}
		log.info("Result of merging is " + ret.getQuantifiedSitesByKey().size() + " quantified sites from "
				+ quantSites1.getQuantifiedSitesByKey().size() + " + " + quantSites2.getQuantifiedSitesByKey().size());
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
