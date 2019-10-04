package edu.scripps.yates.pcq.quantsite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotEntryUtil;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotGeneMapping;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class QuantSiteOutputComparatorInputFileGenerator {
	private final File mutationsFile;
	private final File cellLineCodeTableFile;
	private final File outputFile;
	private final File uniprotPath;
	private final UniprotProteinLocalRetriever uplr;

	public static void main(String[] args) {
		if (args.length < 4) {
			throw new IllegalArgumentException(
					"Some argument is missing.\n- mutations file\n- cell line code table file\n- output file\n- uniprot folder");
		}
		final File mutationsFile = new File(args[0]);
		final File cellLineCodeTableFile = new File(args[1]);
		final File outputFile = new File(args[2]);
		final File uniprotPath = new File(args[3]);
		final QuantSiteOutputComparatorInputFileGenerator inputFileGenerator = new QuantSiteOutputComparatorInputFileGenerator(
				mutationsFile, cellLineCodeTableFile, outputFile, uniprotPath);
		try {
			final File file = inputFileGenerator.run();
			System.out.println("File created at: " + file.getAbsolutePath());

		} catch (final Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}

	}

	public QuantSiteOutputComparatorInputFileGenerator(File mutationsFile, File cellLineCodeTableFile, File outputFile,
			File uniprotPath) {
		this.mutationsFile = mutationsFile;
		this.cellLineCodeTableFile = cellLineCodeTableFile;
		this.outputFile = outputFile;
		this.uniprotPath = uniprotPath;
		uplr = new UniprotProteinLocalRetriever(uniprotPath, true);
	}

	private File run() throws IOException {
		final UniprotGeneMapping geneMapping = UniprotGeneMapping.getInstance(uniprotPath, "Human");
		final Map<String, String> cellLineCodeByName = getCellLineCodeByName();
		final List<Mutation> mutations = getMutations();
		preLoadAnnotations(mutations, geneMapping);
		final FileWriter fw = new FileWriter(outputFile);
		// get mutations by protein names
		final Map<String, List<Mutation>> mutationsByGeneName = getMutationsByGeneNames(mutations);
		for (final String geneName : mutationsByGeneName.keySet()) {
			fw.write(geneName + "\t");
			final List<Mutation> mutations2 = mutationsByGeneName.get(geneName);
			final Set<String> accs = geneMapping.mapGeneToUniprotACC(geneName);
			String acc = "-";
			if (!accs.isEmpty()) {
				if (accs.size() > 1) {
					// try if I only have a reviewed (swissprot) one
					final String reviewedACC = getReviewedEntry(accs);
					if (reviewedACC != null) {
						acc = reviewedACC;
					} else {
						final StringBuilder sb = new StringBuilder();
						for (final String acc2 : accs) {
							if (!"".equals(sb.toString())) {
								sb.append(",");
							}
							sb.append(acc2);
						}
						acc = sb.toString();
					}
				} else {
					acc = accs.iterator().next();
				}
			}
			fw.write(acc + "\t");
			final StringBuilder sb = new StringBuilder();
			final Set<String> set = new THashSet<String>();
			for (final Mutation mutation : mutations2) {

				if (!cellLineCodeByName.containsKey(mutation.getCellLine())) {
					throw new IllegalArgumentException(
							"No cell line code for cell line name " + mutation.getCellLine());
				}
				final String cellLineCode = cellLineCodeByName.get(mutation.getCellLine());
				if (!set.contains(cellLineCode)) {
					if (!"".equals(sb.toString())) {
						sb.append(",");
					}
					sb.append("Line_" + cellLineCode);
					set.add(cellLineCode);
				}
			}
			fw.write(sb.toString() + "\n");
		}
		fw.close();
		return outputFile;
	}

	private void preLoadAnnotations(List<Mutation> mutations, UniprotGeneMapping geneMapping) throws IOException {
		final Set<String> totalAccs = new THashSet<String>();
		for (final Mutation mutation : mutations) {

			final Set<String> accs = geneMapping.mapGeneToUniprotACC(mutation.getGeneName());
			totalAccs.addAll(accs);
		}
		uplr.getAnnotatedProteins(null, totalAccs);
	}

	private String getReviewedEntry(Set<String> accs) {

		final Map<String, Entry> annotatedProteins = uplr.getAnnotatedProteins(null, accs);
		String reviewed = null;
		for (final String acc : accs) {
			if (annotatedProteins.containsKey(acc)) {
				final Entry entry = annotatedProteins.get(acc);
				if (UniprotEntryUtil.isSwissProt(entry)) {
					if (reviewed != null) {
						return null;
					} else {
						reviewed = acc;
					}
				}
			}
		}
		return reviewed;
	}

	private Map<String, List<Mutation>> getMutationsByGeneNames(List<Mutation> mutations) {
		final Map<String, List<Mutation>> ret = new THashMap<String, List<Mutation>>();
		for (final Mutation mutation : mutations) {
			if (!ret.containsKey(mutation.getGeneName())) {
				ret.put(mutation.getGeneName(), new ArrayList<Mutation>());
			}
			ret.get(mutation.getGeneName()).add(mutation);
		}
		return ret;
	}

	private List<Mutation> getMutations() throws IOException {
		final List<Mutation> ret = new ArrayList<Mutation>();
		final List<String> lines = Files.readAllLines(mutationsFile.toPath());
		for (final String line : lines) {
			final Mutation mutation = new Mutation(line);
			ret.add(mutation);
		}
		return ret;
	}

	private Map<String, String> getCellLineCodeByName() throws IOException {
		final Map<String, String> ret = new THashMap<String, String>();
		final List<String> lines = Files.readAllLines(cellLineCodeTableFile.toPath());
		lines.stream().skip(1).forEach(line -> ret.put(line.split("\t")[0], line.split("\t")[1]));
		return ret;
	}

}
