package edu.scripps.yates.pcq.quantsite.custom;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.BooleanUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.pcq.PCQBatchRunner;
import edu.scripps.yates.pcq.ProteinClusterQuant;
import edu.scripps.yates.pcq.params.PropertiesReader;
import edu.scripps.yates.pcq.params.ProteinClusterQuantParameters;
import edu.scripps.yates.pcq.quantsite.QuantSiteOutputComparator;
import edu.scripps.yates.pcq.quantsite.groups.GroupComparison;
import edu.scripps.yates.utilities.appversion.AppVersion;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.maths.PValueCorrectionType;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.THashSet;
import gnu.trove.set.hash.TIntHashSet;

public class QuantSiteCustomAnalysis {
	private final static Logger log = Logger.getLogger(QuantSiteCustomAnalysis.class);
	private static String path = "Z:\\share\\Salva\\data\\cbamberg\\Alzheimer\\60_samples_run_021019\\";
	private final static String[] files = { "Batch_persite.properties", "Batch_input_NCI60replicates.txt",
			"groupComparisonInputFile.txt", "CellLineID_to_LineNumber.txt" };
	private static final int MAX_ITERATIONS = 1000;

	public static void main(String[] args) {
		final AppVersion version = ProteinClusterQuant.getVersion();
		System.out.println("Running ProteinClusterQuant (PCQ) version " + version.toString() + " in batch mode");

		try {
			path = args[0];

			final File paramsFile = new File(path + files[0]);
			if (!paramsFile.exists()) {
				throw new FileNotFoundException(paramsFile.getAbsolutePath() + " doesn't exist");
			}
			final File batchFile = new File(path + files[1]);
			if (!batchFile.exists()) {
				throw new FileNotFoundException(batchFile.getAbsolutePath() + " doesn't exist");
			}
			final File groupFile = new File(path + files[2]);
			if (!groupFile.exists()) {
				throw new FileNotFoundException(groupFile.getAbsolutePath() + " doesn't exist");
			}
			final File cellLinesFile = new File(path + files[3]);
			if (!cellLinesFile.exists()) {
				throw new FileNotFoundException(cellLinesFile.getAbsolutePath() + " doesn't exist");
			}
			Boolean performRandomization = BooleanUtils.toBooleanObject(args[1]);
			if (performRandomization == null || !performRandomization) {
				performRandomization = false;
				System.out.println("Not performing randomization");
			} else {
				System.out.println("Performing randomization");
			}
			run(batchFile, paramsFile, groupFile, cellLinesFile, PValueCorrectionType.BY, 0.05, 2, 1,
					performRandomization);
			System.out.println("Everything finished correctly.");
			System.exit(0);
		} catch (final Exception e) {
			e.printStackTrace();
			ProteinClusterQuant.errorInParameters(true);
		}
		System.err.println("Some error occurred.");
		System.exit(-1);
	}

	private static void run(File batchFile, File paramaterFile, File groupFile, File cellLinesFile,
			PValueCorrectionType pvalueCorrectionType, double qValueThreshold, int numberSigmas,
			int minNumberOfDiscoveries, boolean performRandomization) throws IOException {
		final List<String> totalSampleList = getTotalSampleList(cellLinesFile);
		final Map<String, InputParams> inputDataFilesPerCellLine = getInputDataFilesPerCellLine(batchFile);
		final List<GroupComparison> comparisons = getGroupComparisons(groupFile);
		File pcqComparatorFile = PCQBatchRunner.getRunPCQResultingLogFile(paramaterFile);
		Map<String, String> peptideNodeTableFilePathsByID = getPeptideNodeTableFilePathsByID(pcqComparatorFile);
		final File generalOutputFile = new File(
				batchFile.getParent() + File.separator + "60Samples_custom_analysis.tsv");
		final Set<String> comparisonsPerformed = new THashSet<String>();
		Files.readAllLines(generalOutputFile.toPath()).stream()
				.forEach(line -> comparisonsPerformed.add(line.split("\t")[0]));
		for (final GroupComparison comparison : comparisons) {

			final String mutantComparisonID = comparison.getComparisonID() + "_mutant";
			final String comparisonID = comparison.getComparisonID();
			final String outputFileName = comparisonID + " vs " + mutantComparisonID;
			if (comparisonsPerformed.contains(outputFileName)) {
				log.info(outputFileName + " was already done. Skipping it...");
				continue;
			}
			if (!peptideNodeTableFilePathsByID.containsKey(mutantComparisonID)) {
				log.info("Running PCQ for " + mutantComparisonID + " which have " + comparison.getSampleList().size()
						+ " cell lines");
				final List<String> inputQuantFileList = getInputFileListFromSampleList(comparison.getSampleList(),
						inputDataFilesPerCellLine);
				final String inputQuantFilesLine = getInputFileLine(inputQuantFileList, mutantComparisonID);
				pcqComparatorFile = PCQBatchRunner.runPCQ(paramaterFile, inputQuantFilesLine, null, mutantComparisonID,
						false);
			} else {
				log.info("Comparison '" + mutantComparisonID + "' was already performed. Skipping it...");
			}

			if (!peptideNodeTableFilePathsByID.containsKey(comparisonID)) {
				final List<String> inputQuantFileList2 = getInputFileListFromNonSampleList(comparison.getSampleList(),
						totalSampleList, inputDataFilesPerCellLine);
				log.info("Running PCQ for " + comparisonID + " which have " + inputQuantFileList2.size()
						+ " cell lines");
				final String inputQuantFilesLine2 = getInputFileLine(inputQuantFileList2, comparison.getComparisonID());
				pcqComparatorFile = PCQBatchRunner.runPCQ(paramaterFile, inputQuantFilesLine2, null, comparisonID,
						false);
			} else {
				log.info("Comparison '" + comparisonID + "' was already performed. Skipping it...");
			}
			peptideNodeTableFilePathsByID = getPeptideNodeTableFilePathsByID(pcqComparatorFile);

			QuantSiteOutputComparator.sampleNamesByFiles.clear();
			final List<File> filesToCompare = new ArrayList<File>();
			final File fileMutant = new File(peptideNodeTableFilePathsByID.get(mutantComparisonID));
			filesToCompare.add(fileMutant);
			QuantSiteOutputComparator.sampleNamesByFiles.put(fileMutant, mutantComparisonID);
			final File fileRest = new File(peptideNodeTableFilePathsByID.get(comparisonID));
			filesToCompare.add(fileRest);
			QuantSiteOutputComparator.sampleNamesByFiles.put(fileRest, comparisonID);

			final QuantSiteOutputComparator outputComparator = new QuantSiteOutputComparator(filesToCompare, null,
					1000.0, outputFileName, pvalueCorrectionType, qValueThreshold, numberSigmas, minNumberOfDiscoveries,
					false, null, null, false);
			outputComparator.run();
			// now, read output file of the comparator,
			final File qvalueTableFile = new File(
					outputComparator.getOutputFolder() + File.separator + outputFileName + ".tsv");
			final List<String> lines = Files.readAllLines(qvalueTableFile.toPath());
			final String[] header = lines.get(0).split("\t");
			int numDiscoveriesIndex = -1;
			for (int index = 0; index < header.length; index++) {
				if (header[index].contains("# discoveries")) {
					numDiscoveriesIndex = index;
					break;
				}
			}
			if (numDiscoveriesIndex == -1) {
				throw new IllegalArgumentException(
						"'# discoveries' column has not found in output file " + qvalueTableFile.getAbsolutePath());
			}
			// now keep all lines with qvalue<threshold
			final List<String> linesToKeep = new ArrayList<String>();
			for (int i = 1; i < lines.size(); i++) {
				final String line = lines.get(i);
				final int numDiscoveries = Integer.valueOf(line.split("\t")[numDiscoveriesIndex]);
				if (numDiscoveries >= minNumberOfDiscoveries) {
					linesToKeep.add(line);
				}
			}
			final int numSignificantSites = linesToKeep.size();
			log.info(numSignificantSites + " sites where significantly changing in " + outputFileName);

			Double numRandomlySignificant = null;
			if (performRandomization && numSignificantSites > 0) {
				numRandomlySignificant = performRandomizationAnalysis(comparison, inputDataFilesPerCellLine,
						paramaterFile, totalSampleList, pvalueCorrectionType, qValueThreshold, numberSigmas,
						minNumberOfDiscoveries);

			}
			// add lines to outputfile

			final FileWriter fw = new FileWriter(generalOutputFile, true);
			if (!linesToKeep.isEmpty()) {
				for (final String line : linesToKeep) {
					fw.write(outputFileName + "\t" + comparison.getComparisonID() + "\t"
							+ comparison.getSampleList().size() + "\t" + line + "\t");
					if (numRandomlySignificant != null) {
						fw.write(numRandomlySignificant + " ");
					}
					fw.write("\n");

				}
			} else {
				fw.write(outputFileName + "\t" + comparison.getComparisonID() + "\t" + comparison.getSampleList().size()
						+ "\n");
			}
			fw.close();
		}

	}

	private static double performRandomizationAnalysis(GroupComparison comparison,
			Map<String, InputParams> inputDataFilesPerCellLine, File paramaterFile, List<String> totalSampleList,
			PValueCorrectionType pvalueCorrectionType, double qValueThreshold, int numberSigmas,
			int minNumberOfDiscoveries) throws IOException {
		final int numRandomSamples = comparison.getSampleList().size();
		final String mutantComparisonID = "random_" + comparison.getComparisonID() + "_mutant";
		final String comparisonID = "random_" + comparison.getComparisonID();
		final String outputFileName = comparisonID + " vs " + mutantComparisonID;
		int iteration = 0;
		final TIntList numSignificantSites = new TIntArrayList();
		while (iteration < MAX_ITERATIONS) {

			final List<String> inputQuantFileList = getRandomInputFileListFromSampleList(numRandomSamples,
					inputDataFilesPerCellLine);
			// set output folder to RANDOMIZATION subfolder
			final File outputFileFolder = new File(
					ProteinClusterQuantParameters.getInstance().getOutputFileFolder().getAbsolutePath() + File.separator
							+ "RANDOMIZATION");
			final String inputQuantFilesLine = getInputFileLine(inputQuantFileList, mutantComparisonID);
			File pcqComparatorFile = PCQBatchRunner.runPCQ(paramaterFile, inputQuantFilesLine, null, mutantComparisonID,
					false, outputFileFolder);

			final List<String> inputQuantFileList2 = getInputFileListFromNonSampleList(inputQuantFileList,
					totalSampleList, inputDataFilesPerCellLine);

			final String inputQuantFilesLine2 = getInputFileLine(inputQuantFileList2, comparison.getComparisonID());
			pcqComparatorFile = PCQBatchRunner.runPCQ(paramaterFile, inputQuantFilesLine2, null, comparisonID, false,
					outputFileFolder);

			final Map<String, String> peptideNodeTableFilePathsByID = getPeptideNodeTableFilePathsByID(
					pcqComparatorFile);

			QuantSiteOutputComparator.sampleNamesByFiles.clear();
			final List<File> filesToCompare = new ArrayList<File>();
			final File fileMutant = new File(peptideNodeTableFilePathsByID.get(mutantComparisonID));
			filesToCompare.add(fileMutant);
			QuantSiteOutputComparator.sampleNamesByFiles.put(fileMutant, mutantComparisonID);
			final File fileRest = new File(peptideNodeTableFilePathsByID.get(comparisonID));
			filesToCompare.add(fileRest);
			QuantSiteOutputComparator.sampleNamesByFiles.put(fileRest, comparisonID);

			final QuantSiteOutputComparator outputComparator = new QuantSiteOutputComparator(filesToCompare, null,
					1000.0, outputFileName, pvalueCorrectionType, qValueThreshold, numberSigmas, minNumberOfDiscoveries,
					false, null, null, false);
			outputComparator.run();
			// now, read output file of the comparator,
			final File qvalueTableFile = new File(
					outputComparator.getOutputFolder() + File.separator + outputFileName + ".tsv");
			final List<String> lines = Files.readAllLines(qvalueTableFile.toPath());
			final String[] header = lines.get(0).split("\t");
			int numDiscoveriesIndex = -1;
			for (int index = 0; index < header.length; index++) {
				if (header[index].contains("# discoveries")) {
					numDiscoveriesIndex = index;
					break;
				}
			}
			if (numDiscoveriesIndex == -1) {
				throw new IllegalArgumentException(
						"'# discoveries' column has not found in output file " + qvalueTableFile.getAbsolutePath());
			}
			// now keep all lines with qvalue<threshold
			final List<String> linesToKeep = new ArrayList<String>();
			for (int i = 1; i < lines.size(); i++) {
				final String line = lines.get(i);
				final int numDiscoveries = Integer.valueOf(line.split("\t")[numDiscoveriesIndex]);
				if (numDiscoveries >= minNumberOfDiscoveries) {
					linesToKeep.add(line);
				}
			}
			numSignificantSites.add(linesToKeep.size());
			iteration++;
			log.info("Iteration: " + iteration + ", average number of significant sites: "
					+ Maths.mean(numSignificantSites));

			// remove pcq run folders
			if (fileMutant.getParentFile().exists()) {
				FileUtils.deleteDirectory(fileMutant.getParentFile());
			}
			if (fileRest.getParentFile().exists()) {
				FileUtils.deleteDirectory(fileRest.getParentFile());
			}
		}
		final double average = Maths.mean(numSignificantSites);
		return average;

	}

	private static List<String> getRandomInputFileListFromSampleList(int numRandomSamples,
			Map<String, InputParams> inputDataFilesPerCellLine) {
		final Random rn = new Random();
		final List<String> ret = new ArrayList<String>();
		final List<InputParams> inputParams = new ArrayList<InputParams>();
		inputParams.addAll(inputDataFilesPerCellLine.values());

		final TIntSet indexes = new TIntHashSet();
		while (indexes.size() < numRandomSamples) {
			int index = rn.nextInt(inputParams.size());
			while (indexes.contains(index)) {
				index = rn.nextInt(inputParams.size());
			}
			indexes.add(index);
			final InputParams inputParam = inputParams.get(index);
			ret.addAll(inputParam.getInputFiles());
		}

		return ret;
	}

	private static Map<String, String> getPeptideNodeTableFilePathsByID(File pcqComparatorFile) throws IOException {
		final Map<String, String> ret = new THashMap<String, String>();
		Files.readAllLines(pcqComparatorFile.toPath())
				.forEach(line -> ret.put(line.split("\t")[0], line.split("\t")[1]));
		return ret;
	}

	private static List<String> getTotalSampleList(File cellLinesFile) throws IOException {

		final List<String> totalSamples = Files.readAllLines(cellLinesFile.toPath()).stream().skip(1l)
				.map(l -> "Line_" + l.split("\t")[1]).distinct().collect(Collectors.toList());
		log.info(totalSamples.size() + " samples in total");
		return totalSamples;
	}

	private static String getInputFileLine(List<String> inputFileList, String comparisonID) {
		final StringBuilder sb = new StringBuilder();
		sb.append(comparisonID).append("[");
		for (int i = 0; i < inputFileList.size(); i++) {
			if (i > 0) {
				sb.append(",");
			}
			sb.append(inputFileList.get(i));
		}
		sb.append("]");
		return sb.toString();
	}

	private static List<String> getInputFileListFromNonSampleList(List<String> sampleList, List<String> totalSampleList,
			Map<String, InputParams> inputDataFilesPerCellLine) {
		final List<String> ret = new ArrayList<String>();
		for (final String sample : totalSampleList) {
			if (sampleList.contains(sample)) {
				continue;
			}
			final InputParams inputParams = inputDataFilesPerCellLine.get(sample);
			if (inputParams == null) {
				throw new IllegalArgumentException("no mapping for sample " + sample);
			}
			ret.addAll(inputParams.getInputFiles());
		}
		return ret;
	}

	private static List<String> getInputFileListFromSampleList(List<String> sampleList,
			Map<String, InputParams> inputDataFilesPerCellLine) {
		final List<String> ret = new ArrayList<String>();
		for (final String sample : sampleList) {
			final InputParams inputParams = inputDataFilesPerCellLine.get(sample);
			ret.addAll(inputParams.getInputFiles());
		}
		return ret;
	}

	private static List<GroupComparison> getGroupComparisons(File groupComparisonsFile) throws IOException {
		final List<GroupComparison> ret = new ArrayList<GroupComparison>();
		final List<String> comparisonLines = Files.readAllLines(groupComparisonsFile.toPath());
		for (final String comparisonLine : comparisonLines) {

			final String[] split = comparisonLine.split("\t");
			final String comparisonID = split[0];
			final String subsetString = split[1];
//			if (!comparisonID.equals("PIK3CA")) {
//				continue;
//			}
			final String sampleListString = split[2];
			final GroupComparison groupComparison = new GroupComparison(comparisonID, subsetString, sampleListString,
					true, null, Double.NaN, Double.NaN, 0, false);
			ret.add(groupComparison);
		}
		// sort comparisons by the ones having more mutated cell lines first
		Collections.sort(ret, new Comparator<GroupComparison>() {

			@Override
			public int compare(GroupComparison o1, GroupComparison o2) {
				return Integer.compare(o2.getSampleList().size(), o1.getSampleList().size());
			}
		});
		return ret;
	}

	private static Map<String, InputParams> getInputDataFilesPerCellLine(File batchFile) throws IOException {
		final Map<String, InputParams> ret = new THashMap<String, InputParams>();

		final List<String> lines = Files.readAllLines(Paths.get(batchFile.toURI()));
		String inputFilesLine = null;
		String inputIDFilesLine = null;
		String suffixLine = null;
		int numLine = 0;
		for (final String line : lines) {
			numLine++;
			final List<String> sublines = new ArrayList<String>();
			if (line.contains(";")) {
				final String[] split = line.split(";");
				for (final String subline : split) {
					sublines.add(subline);
				}
			} else {
				sublines.add(line);
			}
			for (final String subline : sublines) {

				if (subline.contains("=")) {
					final String[] split = subline.split("=");
					final String property = split[0].trim();
					final String value = split[1].trim();
					if (property.equals("inputFiles")) {
						if (inputFilesLine != null) {
							final InputParams inputParam = new InputParams(suffixLine, inputFilesLine);
							if (numLine == 60) {
								log.info("sdf");
							}
							ret.put(inputParam.getId(), inputParam);

							// reset params
							inputFilesLine = value;
							inputIDFilesLine = null;
							suffixLine = null;
						} else {
							inputFilesLine = value;
						}
					} else if (property.equals("inputIDFiles")) {
						inputIDFilesLine = value;
					} else if (property.equals("outputSuffix")) {
						suffixLine = value;
					} else {
						log.info("line " + numLine + " ignored ('" + subline + "')");
					}
				} else {
					log.info("line " + numLine + " ignored ('" + subline + "')");
				}
			}
		}
		if (inputFilesLine != null) {
			final InputParams inputParam = new InputParams(suffixLine, inputFilesLine);
			if (numLine == 60) {
				log.info("sdf");
			}
			ret.put(inputParam.getId(), inputParam);

		}
		return ret;
	}
}

class InputParams {
	private final String id;
	private final List<String> inputFiles = new ArrayList<String>();
	private final String suffix;

	public InputParams(String suffix, String inputFilesLine) {
		super();
		this.id = getExperimentNameFromLines(inputFilesLine);
		this.suffix = suffix;
		inputFiles.addAll(getInputFileNamesFromLines(inputFilesLine));
	}

	private List<String> getInputFileNamesFromLines(String inputFilesLine) {
		final List<String> ret = new ArrayList<String>();
		final ProteinClusterQuantParameters params = ProteinClusterQuantParameters.getInstance();
		params.clearInputQuantificationFiles();
		PropertiesReader.parseInputFiles(inputFilesLine, params);
		params.getInputQuantificationFileNames().stream().forEach(f -> ret.addAll(f.getRelicateFileNames()));
		return ret;
	}

	private String getExperimentNameFromLines(String inputFilesLine) {
		final ProteinClusterQuantParameters params = ProteinClusterQuantParameters.getInstance();
		params.clearInputQuantificationFiles();
		PropertiesReader.parseInputFiles(inputFilesLine, params);
		return params.getInputQuantificationFileNames().get(0).getExperimentName();

	}

	public String getId() {
		return id;
	}

	public List<String> getInputFiles() {
		return inputFiles;
	}

	public String getSuffix() {
		return suffix;
	}

}
