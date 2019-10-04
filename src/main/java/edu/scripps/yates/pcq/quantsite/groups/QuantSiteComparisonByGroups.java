package edu.scripps.yates.pcq.quantsite.groups;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Random;

import org.apache.log4j.Logger;

import edu.scripps.yates.pcq.compare.model.QuantifiedSite;
import edu.scripps.yates.pcq.compare.model.QuantifiedSiteSet;
import edu.scripps.yates.pcq.util.PCQUtils;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.maths.PValueCorrectionType;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

public class QuantSiteComparisonByGroups {

	private static int MAX_ITERATIONS = 1000;
	private static final Logger log = Logger.getLogger(QuantSiteComparisonByGroups.class);

	public static void performGroupComparisons(QuantifiedSiteSet quantSites, File groupComparisonsFile,
			boolean subsetSites, File outputFile, double qvalueThreshold, PValueCorrectionType pValueCorrectionMethod,
			double distributionAvg, double distributionSigma, int numberSigmas, boolean useMayorityRule)
			throws IOException {
		final List<String> totalSampleNames = quantSites.getUniqueSampleNames();
		final List<String> comparisonLines = Files.readAllLines(groupComparisonsFile.toPath());
		final FileWriter fw = new FileWriter(outputFile);
		fw.write(
				"comparison ID\tprotein(s)-protein site(s)\tSample(s)\tSite\tAvg ratio 1\tStdev ratio 1\tAvg ratio 2\tStdev ratio 2\tp-value\t# random significant sites\n");
		for (final String comparisonLine : comparisonLines) {

			final String[] split = comparisonLine.split("\t");
			final String comparisonID = split[0];
			final String subsetString = split[1];
			if (!comparisonID.equals("PIK3CA")) {
				continue;
			}
			final String sampleListString = split[2];
			int numSignificative = 0;
			final GroupComparison groupComparison = new GroupComparison(comparisonID, subsetString, sampleListString,
					subsetSites, pValueCorrectionMethod, distributionAvg, distributionSigma, numberSigmas,
					useMayorityRule);
			groupComparison.setSites(quantSites);
			final TObjectDoubleMap<QuantifiedSite> analyzedSites = groupComparison.run(qvalueThreshold);
			for (final QuantifiedSite quantSite : analyzedSites.keySet()) {
				final double pvalue = analyzedSites.get(quantSite);
				if (pvalue < qvalueThreshold) {
					numSignificative++;
				}
				fw.write(groupComparison.getComparisonID() + "\t" + groupComparison.getSubsetString() + "\t"
						+ groupComparison.getSampleListString() + "\t" + quantSite.getNodeKey() + "\t");
				fw.write(PCQUtils.escapeInfinity(groupComparison.getGroup1AvgRatio(quantSite)) + "\t"
						+ PCQUtils.escapeInfinity(groupComparison.getGroup1Stdev(quantSite)) + "\t");

				fw.write(PCQUtils.escapeInfinity(groupComparison.getGroup2AvgRatio(quantSite)) + "\t"
						+ PCQUtils.escapeInfinity(groupComparison.getGroup2Stdev(quantSite)) + "\t" + pvalue + "\n");
			}
			if (numSignificative > 0) {
				final TIntArrayList numSignificantsAfterRandomizations = new TIntArrayList();
				final ProgressCounter counter = new ProgressCounter(MAX_ITERATIONS,
						ProgressPrintingType.PERCENTAGE_STEPS, 0);
				counter.setShowRemainingTime(true);
				for (int iteration = 1; iteration < MAX_ITERATIONS; iteration++) {
					counter.increment();
					final String printIfNecessary = counter.printIfNecessary();
					if (!"".equals(printIfNecessary)) {
						log.info(printIfNecessary);
					}
					final GroupComparison randomizedGroupComparison = getRandomizedGroupComparison(groupComparison,
							totalSampleNames, groupComparison.getSampleList().size(), pValueCorrectionMethod,
							distributionAvg, distributionSigma, numberSigmas, useMayorityRule);
					randomizedGroupComparison.setSites(quantSites);
					final TObjectDoubleMap<QuantifiedSite> qValuesByQuantSite = randomizedGroupComparison
							.run(qvalueThreshold);
					int numSignificativeInRandomization = 0;
					for (final double qvalue : qValuesByQuantSite.values()) {
						if (qvalue < qvalueThreshold) {
							numSignificativeInRandomization++;
						}
					}
					numSignificantsAfterRandomizations.add(numSignificativeInRandomization);
				}
				final double avgNumSignificant = Maths.mean(numSignificantsAfterRandomizations);
				fw.write("Num significant\t" + numSignificative + "\tNum significant by random:\t" + avgNumSignificant);
			}
			fw.write("\n");
		}
		fw.close();
	}

	private static GroupComparison getRandomizedGroupComparison(GroupComparison groupComparison,
			List<String> totalSampleNames, int numSamplesToPick, PValueCorrectionType pValueCorrectionMethod,
			double distributionAvg, double distributionSigma, int numberSigmas, boolean useMayorityRule) {
		final int max = totalSampleNames.size() - 1;

		// take numSamplesToPick out of the list of total samples
		final Random rn = new Random();
		final TIntSet indexes = new TIntHashSet();

		while (indexes.size() < numSamplesToPick) {
			final int sampleIndex = rn.nextInt(max + 1);
			indexes.add(sampleIndex);
		}
		final StringBuilder newSampleListString = new StringBuilder();
		for (final int index : indexes.toArray()) {
			if (!"".equals(newSampleListString.toString())) {
				newSampleListString.append(",");
			}
			newSampleListString.append(totalSampleNames.get(index));
		}
		final GroupComparison ret = new GroupComparison(groupComparison.getComparisonID() + "_randomized", null,
				newSampleListString.toString(), false, pValueCorrectionMethod, distributionAvg, distributionSigma,
				numberSigmas, useMayorityRule);

		return ret;
	}

}
