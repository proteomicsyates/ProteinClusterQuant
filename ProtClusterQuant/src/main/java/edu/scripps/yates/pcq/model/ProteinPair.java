package edu.scripps.yates.pcq.model;

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.pcq.ProteinClusterQuantParameters;
import edu.scripps.yates.pcq.cases.Classification1Case;
import edu.scripps.yates.pcq.cases.Classification2Case;
import edu.scripps.yates.pcq.util.PCQUtils;
import edu.scripps.yates.pcq.util.ProteinPairPValue;
import edu.scripps.yates.utilities.maths.Maths;

public class ProteinPair {
	private static final Logger log = Logger.getLogger(ProteinPair.class);
	private final QuantifiedProteinInterface protein1;
	private final QuantifiedProteinInterface protein2;
	private final ProteinCluster proteinCluster;

	private boolean uniquePeptidesProt1Inconsistent = false;
	private boolean uniquePeptidesProt2Inconsistent = false;
	private boolean sharedPeptidesInconsistent = false;

	private final Map<String, Classification2Case> classification2Cases = new HashMap<String, Classification2Case>();
	private final Map<String, Classification1Case> classification1Cases = new HashMap<String, Classification1Case>();
	private Writer outputKMeans;
	ProteinPairPValue firstCase;
	ProteinPairPValue secondCase;
	private final ProteinClusterQuantParameters params = ProteinClusterQuantParameters.getInstance();

	public ProteinPair(PCQProteinNode proteinNode1, PCQProteinNode proteinNode2, ProteinCluster proteinCluster) {
		protein1 = proteinNode1;
		protein2 = proteinNode2;
		this.proteinCluster = proteinCluster;
	}

	public ProteinPair(QuantifiedProteinInterface protein1, QuantifiedProteinInterface protein2,
			ProteinCluster proteinCluster) {
		this.protein1 = protein1;
		this.protein2 = protein2;
		this.proteinCluster = proteinCluster;
	}

	/**
	 * Returns true if the {@link ProteinPair} is inconsistent according to a
	 * QTest
	 *
	 * @param cond1
	 * @param cond2
	 * @return
	 */
	public boolean isInconsistenceWithQTest(QuantCondition cond1, QuantCondition cond2) {
		// here we will have the methods for detecting inconsistencies

		List<Double> listRatioPep = getCountRatios(cond1, cond2);

		if (qTest(listRatioPep)) {
			return true;
			// or false depending on something
			// if true, flag the cluster: proteinCluster.setFlag()
		}
		return false;
	}

	/**
	 * Gets the ratio values from Protein1 in the {@link ProteinPair}
	 *
	 * @param cond1
	 * @param cond2
	 * @return
	 */
	private List<Double> getCountRatiosFromProtein1(QuantCondition cond1, QuantCondition cond2) {
		return PCQUtils.getConsensusRatioValues(protein1, cond1, cond2);
	}

	/**
	 * Gets the ratio values from Protein2 in the {@link ProteinPair}
	 *
	 * @param cond1
	 * @param cond2
	 * @return
	 */
	private List<Double> getCountRatiosFromProtein2(QuantCondition cond1, QuantCondition cond2) {
		return PCQUtils.getConsensusRatioValues(protein2, cond1, cond2);
	}

	/**
	 * Gets the ratio values from Protein1 and Protein2 in the
	 * {@link ProteinPair}
	 *
	 * @param cond1
	 * @param cond2
	 * @return
	 */
	private List<Double> getCountRatios(QuantCondition cond1, QuantCondition cond2) {
		List<Double> countList = new ArrayList<Double>();
		countList.addAll(getCountRatiosFromProtein1(cond1, cond2));
		countList.addAll(getCountRatiosFromProtein2(cond1, cond2));
		return countList;
	}

	/**
	 * Perform a Dixon's qTest over a {@link Collection} of double numbers
	 *
	 * @param numbers
	 * @return
	 */
	public boolean qTest(Collection<Double> numbers) {
		if (numbers.size() < 2) {
			return false;
		}
		List<Double> numberList = new ArrayList<Double>();
		numberList.addAll(numbers);
		Set<Double> numberSet = new HashSet<Double>();
		numberSet.addAll(numbers);
		if (numberSet.size() < 3) {
			return false;
		}

		// sorts the ratios to get the largest and smallest
		Collections.sort(numberList);

		Double largest = numberList.get(numberList.size() - 1);
		Double smallest = numberList.get(0);

		// finds the nearest number to the largest number
		int substractor = 2;
		Double nearestNumLarge = numberList.get(numbers.size() - substractor);
		// takes into account for duplicates
		while (Double.compare(nearestNumLarge, largest) == 0) {
			substractor++;
			if (numbers.size() <= substractor) {
				return false;
			}
			nearestNumLarge = numberList.get(numbers.size() - substractor);
		}

		// finds the nearest number to the smallest number
		int summator = 1;
		Double nearestNumSmall = numberList.get(summator);
		// takes into account for duplicates
		while (Double.compare(nearestNumSmall, smallest) == 0) {
			summator++;
			if (numbers.size() - 2 == summator) {
				return false;
			}
			nearestNumSmall = numberList.get(summator);
		}

		// qTest calculation
		Double qResult1 = Math.abs(smallest - nearestNumSmall) / Math.abs(largest - smallest);
		Double qResult2 = Math.abs(largest - nearestNumLarge) / Math.abs(largest - smallest);

		// sees if the qtest result falls below or above the threshold
		boolean flag = PCQUtils.qTestThreshold(qResult1, qResult2, numberList);

		// true or false, is it inconsistent of consistent
		return flag;
	}

	public String getAccProt1() {
		return protein1.getAccession();
	}

	public String getAccProt2() {
		return protein2.getAccession();
	}

	// gets the number of peptides associated with the protein pair
	public String getNumDataPoints(QuantCondition cond1, QuantCondition cond2) {
		List<Double> ratioList = getCountRatios(cond1, cond2);
		String countString = String.valueOf(ratioList.size());
		return countString;
	}

	// gets the mean of the ratios from the protein pair
	public String getMean(QuantCondition cond1, QuantCondition cond2) {
		List<Double> ratioList = getCountRatios(cond1, cond2);
		double mean = Maths.mean(ratioList.toArray(new Double[0]));
		String meanString = Double.toString(mean);
		return meanString;
	}

	// gets the standard deviation of the ratios from the protein pair
	public String getStandardDeviation(QuantCondition cond1, QuantCondition cond2) {
		List<Double> ratioList = getCountRatios(cond1, cond2);
		double stddev = Maths.stddev(ratioList.toArray(new Double[0]));
		String stddevString = Double.toString(stddev);
		return stddevString;
	}

	// gets the max ratio from the protein pair
	public String getMaxRatio(QuantCondition cond1, QuantCondition cond2) {
		List<Double> ratioList = getCountRatios(cond1, cond2);
		Double max = PCQUtils.findLargestNumber(ratioList);
		String maxString = Double.toString(max);
		return maxString;
	}

	// gets the min ratio from the protein pair
	public String getMinRatio(QuantCondition cond1, QuantCondition cond2) {
		List<Double> ratioList = getCountRatios(cond1, cond2);
		Double min = PCQUtils.findSmallestNumber(ratioList);
		String minString = Double.toString(min);
		return minString;
	}

	public void compareAverages(QuantCondition cond1, QuantCondition cond2) throws IOException {

		List<Double> threeCombined = new ArrayList<Double>();

		// Protein 1's unique peptides
		// has only ratios or inf
		QuantRatio consensusRatio1 = PCQUtils.getUniqueConsensusPeptideRatio(protein1, protein2, cond1, cond2,
				params.isUniquePepOnly());

		Double log2CountRatio1 = PCQUtils.getRatioValue(consensusRatio1, cond1, cond2);
		threeCombined.add(log2CountRatio1);
		List<Double> pepRatProt1 = PCQUtils.getUniquePepRatioValues(protein1, protein2, cond1, cond2,
				params.isUniquePepOnly());

		// shared peptides
		QuantRatio consensusRatio2 = PCQUtils.getSharedConsensusPeptideRatio(protein1, protein2, cond1, cond2);
		Double log2CountRatioShared = PCQUtils.getRatioValue(consensusRatio2, cond1, cond2);
		threeCombined.add(log2CountRatioShared);
		List<Double> pepRatShared = PCQUtils.getSharedPepRatioValues(protein1, protein2, cond1, cond2);

		// protein 2's peptides
		QuantRatio consensusRatio3 = PCQUtils.getUniqueConsensusPeptideRatio(protein2, protein1, cond1, cond2,
				params.isUniquePepOnly());
		Double log2CountRatio3 = PCQUtils.getRatioValue(consensusRatio3, cond1, cond2);
		threeCombined.add(log2CountRatio3);
		List<Double> pepRatProt2 = PCQUtils.getUniquePepRatioValues(protein2, protein1, cond1, cond2,
				params.isUniquePepOnly());

		final Map<String, Set<QuantifiedPeptideInterface>> sharedPeptidesMap = PCQUtils.getSharedPeptidesMap(protein1,
				protein2, false);
		for (String sharedPeptidesProteinKey : sharedPeptidesMap.keySet()) {
			final Set<QuantifiedPeptideInterface> sharedPeptides = sharedPeptidesMap.get(sharedPeptidesProteinKey);
			final Set<QuantifiedPeptideInterface> uniques1 = PCQUtils.getUniquePeptides(protein1, protein2,
					params.isUniquePepOnly());
			final Set<QuantifiedPeptideInterface> uniques2 = PCQUtils.getUniquePeptides(protein2, protein1,
					params.isUniquePepOnly());
			classifyPairByClassification2(params.getThresholdForSignificance(), uniques1, uniques2, sharedPeptides,
					sharedPeptidesProteinKey, cond1, cond2);

			classifyPairByclassification1(threeCombined, pepRatProt1, pepRatShared, pepRatProt2,
					sharedPeptidesProteinKey, cond1, cond2);
		}

	}

	// public Double getConsensusRatio(List<Ratio> ratios, QuantCondition cond1,
	// QuantCondition cond2) {
	//
	// List<Ratio> INFRatios = Utils.getINFRatios(ratios, cond1, cond2);
	// List<Ratio> nonINFRatios = Utils.getNonINFRatios(ratios, cond1, cond2);
	// if (INFRatios.size() > nonINFRatios.size()) {
	// Double infinities = Utils.areAllINFSame(INFRatios, cond1, cond2);
	// return infinities;
	// } else if (nonINFRatios.size() > 0 && (nonINFRatios.size() >=
	// INFRatios.size())
	// || (INFRatios.isEmpty() && !nonINFRatios.isEmpty())) {
	// double mean = Maths.mean(nonINFRatios.toArray(new Double[0]));
	// return mean;
	// } else {
	// return null;
	// }
	//
	// }
	// CASIMIR

	private void classifyPairByClassification2(double thresholdForSignificance,
			Collection<QuantifiedPeptideInterface> uniques1, Collection<QuantifiedPeptideInterface> uniques2,
			Collection<QuantifiedPeptideInterface> shared, String sharedPeptidesProteinKey, QuantCondition cond1,
			QuantCondition cond2) throws IOException {
		// if (protein1.getAccession().equals("P29845") ||
		// protein2.getAccession().equals("B4LMH0")) {
		// System.out.println(this);
		// System.out.println("[" + protein1.getAccession() + ":" +
		// Utils.getPeptidesSequenceString(uniques1) + "] ["
		// + Utils.getPeptidesSequenceString(shared) + "] [" +
		// protein2.getAccession() + ":"
		// + Utils.getPeptidesSequenceString(uniques2) + "] ");
		// }
		List<Double> threeCombined = new ArrayList<Double>();
		final QuantRatio consensusRatio1 = PCQUtils.getConsensusRatio(uniques1, cond1, cond2, null);
		Double value1 = PCQUtils.getRatioValue(consensusRatio1, cond1, cond2);
		threeCombined.add(value1);
		final QuantRatio consensusRatioShared = PCQUtils.getConsensusRatio(shared, cond1, cond2, null);
		Double valueShared = PCQUtils.getRatioValue(consensusRatioShared, cond1, cond2);

		threeCombined.add(valueShared);
		final QuantRatio consensusRatio2 = PCQUtils.getConsensusRatio(uniques2, cond1, cond2, null);
		Double value2 = PCQUtils.getRatioValue(consensusRatio2, cond1, cond2);
		threeCombined.add(value2);

		// variables
		// final List<String> [non_regulated, middle_up, middle_down, first_up,
		// last_up, both_up, non_classifed, error];
		double firstUnique = 0;
		double secondUnique = 0;
		double foldChange = 0;

		// DO pair-wise comparison of ratios; the two unique and the shared

		// ratios with threshold or 2* standard deviation

		// use value defined in thresholdForSignificance
		foldChange = thresholdForSignificance;

		// System.out.println("This is the fold change cutoff in log2: "
		// + foldChange);

		if (threeCombined.contains(null) || threeCombined.contains(Double.NaN)) {
			// unclassified
			classification2Cases.put(sharedPeptidesProteinKey, Classification2Case.CASE7);
			// System.out.println("unclassified");
		}
		// make sure the unique1 is larger or equal to unique2
		else {

			// for (int e =0 ; e < 3; e++){
			// if (Double.isInfinite(threeCombined.get(e))){
			// if(Double.compare(Double.POSITIVE_INFINITY,threeCombined.get(e))==0){
			// threeCombined.set(e, Utils.findLargestNumber(totalRatios));
			// }
			// }
			// }
			boolean swapped = false;
			if (Double.compare(threeCombined.get(0), threeCombined.get(2)) > 0) {
				swapped = true;
				firstUnique = threeCombined.get(2);
				secondUnique = threeCombined.get(0);

				threeCombined.set(0, firstUnique);
				threeCombined.set(2, secondUnique);
			}
			final Double unique1Ratio = threeCombined.get(0);
			final Double sharedRatio = threeCombined.get(1);
			final Double unique2Ratio = threeCombined.get(2);
			// perform classification

			StringBuilder logString = new StringBuilder();
			// make sure order is correct
			if (Double.compare(unique1Ratio, unique2Ratio) > 0) {
				classification2Cases.put(sharedPeptidesProteinKey, Classification2Case.CASE5);
				log.error("error\t");
				throw new IllegalArgumentException("Error sorting the ratios for classifying them.");
			} else {
				// if sharedRatio is between unique1ratio and unique2ratio
				double distanceBetweenRatio1Ratio2 = unique2Ratio - unique1Ratio;
				if (Double.isInfinite(unique1Ratio) && Double.isInfinite(unique2Ratio)) {
					if (Double.compare(unique1Ratio, unique2Ratio) == 0) {
						distanceBetweenRatio1Ratio2 = 0;
					}
				}
				if ((Double.compare(sharedRatio, unique1Ratio) >= 0)
						&& (Double.compare(sharedRatio, unique2Ratio) <= 0)) {
					// check if difference between unique values is larger than
					// foldChange
					if (Double.compare(distanceBetweenRatio1Ratio2, foldChange) > 0) {

						// check if shared value is larger than first unique
						// plus half foldChange or shared value is smaller than
						// second unique minus half foldChange
						if (((Double.compare(sharedRatio, unique1Ratio + distanceBetweenRatio1Ratio2 / 4) > 0)
								&& (Double.compare(sharedRatio, unique2Ratio - distanceBetweenRatio1Ratio2 / 4) < 0))
								// Casimir change for classification
								|| (Double.isInfinite(unique1Ratio) && Double.isInfinite(unique2Ratio)
										&& !Double.isInfinite(distanceBetweenRatio1Ratio2)
										// ||
										// (Double.isInfinite(distanceBetweenRatio1Ratio2)
										&& !Double.isInfinite(sharedRatio))) {
							classification2Cases.put(sharedPeptidesProteinKey, Classification2Case.CASE4);
							uniquePeptidesProt1Inconsistent = true;
							uniquePeptidesProt2Inconsistent = true;
							logString.append("sign_unique_both\t");
						} else {
							classification2Cases.put(sharedPeptidesProteinKey, Classification2Case.CASE3);
							logString.append("sign_unique_1\t");
							if (Double.compare(sharedRatio, unique1Ratio + distanceBetweenRatio1Ratio2 / 4) > 0) {
								if (swapped) {
									uniquePeptidesProt2Inconsistent = true;
								} else {
									uniquePeptidesProt1Inconsistent = true;
								}
							} else if (Double.compare(sharedRatio,
									unique2Ratio - distanceBetweenRatio1Ratio2 / 4) < 0) {
								if (swapped) {
									uniquePeptidesProt1Inconsistent = true;
								} else {
									uniquePeptidesProt2Inconsistent = true;
								}
							}
						}
					} else {
						classification2Cases.put(sharedPeptidesProteinKey, Classification2Case.CASE6);
						logString.append("no_difference_1\t");
					}
				} else {
					// here the shared ratio is outside of the limits between
					// the two unique ratios.

					// check if shared value is smaller than first unique minus
					// fold change or shared value is larger than second unique
					// plus fold change
					if ((Double.compare(sharedRatio, unique1Ratio - foldChange) < 0)
							|| Double.compare(sharedRatio, unique2Ratio + foldChange) > 0) {
						boolean isSharedPeptideSharedByOtherProteinOutOfThePair = isSharedPeptideSharedByOtherProteinOutOfThePair(
								cond1, cond2);
						if (isSharedPeptideSharedByOtherProteinOutOfThePair) {
							classification2Cases.put(sharedPeptidesProteinKey, Classification2Case.CASE2);
							logString.append("sign_shared_with_shared_peptide\t");

						} else {
							classification2Cases.put(sharedPeptidesProteinKey, Classification2Case.CASE1);
							logString.append("sign_shared_1\t");
						}
						sharedPeptidesInconsistent = true;
					} else {
						// check if difference of unique is larger than fold
						// change
						if (Double.compare(distanceBetweenRatio1Ratio2, foldChange) > 0) {
							classification2Cases.put(sharedPeptidesProteinKey, Classification2Case.CASE3);
							logString.append("sign_unique_2\t");
							if (Double.compare(sharedRatio, unique1Ratio) < 0) {
								if (swapped) {
									uniquePeptidesProt1Inconsistent = true;
								} else {
									uniquePeptidesProt2Inconsistent = true;
								}
							} else {// if sharedratio is greater than
									// unique2ratio
								if (swapped) {
									uniquePeptidesProt2Inconsistent = true;
								} else {
									uniquePeptidesProt1Inconsistent = true;
								}
							}
						} else {
							classification2Cases.put(sharedPeptidesProteinKey, Classification2Case.CASE6);
							logString.append("no_difference_2\t");
						}
					}
				}
			}

			// log.debug(logString.toString());

		}

	}

	private void writeOnOutputKMeans(String string) throws IOException {
		if (ProteinClusterQuantParameters.getInstance().isGenerateMiscellaneousFiles()) {
			if (outputKMeans != null) {
				outputKMeans.append(string);
			}
		}
	}

	/**
	 * Classifies all three ratios into the three different cases
	 *
	 * @param threeCombined
	 * @param pepRatProt2
	 * @param pepRatShared
	 * @param pepRatProt1
	 * @param gAM
	 * @throws IOException
	 */
	private void classifyPairByclassification1(List<Double> threeCombined, List<Double> pepRatProt1,
			List<Double> pepRatShared, List<Double> pepRatProt2, String sharedPeptidesProteinKey, QuantCondition cond1,
			QuantCondition cond2) throws IOException {

		boolean hasRatiosAndINF = PCQUtils.hasRatiosAndINF(threeCombined);
		Classification1Case classification1Case = Classification1Case.CASE4;

		if (threeCombined.contains(null) || threeCombined.contains(Double.NaN)) {
			classification1Case = Classification1Case.CASE5;
		} else if (hasRatiosAndINF) {
			classification1Case = Classification1Case.CASE2;
			writeOnOutputKMeans(getAccProt1() + "\t" + classification1Case.getCaseID() + "\n");
			writeOnOutputKMeans(getAccProt2() + "\t" + classification1Case.getCaseID() + "\n");
			if (params.isPrintKMeans()) {
				writeOnOutputKMeans(classification1Case.getCaseID() + "\t" + getAccProt1() + "\t" + getAccProt2() + "\t"
						+ threeCombined.get(0) + "\t" + threeCombined.get(1) + "\t" + threeCombined.get(2) + "\n");
			}
		} else {
			// at this point we have INF OR RATIOS, not both
			List<Double> INFRatios = PCQUtils.getINFRatiosValues(threeCombined);
			List<Double> nonINFRatios = PCQUtils.getNonINFRatiosValues(threeCombined);

			// if we only have RATIO :
			if (!nonINFRatios.isEmpty()) {
				final int uniquePeptideOutlier = outlierTestUniquePeptides(threeCombined, pepRatProt1, pepRatShared,
						pepRatProt2);
				boolean sharedPeptideOutlier = outlierTest(threeCombined.get(1), pepRatProt1, pepRatProt2);
				if (uniquePeptideOutlier != 0) { // there are some outliers
					if (uniquePeptideOutlier == 1 || uniquePeptideOutlier == 3) {
						uniquePeptidesProt1Inconsistent = true;
					} else if (uniquePeptideOutlier == 2 || uniquePeptideOutlier == 3) {
						uniquePeptidesProt2Inconsistent = true;
					}
					classification1Case = Classification1Case.CASE3;
				} else if (sharedPeptideOutlier) {
					sharedPeptideOutlier = true;
					// shared peptide is significant
					// now we want to know if the shared peptide is significant
					// because an external protein or not.
					// whether one of the shared peptides are related to another
					// protein not being the ones in this {@link ProteinPair}
					boolean isSharedPeptideSharedByOtherProteinOutOfThePair = isSharedPeptideSharedByOtherProteinOutOfThePair(
							cond1, cond2);
					sharedPeptidesInconsistent = true;
					if (isSharedPeptideSharedByOtherProteinOutOfThePair) {
						classification1Case = Classification1Case.CASE6;
					} else {
						classification1Case = Classification1Case.CASE7;
					}
				} else {
					// nothing is significant
					classification1Case = Classification1Case.CASE0;
				}
				writeOnOutputKMeans(getAccProt1() + "\t" + classification1Case.getCaseID() + "\n");
				writeOnOutputKMeans(getAccProt2() + "\t" + classification1Case.getCaseID() + "\n");
				if (params.isPrintKMeans()) {
					writeOnOutputKMeans(classification1Case.getCaseID() + "\t" + getAccProt1() + "\t" + getAccProt2()
							+ "\t" + threeCombined.get(0) + "\t" + threeCombined.get(1) + "\t" + threeCombined.get(2)
							+ "\n");
				}
			}
			// if we only have INF or -INF:
			else if (!INFRatios.isEmpty()) {
				Double infinities = PCQUtils.areAllINFValuesSame(INFRatios);
				if (infinities == null) {
					classification1Case = Classification1Case.CASE1;
					uniquePeptidesProt1Inconsistent = true;
					uniquePeptidesProt2Inconsistent = true;
					sharedPeptidesInconsistent = true;
				} else {
					classification1Case = Classification1Case.CASE4;
				}
				writeOnOutputKMeans(getAccProt1() + "\t" + classification1Case.getCaseID() + "\n");
				writeOnOutputKMeans(getAccProt2() + "\t" + classification1Case.getCaseID() + "\n");
				if (params.isPrintKMeans()) {
					writeOnOutputKMeans(classification1Case.getCaseID() + "\t" + getAccProt1() + "\t" + getAccProt2()
							+ "\t" + threeCombined.get(0) + "\t" + threeCombined.get(1) + "\t" + threeCombined.get(2)
							+ "\n");
				}
			}
		}
		if (ProteinClusterQuantParameters.getInstance().isGenerateMiscellaneousFiles() && outputKMeans != null) {
			outputKMeans.flush();
		}

		classification1Cases.put(sharedPeptidesProteinKey, classification1Case);
	}

	private boolean isSharedPeptideSharedByOtherProteinOutOfThePair(QuantCondition cond1, QuantCondition cond2) {
		final Map<String, Set<QuantifiedPeptideInterface>> sharedPeptides = PCQUtils
				.getSharedPeptidesByProtein(protein1, protein2, false);
		if (sharedPeptides.size() > 2) {
			return true;
		}
		return false;
	}

	/**
	 * This function will apply the outlierTest for the forst and the third
	 * elements in threeCombined, taking appropiatelly the population towards
	 * they are compared with
	 *
	 * @param threeCombined
	 * @param pepRatProt1
	 * @param pepRatShared
	 * @param pepRatProt2
	 * @return 1 if unique peptides of protein1 are outliers, 0 if there is no
	 *         outliers and 2 if unique peptides of protein2 are outliers, or 3
	 *         if both unique peptides are outliers
	 */
	private int outlierTestUniquePeptides(List<Double> threeCombined, List<Double> pepRatProt1,
			List<Double> pepRatShared, List<Double> pepRatProt2) {
		List<Double> populationValues = new ArrayList<Double>();
		double valueToTest;
		int ret = 0;
		// first case
		valueToTest = threeCombined.get(0);
		if (outlierTest(valueToTest, pepRatShared, pepRatProt2)) {
			double pValueFirstCase = Maths.iglewiczHoaglinTest(valueToTest, convertingList(populationValues));
			firstCase = new ProteinPairPValue(pValueFirstCase, this);
			ret = 1;
		}

		// third case
		valueToTest = threeCombined.get(2);
		if (outlierTest(valueToTest, pepRatProt1, pepRatShared)) {
			double pValueSecondCase = Maths.iglewiczHoaglinTest(valueToTest, convertingList(populationValues));
			secondCase = new ProteinPairPValue(pValueSecondCase, this);
			if (ret == 1) {
				ret = 3;
			} else {
				ret = 2;
			}
		}
		return ret;
	}

	/**
	 * A generic function to pass the needed values for an outlier test.<br>
	 * Note that the parameter populationSets with "..." means that you can pass
	 * as many List<Double> lists as you want.
	 *
	 * @param valueToTest
	 * @param populationSets
	 * @return
	 */
	@SafeVarargs
	final private boolean outlierTest(Double valueToTest, List<Double>... populationSets) {

		List<Double> populationValues = new ArrayList<Double>();
		boolean outlier = false;

		populationValues.clear();
		populationValues.add(valueToTest);
		for (List<Double> populationSet : populationSets) {
			for (Double pepRatProt1Value : populationSet) {
				if (!Double.isInfinite(pepRatProt1Value) && !Double.isNaN(pepRatProt1Value)) {
					populationValues.add(pepRatProt1Value);
				}
			}
		}
		if (Maths.iglewiczHoaglinTest(valueToTest, convertingList(populationValues)) >= ProteinClusterQuantParameters
				.getInstance().getIglewiczHoaglinTestThreshold()) {
			outlier = true;
		}

		return outlier;
	}

	private double[] convertingList(List<Double> toBeConverted) {
		double[] ret = new double[toBeConverted.size()];
		for (int i = 0; i < toBeConverted.size(); i++) {
			ret[i] = toBeConverted.get(i);
		}
		return ret;
	}

	/**
	 * @return the classification2Case
	 */
	public Map<String, Classification2Case> getClassification2Cases() {
		return classification2Cases;
	}

	/**
	 * @return the classification1Case
	 */
	public Map<String, Classification1Case> getClassification1Case() {
		return classification1Cases;
	}

	public void setOutputKMeans(Writer output) {
		outputKMeans = output;
	}

	public ProteinPairPValue getFirstCase() {
		return firstCase;
	}

	public ProteinPairPValue getSecondCase() {
		return secondCase;
	}

	/**
	 * @return the proteinCluster
	 */
	public ProteinCluster getProteinCluster() {
		return proteinCluster;
	}

	/**
	 * [PROTEIN1: UniqueTo1Peptides] [SharedPeptides] [PROTEIN2:
	 * UniqueTo2Peptides]
	 */
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		// protein1
		List<QuantifiedPeptideInterface> uniquePeptides1 = PCQUtils
				.getSortedPeptidesBySequence(PCQUtils.getUniquePeptides(protein1, protein2, params.isUniquePepOnly()));
		StringBuilder sb1 = new StringBuilder();
		for (QuantifiedPeptideInterface quantifiedPeptide : uniquePeptides1) {
			if (!"".equals(sb1.toString())) {
				sb1.append(",");
			}
			sb1.append(quantifiedPeptide.getSequence());
		}
		sb.append("[" + protein1.getAccession() + ": " + sb1.toString() + "] ");
		// shared
		// protein1
		String shared = PCQUtils.getPeptidesSequenceString(PCQUtils.getSharedPeptidesMap(protein1, protein2, false));
		sb.append("[SHARED: " + shared + "] ");
		// protein2
		List<QuantifiedPeptideInterface> uniquePeptides2 = PCQUtils
				.getSortedPeptidesBySequence(PCQUtils.getUniquePeptides(protein2, protein1, params.isUniquePepOnly()));
		StringBuilder sb2 = new StringBuilder();
		for (QuantifiedPeptideInterface quantifiedPeptide : uniquePeptides2) {
			if (!"".equals(sb2.toString())) {
				sb2.append(",");
			}
			sb2.append(quantifiedPeptide.getSequence());
		}
		sb.append("[" + protein2.getAccession() + ": " + sb2.toString() + "] ");
		return sb.toString();
	}

	/**
	 * @return the protein1
	 */
	public QuantifiedProteinInterface getProtein1() {
		return protein1;
	}

	/**
	 * @return the protein2
	 */
	public QuantifiedProteinInterface getProtein2() {
		return protein2;
	}

	/**
	 * @return the uniquePeptidesProt1Inconsistent
	 */
	public boolean isUniquePeptidesProt1Inconsistent() {
		return uniquePeptidesProt1Inconsistent;
	}

	/**
	 * @return the uniquePeptidesProt2Inconsistent
	 */
	public boolean isUniquePeptidesProt2Inconsistent() {
		return uniquePeptidesProt2Inconsistent;
	}

	/**
	 * @return the sharedPeptidesInconsistent
	 */
	public boolean isSharedPeptidesInconsistent() {
		return sharedPeptidesInconsistent;
	}

	public static String getSummaryLinesHeader() {
		StringBuilder sb = new StringBuilder();

		sb.append("Unique Pep1 (U1)").append("\t").append("Ratio (U1)").append("\t").append("Protein 1").append("\t")
				.append("Shared Pep (S)").append("\t").append("Ratio (S)").append("\t").append("Protein 2")
				.append("\tUnique Pep2 (U2)").append("\t").append("Ratio (U2)").append("\t")
				.append("ClassificationCase1").append("\t").append("ClassificationCase2");
		return sb.toString();
	}

	/**
	 * Returns a list of strings describing the classification cases of this
	 * protein pair as:<br>
	 * UniquePep1 RatioU1 Protein1 SharedPep RatioS Protein2 UniquePep2 RatioU2
	 * ClassificationCase1 ClassificationCase2
	 *
	 * @return
	 */
	public List<String> getSummaryLines(QuantCondition cond1, QuantCondition cond2) {
		List<String> ret = new ArrayList<String>();
		final Map<String, Set<QuantifiedPeptideInterface>> sharedPeptidesMap = PCQUtils.getSharedPeptidesMap(protein1,
				protein2, false);
		for (String proteinKey : sharedPeptidesMap.keySet()) {
			final Set<QuantifiedPeptideInterface> sharedPeptides = sharedPeptidesMap.get(proteinKey);
			final Classification2Case classification2Case = classification2Cases.get(proteinKey);
			final Classification1Case classification1Case = classification1Cases.get(proteinKey);
			StringBuilder sb = new StringBuilder();
			// UniquePep1
			final Set<QuantifiedPeptideInterface> uniquePeptides1 = PCQUtils.getUniquePeptides(protein1, protein2,
					true);
			sb.append(PCQUtils.getPeptidesSequenceString(uniquePeptides1) + "\t");
			// RatioU1
			Double value1 = Double.NaN;
			final QuantRatio consensusRatio1 = PCQUtils.getConsensusRatio(uniquePeptides1, cond1, cond2, null);
			if (consensusRatio1 != null) {
				Double log2Ratio = consensusRatio1.getLog2Ratio(cond1, cond2);
				if (log2Ratio != null) {
					value1 = log2Ratio;
				}
			}
			sb.append(value1 + "\t");
			// Protein1
			sb.append(protein1.getAccession() + "\t");
			// Shared peptide
			sb.append(PCQUtils.getPeptidesSequenceString(sharedPeptides) + "\t");
			// shared peptide ratio
			Double valueShared = Double.NaN;
			final QuantRatio consensusRatioShared = PCQUtils.getConsensusRatio(sharedPeptides, cond1, cond2, null);
			if (consensusRatioShared != null) {
				Double log2Ratio = consensusRatioShared.getLog2Ratio(cond1, cond2);
				if (log2Ratio != null) {
					valueShared = log2Ratio;
				}
			}
			sb.append(valueShared + "\t");
			// Protein2
			sb.append(protein2.getAccession() + "\t");
			// UniquePep2
			final Set<QuantifiedPeptideInterface> uniquePeptides2 = PCQUtils.getUniquePeptides(protein2, protein1,
					true);
			sb.append(PCQUtils.getPeptidesSequenceString(uniquePeptides2) + "\t");
			// RatioU2
			Double value2 = Double.NaN;
			final QuantRatio consensusRatio2 = PCQUtils.getConsensusRatio(uniquePeptides2, cond1, cond2, null);
			if (consensusRatio2 != null) {
				final Double log2Ratio = consensusRatio2.getLog2Ratio(cond1, cond2);
				if (log2Ratio != null) {
					value2 = log2Ratio;
				}
			}
			sb.append(value2 + "\t");
			// classification case 1
			sb.append(classification1Case.name() + "\t");
			// classification case 2
			sb.append(classification2Case.name());
			ret.add(sb.toString());
		}

		return ret;
	}

}
