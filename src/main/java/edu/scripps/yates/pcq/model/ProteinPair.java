package edu.scripps.yates.pcq.model;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.pcq.cases.Classification1Case;
import edu.scripps.yates.pcq.cases.Classification2Case;
import edu.scripps.yates.pcq.params.ProteinClusterQuantParameters;
import edu.scripps.yates.pcq.util.PCQUtils;
import edu.scripps.yates.pcq.util.ProteinPairPValue;
import edu.scripps.yates.utilities.maths.Maths;
import gnu.trove.map.hash.THashMap;

public class ProteinPair {
	private static final Logger log = Logger.getLogger(ProteinPair.class);
	private PCQProteinNode proteinNode1;
	private PCQProteinNode proteinNode2;
	private boolean uniquePeptidesProt1Inconsistent = false;
	private boolean uniquePeptidesProt2Inconsistent = false;
	private boolean sharedPeptidesInconsistent = false;
	private final boolean containsDiscardedProteinNode;
	private final Map<String, Classification2Case> classification2Cases = new THashMap<String, Classification2Case>();
	private final Map<String, Classification1Case> classification1Cases = new THashMap<String, Classification1Case>();
	ProteinPairPValue firstCase;
	ProteinPairPValue secondCase;
	private final ProteinClusterQuantParameters params = ProteinClusterQuantParameters.getInstance();
	private final static int MIN_RATIOS_TO_TEST = 10;

	public ProteinPair(PCQProteinNode proteinNode1, PCQProteinNode proteinNode2) {
		// sort consistently the protein pairs in order to show the proteinNode1
		// as the node with a unique peptide node and sort them alphabetically
		boolean protein1HasUniquePeptides = false;
		boolean protein2HasUniquePeptides = false;

		if (!PCQUtils.getUniquePeptideNodes(proteinNode1, proteinNode2, false, true).isEmpty()) {
			protein1HasUniquePeptides = true;
		}
		if (!PCQUtils.getUniquePeptideNodes(proteinNode2, proteinNode1, false, true).isEmpty()) {
			protein2HasUniquePeptides = true;
		}
		if (protein1HasUniquePeptides && !protein2HasUniquePeptides) {
			this.proteinNode1 = proteinNode1;
			this.proteinNode2 = proteinNode2;
		} else if (!protein1HasUniquePeptides && protein2HasUniquePeptides) {
			this.proteinNode1 = proteinNode2;
			this.proteinNode2 = proteinNode1;
		} else if (Boolean.compare(protein1HasUniquePeptides, protein2HasUniquePeptides) == 0) {
			// sort by acc
			final int res = proteinNode1.getKey().compareTo(proteinNode2.getKey());
			if (res >= 0) {
				this.proteinNode1 = proteinNode1;
				this.proteinNode2 = proteinNode2;
			} else if (res < 0) {
				this.proteinNode1 = proteinNode2;
				this.proteinNode2 = proteinNode1;
			}
		}

		this.proteinNode1.setProteinPair(this);
		this.proteinNode2.setProteinPair(this);
		containsDiscardedProteinNode = proteinNode1.isDiscarded() || proteinNode2.isDiscarded();

	}

	public String getAccProt1() {
		return proteinNode1.getKey();
	}

	public String getAccProt2() {
		return proteinNode2.getKey();
	}

	public void proteinPairAnalysis(QuantCondition cond1, QuantCondition cond2) throws IOException {
		final Double[] threeCombined = new Double[3];
		// Protein 1's unique peptides
		// has only ratios or inf
		final Set<PCQPeptideNode> uniqueTo1PeptideNodes = PCQUtils.getUniquePeptideNodes(proteinNode1, proteinNode2,
				params.isUniquePepOnly(), true);
		final Double ratioValueUniqueTo1 = PCQUtils
				.getRepresentativeRatioForPeptideNodes(uniqueTo1PeptideNodes, cond1, cond2, null, true)
				.getLog2Ratio(cond1, cond2);
		threeCombined[0] = ratioValueUniqueTo1;

		final List<Double> pepRatProt1 = PCQUtils
				.getIndividualRepresentativeLog2ValuesForEachPeptideForProteinPairAnalysis(uniqueTo1PeptideNodes, cond1,
						cond2, true);

		// protein 2's peptides
		final Set<PCQPeptideNode> uniqueTo2PeptideNodes = PCQUtils.getUniquePeptideNodes(proteinNode2, proteinNode1,
				params.isUniquePepOnly(), true);
		final Double ratioValueUniqueTo2 = PCQUtils
				.getRepresentativeRatioForPeptideNodes(uniqueTo2PeptideNodes, cond1, cond2, null, true)
				.getLog2Ratio(cond1, cond2);
		threeCombined[2] = ratioValueUniqueTo2;
		final List<Double> pepRatProt2 = PCQUtils
				.getIndividualRepresentativeLog2ValuesForEachPeptideForProteinPairAnalysis(uniqueTo2PeptideNodes, cond1,
						cond2, true);

		final Map<String, Set<PCQPeptideNode>> sharedPeptidesMap = PCQUtils.getSharedPeptideNodesMap(proteinNode1,
				proteinNode2, false, true);
		for (final String sharedPeptidesProteinKey : sharedPeptidesMap.keySet()) {
			final Set<PCQPeptideNode> sharedPeptideNodes = sharedPeptidesMap.get(sharedPeptidesProteinKey);

			final QuantRatio representativeRatioForPeptideNodes = PCQUtils
					.getRepresentativeRatioForPeptideNodes(sharedPeptideNodes, cond1, cond2, null, true);
			final Double ratioValueShared = representativeRatioForPeptideNodes.getLog2Ratio(cond1, cond2);
			threeCombined[1] = ratioValueShared;
			final List<Double> pepRatShared = PCQUtils
					.getIndividualRepresentativeLog2ValuesForEachPeptideForProteinPairAnalysis(sharedPeptideNodes,
							cond1, cond2, true);

			classifyPairByClassification2(params.getThresholdForSignificance(), Arrays.asList(threeCombined),
					sharedPeptidesProteinKey, cond1, cond2, pepRatProt1, pepRatShared, pepRatProt2);

			// classifyPairByclassification1(Arrays.asList(threeCombined),
			// pepRatProt1, pepRatShared, pepRatProt2,
			// sharedPeptidesProteinKey, cond1, cond2);
		}

	}

	/**
	 * Classifies all three ratios into the three different cases based on
	 * infinities and using a user defined foldchange threshold when possible.
	 *
	 * @param thresholdForSignificance
	 * @param uniques1
	 * @param uniques2
	 * @param shared
	 * @param sharedPeptidesProteinKey
	 * @param cond1
	 * @param cond2
	 * @throws IOException
	 */
	private void classifyPairByClassification2(double foldChange, List<Double> threeCombinedOriginal,
			String sharedPeptidesProteinKey, QuantCondition cond1, QuantCondition cond2, List<Double> pepRatProt1,
			List<Double> pepRatShared, List<Double> pepRatProt2) throws IOException {

		double firstUnique = 0;
		double secondUnique = 0;
		// copy in another array, otherwise it may get swapped
		final List<Double> threeCombined = new ArrayList<Double>();
		threeCombined.addAll(threeCombinedOriginal);
		// DO pair-wise comparison of ratios; the two unique and the shared
		final double log2ThresholdFoldChange = Math.log(foldChange) / Math.log(2);
		// ratios with threshold or 2* standard deviation

		// use value defined in thresholdForSignificance

		if (threeCombined.contains(null) || threeCombined.contains(Double.NaN)) {
			// unclassified
			classification2Cases.put(sharedPeptidesProteinKey, Classification2Case.CASE6);
			// System.out.println("unclassified");
		}
		// make sure the unique1 is larger or equal to unique2
		else {

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

			final StringBuilder logString = new StringBuilder();
			// make sure order is correct
			if (Double.compare(unique1Ratio, unique2Ratio) > 0) {

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
					// log2ThresholdFoldChange
					if (Double.compare(distanceBetweenRatio1Ratio2, log2ThresholdFoldChange) > 0) {
						if (params.isStatisticalTestForProteinPairApplied()) {
							try {
								classifyPairByclassification1(threeCombinedOriginal, pepRatProt1, pepRatShared,
										pepRatProt2, sharedPeptidesProteinKey, cond1, cond2);
								final Classification1Case case1 = classification1Cases.get(sharedPeptidesProteinKey);
								if (case1 == Classification1Case.CASE3 || case1 == Classification1Case.CASE6
										|| case1 == Classification1Case.CASE7) {
									classification2Cases.put(sharedPeptidesProteinKey, Classification2Case.CASE3);
								} else {
									classification2Cases.put(sharedPeptidesProteinKey, Classification2Case.CASE6);
								}
								return;
							} catch (final NotEnoughMeasurementsException e) {
								// log.info("asdf");
							}
						}
						// check if shared value is larger than first unique
						// plus half log2ThresholdFoldChange or shared value is
						// smaller than
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
							}
							if (Double.compare(sharedRatio, unique2Ratio - distanceBetweenRatio1Ratio2 / 4) < 0) {
								if (swapped) {
									uniquePeptidesProt1Inconsistent = true;
								} else {
									uniquePeptidesProt2Inconsistent = true;
								}
							}
							if (unique1Ratio.isInfinite()) {
								if (swapped) {
									uniquePeptidesProt2Inconsistent = true;
								} else {
									uniquePeptidesProt1Inconsistent = true;
								}
							}
							if (unique2Ratio.isInfinite()) {
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
					if ((Double.compare(sharedRatio, unique1Ratio - log2ThresholdFoldChange) < 0)
							|| Double.compare(sharedRatio, unique2Ratio + log2ThresholdFoldChange) > 0) {
						final boolean isSharedPeptideSharedByOtherProteinOutOfThePair = isSharedPeptideNodeSharedByOtherProteinNodeOutOfThePair(
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
						if (params.isStatisticalTestForProteinPairApplied()) {
							try {
								classifyPairByclassification1(threeCombinedOriginal, pepRatProt1, pepRatShared,
										pepRatProt2, sharedPeptidesProteinKey, cond1, cond2);
								final Classification1Case case1 = classification1Cases.get(sharedPeptidesProteinKey);
								if (case1 == Classification1Case.CASE3 || case1 == Classification1Case.CASE6
										|| case1 == Classification1Case.CASE7) {
									classification2Cases.put(sharedPeptidesProteinKey, Classification2Case.CASE3);
								} else {
									classification2Cases.put(sharedPeptidesProteinKey, Classification2Case.CASE6);
								}
								return;
							} catch (final NotEnoughMeasurementsException e) {
								log.debug(e.getMessage());
							}
						}
						// check if difference of unique is larger than fold
						// change
						if (Double.compare(distanceBetweenRatio1Ratio2, log2ThresholdFoldChange) > 0) {
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
							classification2Cases.put(sharedPeptidesProteinKey, Classification2Case.CASE5);
							logString.append("no_difference_2\t");
						}

					}
				}
			}

			// log.debug(logString.toString());

		}

	}

	/**
	 * Classifies all three ratios into the three different cases based on
	 * infinities and using an statistical test when possible
	 *
	 * @param threeCombined
	 * @param pepRatProt2
	 * @param pepRatShared
	 * @param pepRatProt1
	 * @param minRatiosToTest
	 * @param gAM
	 * @throws IOException
	 */
	private void classifyPairByclassification1(List<Double> threeCombined, List<Double> pepRatProt1,
			List<Double> pepRatShared, List<Double> pepRatProt2, String sharedPeptidesProteinKey, QuantCondition cond1,
			QuantCondition cond2) throws IOException, NotEnoughMeasurementsException {

		final boolean hasRatiosAndINF = PCQUtils.hasRatiosAndINF(threeCombined);
		Classification1Case classification1Case = Classification1Case.CASE4;

		if (threeCombined.contains(null) || threeCombined.contains(Double.NaN)) {
			classification1Case = Classification1Case.CASE5;
		} else if (hasRatiosAndINF) {
			classification1Case = Classification1Case.CASE2;

		} else {
			// at this point we have INF OR RATIOS, not both
			final List<Double> INFRatios = PCQUtils.getINFRatiosValues(threeCombined);
			final List<Double> nonINFRatios = PCQUtils.getNonINFRatiosValues(threeCombined);

			// if we only have RATIO :
			if (!nonINFRatios.isEmpty()) {
				final int uniquePeptideOutlier = outlierTestUniquePeptides(threeCombined, pepRatProt1, pepRatShared,
						pepRatProt2);
				boolean sharedPeptideOutlier = outlierTest(threeCombined.get(1), pepRatProt1, pepRatProt2);
				if (uniquePeptideOutlier != 0) { // there are some outliers
					if (uniquePeptideOutlier == 1 || uniquePeptideOutlier == 3) {
						uniquePeptidesProt1Inconsistent = true;
					}
					if (uniquePeptideOutlier == 2 || uniquePeptideOutlier == 3) {
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
					final boolean isSharedPeptideSharedByOtherProteinOutOfThePair = isSharedPeptideNodeSharedByOtherProteinNodeOutOfThePair(
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

			}
			// if we only have INF or -INF:
			else if (!INFRatios.isEmpty()) {
				final Double infinities = PCQUtils.areAllINFValuesSame(INFRatios);
				if (infinities == null) {
					classification1Case = Classification1Case.CASE1;
					uniquePeptidesProt1Inconsistent = true;
					uniquePeptidesProt2Inconsistent = true;
					sharedPeptidesInconsistent = true;
				} else {
					classification1Case = Classification1Case.CASE4;
				}

			}
		}

		classification1Cases.put(sharedPeptidesProteinKey, classification1Case);
	}

	private boolean isSharedPeptideNodeSharedByOtherProteinNodeOutOfThePair(QuantCondition cond1,
			QuantCondition cond2) {
		final Map<String, Set<PCQPeptideNode>> sharedPeptides = PCQUtils
				.getSharedPeptideNodesByProteinNode(proteinNode1, proteinNode2, false, true);
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
	 * @throws NotEnoughMeasurementsException
	 */
	private int outlierTestUniquePeptides(List<Double> threeCombined, List<Double> pepRatProt1,
			List<Double> pepRatShared, List<Double> pepRatProt2) throws NotEnoughMeasurementsException {
		final List<Double> populationValues = new ArrayList<Double>();
		double valueToTest;
		int ret = 0;
		// first case
		valueToTest = threeCombined.get(0);
		if (outlierTest(valueToTest, pepRatShared, pepRatProt2)) {
			final double pValueFirstCase = Maths.iglewiczHoaglinTest(valueToTest, convertingList(populationValues));
			firstCase = new ProteinPairPValue(pValueFirstCase, this);
			ret = 1;
		}

		// third case
		valueToTest = threeCombined.get(2);
		if (outlierTest(valueToTest, pepRatProt1, pepRatShared)) {
			final double pValueSecondCase = Maths.iglewiczHoaglinTest(valueToTest, convertingList(populationValues));
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
	 * @throws NotEnoughMeasurementsException
	 */

	final private boolean outlierTest(Double valueToTest, List<Double>... populationSets)
			throws NotEnoughMeasurementsException {

		final List<Double> populationValues = new ArrayList<Double>();
		boolean outlier = false;

		populationValues.clear();
		populationValues.add(valueToTest);
		for (final List<Double> populationSet : populationSets) {
			for (final Double value : populationSet) {
				if (value != null && !Double.isInfinite(value) && !Double.isNaN(value)) {
					populationValues.add(value);
				}
			}
		}
		if (populationValues.size() >= MIN_RATIOS_TO_TEST) {
			final double iglewiczHoaglinTest = Maths.iglewiczHoaglinTest(valueToTest, convertingList(populationValues));
			if (iglewiczHoaglinTest >= ProteinClusterQuantParameters.getInstance().getIglewiczHoaglinTestThreshold()) {
				outlier = true;
			}
		} else {
			throw new NotEnoughMeasurementsException();
		}
		return outlier;
	}

	private double[] convertingList(List<Double> toBeConverted) {
		final double[] ret = new double[toBeConverted.size()];
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

	/**
	 * [PROTEIN1: UniqueTo1Peptides] [SharedPeptides] [PROTEIN2:
	 * UniqueTo2Peptides]
	 */
	@Override
	public String toString() {
		final StringBuilder sb = new StringBuilder();
		// protein1
		final List<PCQPeptideNode> uniquePeptideNodes1 = PCQUtils.getSortedPeptideNodesBySequence(
				PCQUtils.getUniquePeptideNodes(proteinNode1, proteinNode2, params.isUniquePepOnly(), true));
		final StringBuilder sb1 = new StringBuilder();
		for (final PCQPeptideNode quantifiedPeptideNode : uniquePeptideNodes1) {
			if (!"".equals(sb1.toString())) {
				sb1.append(",");
			}
			sb1.append(quantifiedPeptideNode.getSequence());
		}
		sb.append("[" + proteinNode1.getKey() + ": " + sb1.toString() + "] ");
		// shared
		// protein1
		final String shared = PCQUtils.getPeptideNodesSequenceString(
				PCQUtils.getSharedPeptideNodeSet(proteinNode1, proteinNode2, false, true));
		sb.append("[SHARED: " + shared + "] ");
		// protein2
		final List<PCQPeptideNode> uniquePeptideNodes2 = PCQUtils.getSortedPeptideNodesBySequence(
				PCQUtils.getUniquePeptideNodes(proteinNode2, proteinNode1, params.isUniquePepOnly(), true));
		final StringBuilder sb2 = new StringBuilder();
		for (final PCQPeptideNode peptideNode : uniquePeptideNodes2) {
			if (!"".equals(sb2.toString())) {
				sb2.append(",");
			}
			sb2.append(peptideNode.getSequence());
		}
		sb.append("[" + proteinNode2.getKey() + ": " + sb2.toString() + "] ");
		return sb.toString();
	}

	/**
	 * @return the protein1
	 */
	public PCQProteinNode getProteinNode1() {
		return proteinNode1;
	}

	/**
	 * @return the protein2
	 */
	public PCQProteinNode getProteinNode2() {
		return proteinNode2;
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
		final StringBuilder sb = new StringBuilder();

		sb.append("Unique Pep1 (U1)\t").append("Ratio (U1)\t").append("Protein 1\t").append("Taxon_P1\t")
				.append("Shared Pep (S)\t").append("Ratio (S)\t").append("Protein 2\t").append("Taxon_P2\t")
				.append("Unique Pep2 (U2)\t").append("Ratio (U2)\t").append("ClassificationCase\t")
				.append("Incomplete pair and significant case 6");
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
		final List<String> ret = new ArrayList<String>();
		final Map<String, Set<PCQPeptideNode>> sharedPeptidesMap = PCQUtils.getSharedPeptideNodesMap(proteinNode1,
				proteinNode2, false, true);
		final String taxonomy1 = PCQUtils.getTaxonomyString(proteinNode1);

		final String taxonomy2 = PCQUtils.getTaxonomyString(proteinNode2);
		for (final String proteinKey : sharedPeptidesMap.keySet()) {
			final Set<PCQPeptideNode> sharedPeptideNodes = sharedPeptidesMap.get(proteinKey);

			final StringBuilder sb = new StringBuilder();
			// UniquePep1
			final Set<PCQPeptideNode> uniquePeptideNodes1 = PCQUtils.getUniquePeptideNodes(proteinNode1, proteinNode2,
					true, true);
			sb.append(PCQUtils.getPeptideNodesSequenceString(uniquePeptideNodes1) + "\t");
			// RatioU1
			Double value1 = Double.NaN;
			final QuantRatio consensusRatio1 = PCQUtils.getRepresentativeRatioForPeptideNodes(uniquePeptideNodes1,
					cond1, cond2, null, true);
			if (consensusRatio1 != null) {
				final Double log2Ratio = consensusRatio1.getLog2Ratio(cond1, cond2);
				if (log2Ratio != null) {
					value1 = log2Ratio;
				}
			}
			sb.append(PCQUtils.escapeInfinity(value1) + "\t");
			// Protein1
			sb.append(proteinNode1.getKey() + "\t");
			// taxonomy of protein 1
			sb.append(taxonomy1 + "\t");
			// Shared peptide
			sb.append(PCQUtils.getPeptideNodesSequenceString(sharedPeptideNodes) + "\t");
			// shared peptide ratio
			Double valueShared = Double.NaN;
			final QuantRatio consensusRatioShared = PCQUtils.getRepresentativeRatioForPeptideNodes(sharedPeptideNodes,
					cond1, cond2, null, true);
			if (consensusRatioShared != null) {
				final Double log2Ratio = consensusRatioShared.getLog2Ratio(cond1, cond2);
				if (log2Ratio != null) {
					valueShared = log2Ratio;
				}
			}
			sb.append(PCQUtils.escapeInfinity(valueShared) + "\t");
			// Protein2
			sb.append(proteinNode2.getKey() + "\t");
			// taxonomy of protein 2
			sb.append(taxonomy2 + "\t");
			// UniquePep2
			final Set<PCQPeptideNode> uniquePeptideNodes2 = PCQUtils.getUniquePeptideNodes(proteinNode2, proteinNode1,
					true, true);
			sb.append(PCQUtils.getPeptideNodesSequenceString(uniquePeptideNodes2) + "\t");
			// RatioU2
			Double value2 = Double.NaN;
			final QuantRatio consensusRatio2 = PCQUtils.getRepresentativeRatioForPeptideNodes(uniquePeptideNodes2,
					cond1, cond2, null, true);
			if (consensusRatio2 != null) {
				final Double log2Ratio = consensusRatio2.getLog2Ratio(cond1, cond2);
				if (log2Ratio != null) {
					value2 = log2Ratio;
				}
			}
			sb.append(PCQUtils.escapeInfinity(value2) + "\t");
			// classification case 1
			final Classification1Case classification1Case = classification1Cases.get(proteinKey);
			if (classification1Case != null) {
				sb.append(classification1Case.name());
			}
			sb.append("\t");
			// classification case 2
			final Classification2Case classification2Case = classification2Cases.get(proteinKey);
			if (classification2Case != null) {
				sb.append(classification2Case.name());
			}

			// show a flag if we have case 6, that is, shared peptide node is
			// significantly different
			if (classification1Case != null && classification1Case == Classification1Case.CASE6
					&& Double.isNaN(value2)) {
				if (Math.abs(valueShared - value1) >= params.getThresholdForSignificance()) {
					sb.append("\tX");
				}
			}

			ret.add(sb.toString());
		}

		return ret;
	}

	/**
	 * @return the containsDiscardedProteinNode
	 */
	public boolean isContainsDiscardedProteinNode() {
		return containsDiscardedProteinNode;
	}

	public String getProtein1Protein2Key() {
		final List<String> list = new ArrayList<String>();
		list.add(proteinNode1.getKey());
		list.add(proteinNode2.getKey());
		Collections.sort(list);
		final StringBuilder sb = new StringBuilder();
		for (final String acc : list) {
			sb.append(acc);
		}
		return sb.toString();
	}

	public int getNumSharedNodes() {
		return PCQUtils.getSharedPeptideNodesByProteinNode(proteinNode1, proteinNode2, false, true).size();
	}
}
