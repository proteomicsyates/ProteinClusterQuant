package edu.scripps.yates.pcq.model;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.CensusRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.pcq.params.ProteinClusterQuantParameters;
import edu.scripps.yates.pcq.util.PCQUtils;
import edu.scripps.yates.utilities.model.enums.AggregationLevel;
import edu.scripps.yates.utilities.sequence.PositionInPeptide;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

/**
 * Wrapper class for a peptide node, that could contain more than one
 * {@link QuantifiedPSMInterface} and even different
 * {@link QuantifiedPeptideInterface}.<br>
 * It will contain a {@link QuantRatio} and a confidenceValue associated.
 *
 * @author Salva
 *
 */
public class PCQPeptideNode extends AbstractNode<QuantifiedPeptideInterface> {
	private final static Logger log = Logger.getLogger(PCQPeptideNode.class);

	public static final String INTEGRATED_PEPTIDE_NODE_RATIO = "Integrated peptide node ratio";

	private final Set<QuantifiedPeptideInterface> peptideSet = new THashSet<QuantifiedPeptideInterface>();
	private Double confidenceValue;
	private final Set<QuantRatio> sanXotRatios = new THashSet<QuantRatio>();

	private final Set<PCQProteinNode> proteinNodes = new THashSet<PCQProteinNode>();

	private final Map<String, QuantRatio> sanXotRatiosByReplicateName = new THashMap<String, QuantRatio>();

	private final ProteinCluster proteinCluster;

	private String key;

	private final Map<QuantifiedPeptideInterface, PositionInPeptide> positionInPeptideByPeptide = new THashMap<QuantifiedPeptideInterface, PositionInPeptide>();

	public PCQPeptideNode(ProteinCluster proteinCluster, Collection<QuantifiedPeptideInterface> peptideCollection) {
		peptideSet.addAll(peptideCollection);
		this.proteinCluster = proteinCluster;
	}

	public PCQPeptideNode(ProteinCluster proteinCluster, String key,
			Pair<QuantifiedPeptideInterface, PositionInPeptide>... peptidesAndPositionInPeptides) {
		for (final Pair<QuantifiedPeptideInterface, PositionInPeptide> pair : peptidesAndPositionInPeptides) {
			final QuantifiedPeptideInterface peptide = pair.getFirstelement();
			peptideSet.add(peptide);
			final PositionInPeptide positionInPeptide = pair.getSecondElement();
			positionInPeptideByPeptide.put(peptide, positionInPeptide);
		}
		this.proteinCluster = proteinCluster;
		this.key = key;
	}

	public PCQPeptideNode(ProteinCluster proteinCluster, QuantifiedPeptideInterface... peptides) {
		for (final QuantifiedPeptideInterface peptide : peptides) {
			peptideSet.add(peptide);
		}
		this.proteinCluster = proteinCluster;

	}

	public Set<PCQProteinNode> getProteinNodes() {
		return proteinNodes;
	}

	/**
	 * @return the confidenceValue
	 */
	public Double getConfidenceValue() {
		return confidenceValue;
	}

	/**
	 * @param confidenceValue
	 *            the confidenceValue to set
	 */
	public void setConfidenceValue(double confidenceValue) {
		this.confidenceValue = confidenceValue;
	}

	public String getSequence() {
		return PCQUtils.getPeptidesSequenceString(peptideSet);
	}

	public String getFullSequence() {
		return QuantUtils.getPeptidesFullSequenceString(peptideSet);
	}

	public Set<String> getTaxonomies() {
		final Set<String> set = new THashSet<String>();
		if (!ProteinClusterQuantParameters.getInstance().ignoreTaxonomies()) {
			final Set<QuantifiedProteinInterface> quantifiedProteins = getQuantifiedProteins();
			for (final QuantifiedProteinInterface protein : quantifiedProteins) {
				final Set<String> taxonomies = protein.getTaxonomies();
				set.addAll(taxonomies);
			}
		}
		return set;
	}

	@Override
	public String getKey() {
		if (key == null) {
			key = getFullSequence();
		}
		return key;
	}

	public Set<QuantRatio> getRatios() {
		final Set<QuantRatio> ratios = new THashSet<QuantRatio>();
		for (final QuantifiedPeptideInterface peptide : peptideSet) {
			ratios.addAll(peptide.getRatios());
			for (final QuantifiedPSMInterface psm : peptide.getQuantifiedPSMs()) {
				ratios.addAll(psm.getRatios());
			}
		}
		return ratios;
	}

	public QuantRatio getSanXotRatio(QuantCondition cond1, QuantCondition cond2) {
		// first look if there is an externally set consensus ratio that can be
		// come from the ratio integration analysis
		for (final QuantRatio ratio : sanXotRatios) {
			if (ratio.getCondition1().equals(cond1) && ratio.getCondition2().equals(cond2)) {
				return ratio;
			}
			if (ratio.getCondition1().equals(cond2) && ratio.getCondition2().equals(cond1)) {
				return ratio;
			}
		}
		// return Nan ratio
		final CensusRatio ratio = new CensusRatio(Double.NaN, false, cond1, cond2, AggregationLevel.PEPTIDE_NODE,
				INTEGRATED_PEPTIDE_NODE_RATIO);
		return ratio;

	}

	/**
	 * Called by ProteinClusterQuant.setIntegrationResultsIntoPeptideNodes
	 * 
	 * @param ratio
	 * @return
	 */
	public boolean addSanXotRatio(QuantRatio ratio) {
		return sanXotRatios.add(ratio);
	}

	/**
	 * Called by
	 * ProteinClusterQuant.setIntegrationResultsPerReplicateIntoPeptideNodes
	 * 
	 * @param ratio
	 * @param replicateName
	 * @return
	 */
	public boolean addSanXotRatio(QuantRatio ratio, String replicateName) {
		if (replicateName != null) {

			return sanXotRatiosByReplicateName.put(replicateName, ratio) != null;
		} else {
			return addSanXotRatio(ratio);
		}
	}

	@Override
	public Set<QuantifiedProteinInterface> getQuantifiedProteins() {
		final Set<QuantifiedProteinInterface> proteins = new THashSet<QuantifiedProteinInterface>();

		for (final PCQProteinNode proteinNode : getProteinNodes()) {
			proteins.addAll(proteinNode.getQuantifiedProteins());
		}
		return proteins;
	}

	@Override
	public Set<QuantifiedPSMInterface> getQuantifiedPSMs() {
		final Set<QuantifiedPSMInterface> ret = new THashSet<QuantifiedPSMInterface>();
		for (final QuantifiedPeptideInterface peptide : getQuantifiedPeptides()) {
			ret.addAll(peptide.getQuantifiedPSMs());
		}
		return ret;
	}

	@Override
	public Set<QuantifiedPeptideInterface> getQuantifiedPeptides() {
		return peptideSet;
	}

	public boolean addQuantifiedPeptide(QuantifiedPeptideInterface peptide) {
		return peptideSet.add(peptide);
	}

	public boolean addQuantifiedPeptide(QuantifiedPeptideInterface peptide, PositionInPeptide positionInPeptide) {
		positionInPeptideByPeptide.put(peptide, positionInPeptide);

		return peptideSet.add(peptide);

	}

	public void removePeptidesFromProteinsInNode() {
		final Iterator<QuantifiedPeptideInterface> peptideSetIterator = peptideSet.iterator();
		while (peptideSetIterator.hasNext()) {
			final QuantifiedPeptideInterface peptide = peptideSetIterator.next();
			QuantUtils.discardPeptide(peptide);
			if (peptide.getQuantifiedProteins().isEmpty() || peptide.getQuantifiedPSMs().isEmpty()) {
				peptideSetIterator.remove();
			}
		}
	}

	@Override
	public String toString() {

		return "{" + getKey() + "'}";
	}

	public boolean addProteinNode(PCQProteinNode proteinNode) {
		final boolean added = proteinNodes.add(proteinNode);
		return added;
	}

	public Set<QuantifiedPeptideInterface> getQuantifiedPeptidesInReplicate(String replicateName) {
		final Set<QuantifiedPeptideInterface> ret = new THashSet<QuantifiedPeptideInterface>();
		final Set<QuantifiedPeptideInterface> quantifiedPeptides = getQuantifiedPeptides();
		for (final QuantifiedPeptideInterface quantifiedPeptideInterface : quantifiedPeptides) {
			if (quantifiedPeptideInterface.getFileNames().contains(replicateName)) {
				ret.add(quantifiedPeptideInterface);
			}
		}
		return ret;
	}

	public QuantRatio getSanXotRatio(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator,
			String replicateName) {
		if (replicateName != null) {
			return sanXotRatiosByReplicateName.get(replicateName);
		} else {
			return getSanXotRatio(quantConditionNumerator, quantConditionDenominator);
		}
	}

	@Override
	public Set<QuantifiedPeptideInterface> getItemsInNode() {
		return peptideSet;
	}

	public ProteinCluster getCluster() {
		return proteinCluster;
	}

	public Set<QuantifiedProteinInterface> getNonDiscardedQuantifiedProteins() {
		final Set<QuantifiedProteinInterface> ret = new THashSet<QuantifiedProteinInterface>();
		for (final QuantifiedProteinInterface protein : getQuantifiedProteins()) {
			if (!protein.isDiscarded()) {
				ret.add(protein);
			}
		}
		return ret;
	}

	public boolean containsPTMs() {
		for (final QuantifiedPeptideInterface peptide : getItemsInNode()) {
			if (peptide.containsPTMs()) {
				return true;
			}
		}
		return false;
	}

	/**
	 * For a given peptide, gets the {@link PositionInPeptide} for which this
	 * node was created.<br>
	 * Note that this only is used for experiments in which isCollapsedBySites()
	 * is TRUE.
	 * 
	 * @param peptide
	 * @return
	 */
	public PositionInPeptide getPositionInPeptide(QuantifiedPeptideInterface peptide) {
		if (!peptideSet.contains(peptide)) {
			throw new IllegalArgumentException(
					"The peptide " + peptide.getSequence() + " (" + peptide.getKey() + ") is not in this peptide node");
		}
		if (positionInPeptideByPeptide.containsKey(peptide)) {
			return positionInPeptideByPeptide.get(peptide);
		} else {
			return null;
		}
	}

	public List<Pair<QuantifiedPeptideInterface, PositionInPeptide>> getPeptidesWithPositionsInPeptide() {
		final List<Pair<QuantifiedPeptideInterface, PositionInPeptide>> peptidesAndPositionsInPeptides = new ArrayList<Pair<QuantifiedPeptideInterface, PositionInPeptide>>();
		final Set<QuantifiedPeptideInterface> peptides = getQuantifiedPeptides();
		for (final QuantifiedPeptideInterface peptide : peptides) {

			final Pair<QuantifiedPeptideInterface, PositionInPeptide> pair = new Pair<QuantifiedPeptideInterface, PositionInPeptide>(
					peptide, getPositionInPeptide(peptide));
			peptidesAndPositionsInPeptides.add(pair);

		}
		Collections.sort(peptidesAndPositionsInPeptides,
				new java.util.Comparator<Pair<QuantifiedPeptideInterface, PositionInPeptide>>() {

					@Override
					public int compare(Pair<QuantifiedPeptideInterface, PositionInPeptide> arg0,
							Pair<QuantifiedPeptideInterface, PositionInPeptide> arg1) {
						return arg0.getFirstelement().getFullSequence()
								.compareTo(arg1.getFirstelement().getFullSequence());
					}
				});
		return peptidesAndPositionsInPeptides;
	}

}
