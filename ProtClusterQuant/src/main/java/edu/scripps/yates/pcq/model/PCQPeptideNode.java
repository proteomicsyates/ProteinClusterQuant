package edu.scripps.yates.pcq.model;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
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
import edu.scripps.yates.pcq.util.PCQUtils;
import edu.scripps.yates.utilities.model.enums.AggregationLevel;

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

	private final Set<QuantifiedPeptideInterface> peptideSet = new HashSet<QuantifiedPeptideInterface>();
	private Double confidenceValue;
	private final Set<QuantRatio> consensusRatios = new HashSet<QuantRatio>();

	private final Set<PCQProteinNode> proteinNodes = new HashSet<PCQProteinNode>();

	private final Map<String, QuantRatio> consensusRatiosByReplicateName = new HashMap<String, QuantRatio>();

	private final ProteinCluster proteinCluster;

	private String key;

	public PCQPeptideNode(ProteinCluster proteinCluster, Collection<QuantifiedPeptideInterface> peptideCollection) {
		peptideSet.addAll(peptideCollection);
		this.proteinCluster = proteinCluster;

	}

	public PCQPeptideNode(ProteinCluster proteinCluster, QuantifiedPeptideInterface... peptides) {
		for (QuantifiedPeptideInterface peptide : peptides) {
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
		Set<String> set = new HashSet<String>();
		final Set<QuantifiedProteinInterface> quantifiedProteins = getQuantifiedProteins();
		for (QuantifiedProteinInterface protein : quantifiedProteins) {
			final Set<String> taxonomies = protein.getTaxonomies();
			set.addAll(taxonomies);
		}
		return set;
	}

	@Override
	public String getKey() {
		if (key != null) {
			return key;
		}
		return getFullSequence();
	}

	@Override
	public void setKey(String key) {
		this.key = key;

	}

	public Set<QuantRatio> getRatios() {
		Set<QuantRatio> ratios = new HashSet<QuantRatio>();
		for (QuantifiedPeptideInterface peptide : peptideSet) {
			ratios.addAll(peptide.getRatios());
			for (QuantifiedPSMInterface psm : peptide.getQuantifiedPSMs()) {
				ratios.addAll(psm.getRatios());
			}
		}
		return ratios;
	}

	public QuantRatio getConsensusRatio(QuantCondition cond1, QuantCondition cond2) {
		// first look if there is an externally set consensus ratio that can be
		// come from the ratio integration analysis
		for (QuantRatio ratio : consensusRatios) {
			if (ratio.getCondition1().equals(cond1) && ratio.getCondition2().equals(cond2)) {
				return ratio;
			}
			if (ratio.getCondition1().equals(cond2) && ratio.getCondition2().equals(cond1)) {
				return ratio;
			}
		}
		// return Nan ratio
		CensusRatio ratio = new CensusRatio(Double.NaN, false, cond1, cond2, AggregationLevel.PEPTIDE_NODE,
				INTEGRATED_PEPTIDE_NODE_RATIO);
		return ratio;

	}

	public boolean addConsensusRatio(QuantRatio ratio) {
		return consensusRatios.add(ratio);
	}

	public boolean addConsensusRatio(QuantRatio ratio, String replicateName) {
		if (replicateName != null) {

			return consensusRatiosByReplicateName.put(replicateName, ratio) != null;
		} else {
			return addConsensusRatio(ratio);
		}
	}

	public Set<QuantRatio> getNonInfinityRatios() {
		Set<QuantRatio> ratios = new HashSet<QuantRatio>();
		for (QuantifiedPeptideInterface peptide : peptideSet) {
			ratios.addAll(peptide.getNonInfinityRatios());
			for (QuantifiedPSMInterface psm : peptide.getQuantifiedPSMs()) {
				ratios.addAll(psm.getNonInfinityRatios());
			}
		}
		return ratios;
	}

	@Override
	public Set<QuantifiedProteinInterface> getQuantifiedProteins() {
		Set<QuantifiedProteinInterface> proteins = new HashSet<QuantifiedProteinInterface>();

		for (PCQProteinNode proteinNode : getProteinNodes()) {
			proteins.addAll(proteinNode.getQuantifiedProteins());
		}
		return proteins;
	}

	@Override
	public Set<QuantifiedPSMInterface> getQuantifiedPSMs() {
		Set<QuantifiedPSMInterface> ret = new HashSet<QuantifiedPSMInterface>();
		for (QuantifiedPeptideInterface peptide : getQuantifiedPeptides()) {
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

	public void removePeptidesFromProteinsInNode() {
		final Iterator<QuantifiedPeptideInterface> peptideSetIterator = peptideSet.iterator();
		while (peptideSetIterator.hasNext()) {
			QuantifiedPeptideInterface peptide = peptideSetIterator.next();
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
		Set<QuantifiedPeptideInterface> ret = new HashSet<QuantifiedPeptideInterface>();
		final Set<QuantifiedPeptideInterface> quantifiedPeptides = getQuantifiedPeptides();
		for (QuantifiedPeptideInterface quantifiedPeptideInterface : quantifiedPeptides) {
			if (quantifiedPeptideInterface.getFileNames().contains(replicateName)) {
				ret.add(quantifiedPeptideInterface);
			}
		}
		return ret;
	}

	public QuantRatio getConsensusRatio(QuantCondition quantConditionNumerator,
			QuantCondition quantConditionDenominator, String replicateName) {
		if (replicateName != null) {
			return consensusRatiosByReplicateName.get(replicateName);
		} else {
			return getConsensusRatio(quantConditionNumerator, quantConditionDenominator);
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
		Set<QuantifiedProteinInterface> ret = new HashSet<QuantifiedProteinInterface>();
		for (QuantifiedProteinInterface protein : getQuantifiedProteins()) {
			if (!protein.isDiscarded()) {
				ret.add(protein);
			}
		}
		return ret;
	}

	public boolean containsPTMs() {
		for (QuantifiedPeptideInterface peptide : getItemsInNode()) {
			if (peptide.containsPTMs()) {
				return true;
			}
		}
		return false;
	}

}
