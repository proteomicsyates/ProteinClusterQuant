package edu.scripps.yates.pcq.model;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.IsobaricQuantifiedPeptide;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.pcq.util.PCQUtils;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;

/**
 * Wrapper class for a peptide node, that could contain more than one
 * {@link QuantifiedPSMInterface} and even different
 * {@link QuantifiedPeptideInterface}.<br>
 * It will contain a {@link QuantRatio} and a confidenceValue associated.
 *
 * @author Salva
 *
 */
public class PCQPeptideNode implements QuantifiedPeptideInterface {
	private final static Logger log = Logger.getLogger(PCQPeptideNode.class);

	private final Set<QuantifiedPeptideInterface> peptideSet = new HashSet<QuantifiedPeptideInterface>();
	private Double confidenceValue;
	private final Set<QuantRatio> consensusRatios = new HashSet<QuantRatio>();

	private final Set<PCQProteinNode> proteinNodes = new HashSet<PCQProteinNode>();

	private final Map<String, QuantRatio> consensusRatiosByReplicateName = new HashMap<String, QuantRatio>();

	private boolean discarded;

	public PCQPeptideNode(Collection<QuantifiedPeptideInterface> peptideCollection) {
		peptideSet.addAll(peptideCollection);
		for (QuantifiedPeptideInterface quantifiedPeptideInterface : peptideCollection) {

			if (quantifiedPeptideInterface.getSequence().equals("AALCAVHVIR")) {
				log.info(quantifiedPeptideInterface);
			}
		}
	}

	public PCQPeptideNode(QuantifiedPeptideInterface... peptides) {
		for (QuantifiedPeptideInterface peptide : peptides) {
			if (peptide.getSequence().equals("AALCAVHVIR")) {
				log.info(peptide);
			}
			peptideSet.add(peptide);
		}
	}

	public Set<PCQProteinNode> getPCQProteinNodes() {
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

	@Override
	public String getSequence() {
		return PCQUtils.getPeptidesSequenceString(peptideSet);
	}

	@Override
	public String getFullSequence() {
		return PCQUtils.getPeptidesFullSequenceString(peptideSet);
	}

	@Override
	public Float getCalcMHplus() {
		throw new UnsupportedOperationException();
	}

	@Override
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
	public Float getMHplus() {
		throw new UnsupportedOperationException();
	}

	@Override
	public String getKey() {
		return getSequence();
	}

	@Override
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

	@Override
	public QuantRatio getConsensusRatio(QuantCondition cond1, QuantCondition cond2) {
		for (QuantRatio ratio : consensusRatios) {
			if (ratio.getCondition1().equals(cond1) && ratio.getCondition2().equals(cond2)) {
				return ratio;
			}
			if (ratio.getCondition1().equals(cond2) && ratio.getCondition2().equals(cond1)) {
				return ratio;
			}
		}
		// instead of returning a null ratio, return the ion count ratio if
		// possible
		QuantRatio ionCountRatio = getIonCountRatio(cond1, cond2);
		if (ionCountRatio != null) {
			return ionCountRatio;
		}
		return null;
	}

	private QuantRatio getIonCountRatio(QuantCondition cond1, QuantCondition cond2) {
		Set<IsobaricQuantifiedPeptide> isobaricPeptides = new HashSet<IsobaricQuantifiedPeptide>();
		final Set<QuantifiedPeptideInterface> individualPeptides = getIndividualPeptides();
		for (QuantifiedPeptideInterface quantifiedPeptideInterface : individualPeptides) {
			if (quantifiedPeptideInterface instanceof IsobaricQuantifiedPeptide) {
				isobaricPeptides.add((IsobaricQuantifiedPeptide) quantifiedPeptideInterface);
			}
		}
		if (isobaricPeptides.isEmpty()) {
			return null;
		}
		return PCQUtils.getConsensusIonCountRatio(isobaricPeptides, cond1, cond2);
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

	@Override
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
	public void addRatio(QuantRatio ratio) {
		throw new IllegalArgumentException("Not allowed");

	}

	@Override
	public double getMeanRatios(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator) {
		throw new UnsupportedOperationException();
	}

	@Override
	public double getSTDRatios(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator) {
		throw new UnsupportedOperationException();
	}

	@Override
	public Set<Amount> getAmounts() {
		throw new UnsupportedOperationException();
	}

	@Override
	public void addAmount(Amount amount) {
		throw new UnsupportedOperationException();
	}

	/**
	 * Gets all the {@link PCQProteinNode} connected to this
	 * {@link PCQPeptideNode}
	 */
	@Override
	public Set<QuantifiedProteinInterface> getQuantifiedProteins() {
		Set<QuantifiedProteinInterface> proteins = new HashSet<QuantifiedProteinInterface>();

		for (PCQProteinNode proteinNode : proteinNodes) {
			proteins.add(proteinNode);
		}
		return proteins;
	}

	@Override
	public Set<QuantifiedPSMInterface> getQuantifiedPSMs() {
		Set<QuantifiedPSMInterface> ret = new HashSet<QuantifiedPSMInterface>();
		for (QuantifiedPeptideInterface peptide : peptideSet) {
			ret.addAll(peptide.getQuantifiedPSMs());
		}
		return ret;
	}

	@Override
	public Set<String> getRawFileNames() {
		Set<String> rawFileNames = new HashSet<String>();
		for (QuantifiedPSMInterface psm : getQuantifiedPSMs()) {
			rawFileNames.add(psm.getRawFileName());
		}
		return rawFileNames;
	}

	public Set<QuantifiedPeptideInterface> getQuantifiedPeptides() {
		return peptideSet;
	}

	@Override
	public boolean addQuantifiedPSM(QuantifiedPSMInterface psm) {
		throw new IllegalArgumentException("Not allowed");
	}

	public boolean addQuantifiedPeptide(QuantifiedPeptideInterface peptide) {
		return peptideSet.add(peptide);
	}

	public void disconnectPeptidesInNode() {
		for (QuantifiedPeptideInterface peptide : peptideSet) {
			PCQUtils.discardPeptide(peptide);
		}
	}

	@Override
	public String toString() {
		final List<QuantifiedPeptideInterface> quantifiedPeptides = PCQUtils.getSortedPeptidesBySequence(peptideSet);
		StringBuilder sb = new StringBuilder();
		for (QuantifiedPeptideInterface quantifiedPeptideInterface : quantifiedPeptides) {
			if (!"".equals(sb.toString())) {
				sb.append(",");
			}
			sb.append(quantifiedPeptideInterface.getSequence() + "'");
		}
		return "{" + getSequence() + "'}";
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

	@Override
	public void addFileName(String fileName) {
		log.error("addFileName not supported by Peptide Node");
	}

	@Override
	public Set<String> getFileNames() {
		Set<String> fileNames = new HashSet<String>();
		for (QuantifiedPSMInterface psm : getQuantifiedPSMs()) {
			fileNames.addAll(psm.getFileNames());
		}
		return fileNames;
	}

	@Override
	public QuantRatio getConsensusRatio(QuantCondition quantConditionNumerator,
			QuantCondition quantConditionDenominator, String replicateName) {
		if (replicateName != null) {
			return consensusRatiosByReplicateName.get(replicateName);
		} else {
			return getConsensusRatio(quantConditionNumerator, quantConditionDenominator);
		}
	}

	public Set<QuantifiedPeptideInterface> getIndividualPeptides() {
		return peptideSet;
	}

	@Override
	public void setDiscarded(boolean b) {
		discarded = b;

	}

	@Override
	public boolean isDiscarded() {
		return discarded;
	}

}
