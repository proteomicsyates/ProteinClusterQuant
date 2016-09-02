package edu.scripps.yates.pcq.model;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.pcq.util.PCQUtils;
import edu.scripps.yates.utilities.grouping.GroupablePSM;
import edu.scripps.yates.utilities.grouping.ProteinEvidence;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;

/**
 * Wrapper class for a protein node, that could contain several proteins sharing
 * the same set of {@link PCQPeptideNode}
 *
 * @author Salva
 *
 */
public class PCQProteinNode implements QuantifiedProteinInterface {
	// private final static Logger log = Logger.getLogger(PCQProteinNode.class);
	private final Set<QuantifiedProteinInterface> proteinSet = new HashSet<QuantifiedProteinInterface>();
	private final Set<PCQPeptideNode> peptideNodes = new HashSet<PCQPeptideNode>();
	private final Set<String> fileNames = new HashSet<String>();
	private boolean discarded;
	private String accessionString;
	private String descriptionString;
	private String taxonomiesString;
	private HashSet<String> taxonomies;
	private final ProteinCluster proteinCluster;

	public PCQProteinNode(ProteinCluster proteinCluster, Collection<QuantifiedProteinInterface> proteinCollection) {
		proteinSet.addAll(proteinCollection);
		this.proteinCluster = proteinCluster;
	}

	public PCQProteinNode(ProteinCluster proteinCluster, Collection<QuantifiedProteinInterface> proteinCollection1,
			Collection<QuantifiedProteinInterface> proteinCollection2) {
		proteinSet.addAll(proteinCollection1);
		proteinSet.addAll(proteinCollection2);
		this.proteinCluster = proteinCluster;
	}

	public PCQProteinNode(ProteinCluster proteinCluster, QuantifiedProteinInterface... proteinCollection) {
		for (QuantifiedProteinInterface protein : proteinCollection) {
			proteinSet.add(protein);
		}
		this.proteinCluster = proteinCluster;
	}

	public Set<PCQPeptideNode> getPeptideNodes() {
		return peptideNodes;
	}

	/**
	 * Gets the set of {@link PCQPeptideNode} connected to the
	 * {@link PCQProteinNode}
	 */
	@Override
	public Set<QuantifiedPeptideInterface> getQuantifiedPeptides() {
		Set<QuantifiedPeptideInterface> ret = new HashSet<QuantifiedPeptideInterface>();
		for (PCQPeptideNode peptideNode : peptideNodes) {
			ret.add(peptideNode);
		}
		return ret;
	}

	@Override
	public List<GroupablePSM> getGroupablePSMs() {
		List<GroupablePSM> ret = new ArrayList<GroupablePSM>();
		final Set<QuantifiedPSMInterface> quantifiedPSMs = getQuantifiedPSMs();
		for (QuantifiedPSMInterface quantifiedPSMInterface : quantifiedPSMs) {
			ret.add(quantifiedPSMInterface);
		}
		return ret;
	}

	@Override
	public ProteinGroup getProteinGroup() {
		throw new UnsupportedOperationException();
	}

	@Override
	public int getDBId() {
		throw new UnsupportedOperationException();
	}

	@Override
	public String getAccession() {
		if (accessionString == null) {
			accessionString = PCQUtils.getAccessionString(proteinSet);
		}
		return accessionString;
	}

	@Override
	public void setEvidence(ProteinEvidence evidence) {
		throw new UnsupportedOperationException();
	}

	@Override
	public ProteinEvidence getEvidence() {
		throw new UnsupportedOperationException();
	}

	@Override
	public void setProteinGroup(ProteinGroup proteinGroup) {
		throw new UnsupportedOperationException();
	}

	@Override
	public Set<QuantRatio> getRatios() {
		throw new UnsupportedOperationException();
	}

	@Override
	public Set<QuantRatio> getNonInfinityRatios() {
		throw new UnsupportedOperationException();
	}

	@Override
	public void addRatio(QuantRatio ratio) {
		throw new UnsupportedOperationException();
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

	@Override
	public String getKey() {
		return getAccession();
	}

	@Override
	public String getDescription() {
		if (descriptionString == null) {
			descriptionString = PCQUtils.getDescriptionString(proteinSet, true);
		}
		return descriptionString;
	}

	@Override
	public Set<String> getTaxonomies() {
		if (taxonomies == null) {
			taxonomies = new HashSet<String>();
			taxonomies.addAll(PCQUtils.getSortedTaxonomies(proteinSet));

		}
		return taxonomies;
	}

	@Override
	public Set<QuantifiedPSMInterface> getQuantifiedPSMs() {
		Set<QuantifiedPSMInterface> set = new HashSet<QuantifiedPSMInterface>();
		for (QuantifiedProteinInterface protein : proteinSet) {
			set.addAll(protein.getQuantifiedPSMs());
		}
		return set;
	}

	@Override
	public void addPeptide(QuantifiedPeptideInterface peptide) {
		throw new UnsupportedOperationException();
	}

	@Override
	public String getAccessionType() {
		List<String> list = new ArrayList<String>();
		for (QuantifiedProteinInterface protein : proteinSet) {
			if (!list.contains(protein.getAccessionType())) {
				list.add(protein.getAccessionType());
			}
		}
		Collections.sort(list);
		StringBuilder sb = new StringBuilder();
		for (String accType : list) {
			if (!"".equals(sb.toString())) {
				sb.append(",");
			}
			sb.append(accType);
		}
		return sb.toString();
	}

	@Override
	public void addPSM(QuantifiedPSMInterface quantifiedPSM) {
		throw new UnsupportedOperationException();

	}

	@Override
	public void setTaxonomy(String taxonomy) {
		throw new UnsupportedOperationException();
	}

	@Override
	public Set<String> getRawFileNames() {
		Set<String> fileNames = new HashSet<String>();
		for (QuantifiedProteinInterface protein : proteinSet) {
			fileNames.addAll(protein.getRawFileNames());
		}
		return fileNames;
	}

	@Override
	public Integer getLength() {
		throw new UnsupportedOperationException();
	}

	@Override
	public void setAccession(String primaryAccession) {
		throw new UnsupportedOperationException();
	}

	@Override
	public void setDescription(String description) {
		throw new UnsupportedOperationException();
	}

	public void addProtein(QuantifiedProteinInterface protein) {
		proteinSet.add(protein);
		resetStrings();

	}

	private void resetStrings() {
		accessionString = null;
		taxonomiesString = null;
		descriptionString = null;
	}

	public void addProteins(Collection<QuantifiedProteinInterface> proteins) {
		proteinSet.addAll(proteins);
		resetStrings();
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(getAccession() + "': ");
		List<PCQPeptideNode> list = new ArrayList<PCQPeptideNode>();
		list.addAll(getPeptideNodes());
		Collections.sort(list, new Comparator<PCQPeptideNode>() {

			@Override
			public int compare(PCQPeptideNode o1, PCQPeptideNode o2) {
				return o1.getSequence().compareTo(o2.getSequence());
			}
		});
		StringBuilder sb2 = new StringBuilder();
		for (PCQPeptideNode quantifiedPeptide : list) {
			if (!"".equals(sb2.toString()))
				sb2.append(",");
			sb2.append(quantifiedPeptide);
		}
		sb.append(sb2);
		return "{" + sb.toString() + "}";
	}

	public boolean addPeptideNode(PCQPeptideNode peptideNode) {
		return peptideNodes.add(peptideNode);
	}

	public Set<QuantifiedProteinInterface> getQuantifiedProteins() {
		return proteinSet;
	}

	@Override
	public void addFileName(String fileName) {
		fileNames.add(fileName);

	}

	@Override
	public Set<String> getFileNames() {

		return fileNames;
	}

	@Override
	public QuantRatio getConsensusRatio(QuantCondition cond1, QuantCondition cond2) {
		throw new IllegalArgumentException("Protein nodes have no consensus ratios");
	}

	@Override
	public QuantRatio getConsensusRatio(QuantCondition quantConditionNumerator,
			QuantCondition quantConditionDenominator, String replicateName) {
		throw new IllegalArgumentException("Protein nodes have no consensus ratios");
	}

	public Set<QuantifiedProteinInterface> getIndividualProteins() {
		return proteinSet;
	}

	public void disconnectProteinsInNode() {
		for (QuantifiedProteinInterface protein : proteinSet) {
			PCQUtils.discardProtein(protein);
		}
	}

	@Override
	public void setDiscarded(boolean b) {
		discarded = b;

	}

	@Override
	public boolean isDiscarded() {
		return discarded;
	}

	/**
	 * @return the proteinCluster
	 */
	public ProteinCluster getProteinCluster() {
		return proteinCluster;
	}
}
