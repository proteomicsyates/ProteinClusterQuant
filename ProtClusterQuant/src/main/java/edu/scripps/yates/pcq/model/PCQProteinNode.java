package edu.scripps.yates.pcq.model;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.pcq.util.PCQUtils;

/**
 * Wrapper class for a protein node, that could contain several proteins sharing
 * the same set of {@link PCQPeptideNode}
 *
 * @author Salva
 *
 */
public class PCQProteinNode extends AbstractNode<QuantifiedProteinInterface> {
	// private final static Logger log = Logger.getLogger(PCQProteinNode.class);
	private final Set<QuantifiedProteinInterface> proteinSet = new HashSet<QuantifiedProteinInterface>();
	private final Set<PCQPeptideNode> peptideNodes = new HashSet<PCQPeptideNode>();
	private String accessionString;
	private String descriptionString;
	private HashSet<String> taxonomies;
	private final ProteinCluster proteinCluster;
	private ProteinPair proteinPair;

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

	@Override
	public Set<QuantifiedPeptideInterface> getQuantifiedPeptides() {
		Set<QuantifiedPeptideInterface> ret = new HashSet<QuantifiedPeptideInterface>();
		for (PCQPeptideNode peptideNode : getPeptideNodes()) {
			ret.addAll(peptideNode.getQuantifiedPeptides());
		}
		return ret;
	}

	public String getAccession() {
		if (accessionString == null) {
			accessionString = PCQUtils.getAccessionString(getQuantifiedProteins());
		}
		return accessionString;
	}

	public String getKey() {
		return getAccession();
	}

	public String getDescription() {
		if (descriptionString == null) {
			descriptionString = PCQUtils.getDescriptionStringFromIndividualProteins(proteinSet, true);
		}
		return descriptionString;
	}

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

	public void addProtein(QuantifiedProteinInterface protein) {
		proteinSet.add(protein);
		resetStrings();

	}

	private void resetStrings() {
		accessionString = null;
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

	@Override
	public Set<QuantifiedProteinInterface> getQuantifiedProteins() {
		return proteinSet;
	}

	@Override
	public Set<QuantifiedProteinInterface> getItemsInNode() {
		return proteinSet;
	}

	public void removeProteinsFromPeptidesInNode() {
		for (QuantifiedProteinInterface protein : proteinSet) {
			QuantUtils.discardProtein(protein);
		}
	}

	/**
	 * @return the proteinCluster
	 */
	public ProteinCluster getProteinCluster() {
		return proteinCluster;
	}

	public void setProteinPair(ProteinPair proteinPair) {
		this.proteinPair = proteinPair;
	}

	public ProteinPair getProteinPair() {
		return proteinPair;
	}
}
