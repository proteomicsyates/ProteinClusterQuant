package edu.scripps.yates.pcq.model;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.annotations.uniprot.xml.Entry;
import edu.scripps.yates.annotations.util.UniprotEntryUtil;
import edu.scripps.yates.census.read.model.QuantifiedProtein;
import edu.scripps.yates.census.read.model.QuantifiedProteinFromDBIndexEntry;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.pcq.params.ProteinClusterQuantParameters;
import edu.scripps.yates.pcq.util.PCQUtils;
import gnu.trove.set.hash.THashSet;

/**
 * Wrapper class for a protein node, that could contain several proteins sharing
 * the same set of {@link PCQPeptideNode}
 *
 * @author Salva
 *
 */
public class PCQProteinNode extends AbstractNode<QuantifiedProteinInterface> {
	// private final static Logger log = Logger.getLogger(PCQProteinNode.class);
	private final Set<QuantifiedProteinInterface> proteinSet = new THashSet<QuantifiedProteinInterface>();
	private final List<PCQPeptideNode> peptideNodes = new ArrayList<PCQPeptideNode>();
	private Set<String> taxonomies;
	private final ProteinCluster proteinCluster;
	private ProteinPair proteinPair;
	private String key;

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
		for (final QuantifiedProteinInterface protein : proteinCollection) {
			proteinSet.add(protein);
		}
		this.proteinCluster = proteinCluster;
	}

	public List<PCQPeptideNode> getPeptideNodes() {
		return peptideNodes;
	}

	@Override
	public Set<QuantifiedPeptideInterface> getQuantifiedPeptides() {
		final Set<QuantifiedPeptideInterface> ret = new THashSet<QuantifiedPeptideInterface>();
		for (final PCQPeptideNode peptideNode : getPeptideNodes()) {
			ret.addAll(peptideNode.getQuantifiedPeptides());
		}
		return ret;
	}

	private String getAccession() {

		return PCQUtils.getAccessionString(getQuantifiedProteins());

	}

	@Override
	public String getKey() {
		if (key == null) {
			key = PCQUtils.getAccessionString(getQuantifiedProteins());
		}
		return key;
	}

	public String getDescription() {

		return PCQUtils.getDescriptionStringFromIndividualProteins(proteinSet, true);

	}

	/**
	 * Sets description and taxonomy of proteins if they don't have it (it could
	 * be because they have been read from a TSV file with just Accession)
	 * 
	 * @param annotatedProteins
	 */
	public void annotateProteins(Map<String, Entry> annotatedProteins) {
		for (final QuantifiedProteinInterface protein : proteinSet) {
			if (annotatedProteins.containsKey(protein.getAccession())) {
				if (protein instanceof QuantifiedProtein) {

					final QuantifiedProtein quantprotein = (QuantifiedProtein) protein;
					final Entry entry = annotatedProteins.get(protein.getAccession());
					// if (protein.getDescription() == null ||
					// "".equals(protein.getDescription())) {
					quantprotein.setDescription(UniprotEntryUtil.getProteinDescription(entry));
					// }
					// if (protein.getTaxonomies() == null ||
					// protein.getTaxonomies().isEmpty()) {
					quantprotein.setTaxonomy(UniprotEntryUtil.getTaxonomy(entry));
					// }
				} else if (protein instanceof QuantifiedProteinFromDBIndexEntry) {
					final QuantifiedProteinFromDBIndexEntry quantprotein = (QuantifiedProteinFromDBIndexEntry) protein;
					final Entry entry = annotatedProteins.get(protein.getAccession());
					// if (protein.getDescription() == null ||
					// "".equals(protein.getDescription())) {
					quantprotein.setDescription(UniprotEntryUtil.getProteinDescription(entry));
					// }
					// if (protein.getTaxonomies() == null ||
					// protein.getTaxonomies().isEmpty()) {
					quantprotein.setTaxonomy(UniprotEntryUtil.getTaxonomy(entry));
					// }
				}
			}
		}
	}

	public Set<String> getTaxonomies() {
		if (!ProteinClusterQuantParameters.getInstance().ignoreTaxonomies()) {
			if (taxonomies == null) {
				final List<String> sortedTaxonomies = PCQUtils.getSortedTaxonomies(proteinSet);
				for (final String taxonomy : sortedTaxonomies) {
					if (taxonomy != null) {
						if (taxonomies == null) {
							taxonomies = new THashSet<String>();
						}
						taxonomies.add(taxonomy);
					}
				}
			}
		}
		return taxonomies;
	}

	@Override
	public Set<QuantifiedPSMInterface> getQuantifiedPSMs() {
		final Set<QuantifiedPSMInterface> set = new THashSet<QuantifiedPSMInterface>();
		for (final QuantifiedProteinInterface protein : proteinSet) {
			set.addAll(protein.getQuantifiedPSMs());
		}
		return set;
	}

	public void addProtein(QuantifiedProteinInterface protein) {
		proteinSet.add(protein);

	}

	public void addProteins(Collection<QuantifiedProteinInterface> proteins) {
		proteinSet.addAll(proteins);

	}

	@Override
	public String toString() {
		final StringBuilder sb = new StringBuilder();
		sb.append(getKey() + "': ");
		final List<PCQPeptideNode> list = new ArrayList<PCQPeptideNode>();
		list.addAll(getPeptideNodes());
		Collections.sort(list, new Comparator<PCQPeptideNode>() {

			@Override
			public int compare(PCQPeptideNode o1, PCQPeptideNode o2) {
				return o1.getSequence().compareTo(o2.getSequence());
			}
		});
		final StringBuilder sb2 = new StringBuilder();
		for (final PCQPeptideNode peptideNode : list) {
			if (!"".equals(sb2.toString()))
				sb2.append(",");
			sb2.append(peptideNode);
		}
		sb.append(sb2);
		return "{" + sb.toString() + "}";
	}

	public boolean addPeptideNode(PCQPeptideNode peptideNode) {
		if (peptideNodes.contains(peptideNode)) {
			return false;
		}
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
		for (final QuantifiedProteinInterface protein : proteinSet) {
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

	/**
	 * A {@link PCQProteinNode} is discarded when:
	 * <ul>
	 * <li>it is filtered out, OR</li>
	 * <li>all its {@link PCQPeptideNode}s are discarded</li>
	 * </ul>
	 */
	@Override
	public boolean isDiscarded() {
		if (super.isDiscarded()) {
			return true;
		} else {
			// return true if all peptide nodes are discarded
			final List<PCQPeptideNode> peptideNodes2 = getPeptideNodes();
			for (final PCQPeptideNode pcqPeptideNode : peptideNodes2) {
				if (!pcqPeptideNode.isDiscarded()) {
					return false;
				}
			}
			return true;
		}
	}

	public boolean containsPTMs() {
		for (final PCQPeptideNode peptideNode : getPeptideNodes()) {
			if (peptideNode.containsPTMs()) {
				return true;
			}
		}
		return false;
	}
}
