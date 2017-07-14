package edu.scripps.yates.pcq.xgmml.util;

import java.util.Collections;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class AlignmentSet {
	private static final Logger log = Logger.getLogger(AlignmentSet.class);
	private final Map<QuantifiedPeptideInterface, Set<AlignedPeptides>> map = new THashMap<QuantifiedPeptideInterface, Set<AlignedPeptides>>();

	public void addAlignment(AlignedPeptides alignment) {
		addToMap(alignment.getPeptide1(), alignment);
		addToMap(alignment.getPeptide2(), alignment);
	}

	private void addToMap(QuantifiedPeptideInterface peptide, AlignedPeptides alignment) {
		if (map.containsKey(peptide)) {
			map.get(peptide).add(alignment);
		} else {
			Set<AlignedPeptides> set = new THashSet<AlignedPeptides>();
			set.add(alignment);
			map.put(peptide, set);
		}
	}

	/**
	 * Returns the set of {@link AlignedPeptides} for a particular
	 * {@link QuantifiedPeptideInterface}, or an empty set if there is non
	 *
	 * @param peptide
	 * @return
	 */
	public Set<AlignedPeptides> getAlignmentsForPeptide(QuantifiedPeptideInterface peptide) {
		if (map.containsKey(peptide)) {
			return map.get(peptide);
		}
		return Collections.emptySet();
	}
}
