package edu.scripps.yates.pcq.quantsite.groups;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import edu.scripps.yates.pcq.compare.model.QuantifiedSite;
import edu.scripps.yates.pcq.compare.model.QuantifiedSiteSet;
import edu.scripps.yates.utilities.sequence.PositionInProtein;
import gnu.trove.map.hash.THashMap;

public class QuantSiteSubsetIterator implements Iterator<QuantifiedSite> {
	private final List<QuantifiedSite> quantSites = new ArrayList<QuantifiedSite>();
	private int index = 0;
	private final Map<String, Integer> proteinSites = new THashMap<String, Integer>();

	public QuantSiteSubsetIterator(QuantifiedSiteSet quantSites, List<String> subsetlist) {
		this.quantSites.addAll(quantSites);
		for (final String subset : subsetlist) {
			if (subset.contains("#")) {
				final String[] split = subset.split("#");
				proteinSites.put(split[0], Integer.valueOf(split[1]));
			} else {
				proteinSites.put(subset, null);
			}
		}
	}

	@Override
	public boolean hasNext() {
		while (index < quantSites.size() && !isNextValid()) {
			index++;
		}
		if (index < quantSites.size()) {
			return true;
		}
		return false;
	}

	@Override
	public QuantifiedSite next() {
		while (index < quantSites.size() && !isNextValid()) {
			index++;
		}

		if (index < quantSites.size()) {
			try {
				return quantSites.get(index);
			} finally {
				index++;
			}
		}
		return null;

	}

	private boolean isNextValid() {
		final QuantifiedSite next = quantSites.get(index);
		final List<PositionInProtein> matchingSites = next.getPositionInProteinList().stream()
				.filter(p -> proteinSites.containsKey(p.getProteinACC())).collect(Collectors.toList());
		if (matchingSites.isEmpty()) {
			return false;
		}
		// if here, we have a quant site of a protein that is in the subset
		for (final PositionInProtein positionInProtein : matchingSites) {
			final Integer position = proteinSites.get(positionInProtein.getProteinACC());
			if (position == null) {
				return true;
			}
			if (Integer.compare(position, positionInProtein.getPosition()) == 0) {
				return true;
			}
		}
		return false;
	}

}
