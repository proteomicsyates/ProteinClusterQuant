package edu.scripps.yates.pcq.compare.model;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class QuantifiedSiteSet extends THashSet<QuantifiedSite> {
	private final THashMap<String, QuantifiedSite> quantifiedSitesByKey = new THashMap<String, QuantifiedSite>();

	@Override
	public boolean add(QuantifiedSite quantSite) {
		final boolean ret = super.add(quantSite);
		quantifiedSitesByKey.put(quantSite.getNodeKey(), quantSite);
		return ret;
	}

	public THashMap<String, QuantifiedSite> getQuantifiedSitesByKey() {
		return quantifiedSitesByKey;
	}

	public List<QuantifiedSite> getSortedByRatios() {
		final List<QuantifiedSite> ret = new ArrayList<QuantifiedSite>();
		ret.addAll(this);
		final Comparator<QuantifiedSite> comparator = new Comparator<QuantifiedSite>() {

			@Override
			public int compare(QuantifiedSite o1, QuantifiedSite o2) {

				Double num1 = o1.getLog2Ratio();
				Double num2 = o2.getLog2Ratio();
				if (Double.isInfinite(num1) && !Double.isInfinite(num2)) {
					return 1;
				}

				if (Double.isInfinite(num2) && !Double.isInfinite(num1)) {
					return -1;
				}
				int ret = Double.compare(num1, num2);
				if (ret != 0) {
					return ret;
				}
				num1 = o1.getLog2Ratio2();
				num2 = o2.getLog2Ratio2();
				ret = Double.compare(num1, num2);
				return ret;
			}
		};

		ret.sort(comparator);
		return ret;
	}
}
