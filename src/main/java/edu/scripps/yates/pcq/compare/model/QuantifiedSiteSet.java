package edu.scripps.yates.pcq.compare.model;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class QuantifiedSiteSet extends THashSet<QuantifiedSite> {
	private final THashMap<String, QuantifiedSite> quantifiedSitesByKey = new THashMap<String, QuantifiedSite>();
	private final List<String> sampleNames = new ArrayList<String>();

	@Override
	public boolean add(QuantifiedSite quantSite) {
		final boolean ret = super.add(quantSite);
		quantifiedSitesByKey.put(quantSite.getNodeKey(), quantSite);
		return ret;
	}

	public THashMap<String, QuantifiedSite> getQuantifiedSitesByKey() {
		return quantifiedSitesByKey;
	}

	public void addSampleName(String sampleName) {
		sampleNames.add(sampleName);
	}

	public int getNumExperiments() {
		return iterator().next().getNumExperiments();
	}

	public List<QuantifiedSite> getSortedByRatios() {
		final List<QuantifiedSite> ret = new ArrayList<QuantifiedSite>();
		ret.addAll(this);
		final Comparator<QuantifiedSite> comparator = new Comparator<QuantifiedSite>() {

			@Override
			public int compare(QuantifiedSite o1, QuantifiedSite o2) {

				for (int i = 0; i < o1.getNumExperiments(); i++) {
					final Double num1 = o1.getLog2Ratio(i);
					final Double num2 = o2.getLog2Ratio(i);
					if (Double.isInfinite(num1) && !Double.isInfinite(num2)) {
						return 1;
					}

					if (Double.isInfinite(num2) && !Double.isInfinite(num1)) {
						return -1;
					}
					final int ret = Double.compare(num1, num2);
					if (ret != 0) {
						return ret;
					}
				}
				return 0;
			}
		};

		ret.sort(comparator);
		return ret;
	}

	public List<String> getSampleNames() {
		return sampleNames;
	}
}
