package edu.scripps.yates.pcq.model;

import java.util.Set;

import edu.scripps.yates.census.read.model.interfaces.QuantifiedItem;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.utilities.proteomicsmodel.HasKey;
import gnu.trove.set.hash.THashSet;

public abstract class AbstractNode<T extends QuantifiedItem> implements QuantifiedItem, HasKey {
	private boolean discarded;

	public abstract Set<QuantifiedPSMInterface> getQuantifiedPSMs();

	public abstract Set<QuantifiedPeptideInterface> getQuantifiedPeptides();

	public abstract Set<QuantifiedProteinInterface> getQuantifiedProteins();

	public abstract Set<T> getItemsInNode();

	@Override
	public abstract String getKey();

	@Override
	public void setDiscarded(boolean b) {
		discarded = b;

		for (final T item : getItemsInNode()) {
			item.setDiscarded(b);
		}
	}

	@Override
	public boolean isDiscarded() {
		return discarded;
	}

	@Override
	public Set<String> getFileNames() {
		final Set<String> ret = new THashSet<String>();
		for (final T t : getItemsInNode()) {
			ret.addAll(t.getFileNames());
		}
		return ret;
	}

	@Override
	public Set<String> getRawFileNames() {
		final Set<String> rawFileNames = new THashSet<String>();
		for (final T t : getItemsInNode()) {
			rawFileNames.addAll(t.getRawFileNames());
		}
		return rawFileNames;
	}

	public Set<T> getNonDiscardedItems() {
		final Set<T> ret = new THashSet<T>();
		for (final T t : getItemsInNode()) {
			if (!t.isDiscarded()) {
				ret.add(t);
			}
		}
		return ret;

	}

	@Override
	public boolean isQuantified() {

		for (final T t : getItemsInNode()) {
			if (t.isQuantified()) {
				return true;
			}
		}
		return false;
	}

	// @Override
	// public int hashCode() {
	// if (hashCode == -1) {
	// hashCode = Objects.hash(getKey());
	// }
	// return hashCode;
	// }
}
