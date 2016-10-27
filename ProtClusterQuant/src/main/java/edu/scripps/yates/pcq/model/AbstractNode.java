package edu.scripps.yates.pcq.model;

import java.util.HashSet;
import java.util.Set;

import edu.scripps.yates.census.read.model.interfaces.HasKey;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedItem;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;

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

		for (T item : getItemsInNode()) {
			item.setDiscarded(b);
		}
	}

	@Override
	public boolean isDiscarded() {
		return discarded;
	}

	@Override
	public Set<String> getFileNames() {
		Set<String> ret = new HashSet<String>();
		for (T t : getItemsInNode()) {
			ret.addAll(t.getFileNames());
		}
		return ret;
	}

	@Override
	public Set<String> getRawFileNames() {
		Set<String> rawFileNames = new HashSet<String>();
		for (T t : getItemsInNode()) {
			rawFileNames.addAll(t.getRawFileNames());
		}
		return rawFileNames;
	}

	public Set<T> getNonDiscardedItems() {
		Set<T> ret = new HashSet<T>();
		for (T t : getItemsInNode()) {
			if (!t.isDiscarded()) {
				ret.add(t);
			}
		}
		return ret;

	}

	@Override
	public boolean isQuantified() {

		for (T t : getItemsInNode()) {
			if (t.isQuantified()) {
				return true;
			}
		}
		return false;
	}
}
