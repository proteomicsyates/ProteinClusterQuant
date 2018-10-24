package edu.scripps.yates.pcq.compare.model;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import edu.scripps.yates.utilities.sequence.PositionInPeptide;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.hash.THashSet;

public class QuantifiedSite {

	public static final String NODE_KEY = "Node key";
	public static final String SEQUENCE = "Sequence";
	public static final String LOG2RATIO = "Log2Ratio";
	public static final String RATIOSCOREVALUE = "Ratio Score Value";
	public static final String NUMPSMS = "Num PSMs";
	public static final String NUMPEPTIDES = "Num Peptides";
	public static final String NUMMEASUREMENTS = "Num Measurements";
	public static final String PROTEINS = "Protein(s)";
	public static final String GENES = "Genes";
	public static final String POSITIONSINPEPTIDE = "QuantSitePositionInPeptide";

	private final String nodeKey;
	private final String sequence;
	private final Double log2Ratio;
	private Double log2Ratio2;
	private final Double ratioScoreValue;
	private Double ratioScoreValue2;
	private final int numPSMs;
	private int numPSMs2;
	private final int numPeptides;
	private int numPeptides2;
	private final int numMeasurements;
	private int numMeasurements2;
	private final Set<PositionInPeptide> positionsInPeptide = new THashSet<PositionInPeptide>();
	private final String proteins;
	private final String genes;

	public QuantifiedSite(String[] split, TObjectIntHashMap<String> indexesByHeaders) {
		nodeKey = split[indexesByHeaders.get(NODE_KEY)];
		sequence = split[indexesByHeaders.get(SEQUENCE)];
		String numberString = split[indexesByHeaders.get(LOG2RATIO)];
		if (numberString.startsWith("'")) {
			numberString = numberString.substring(1);
		}
		if (!"".equals(numberString)) {
			log2Ratio = Double.valueOf(numberString);
		} else {
			log2Ratio = null;
		}
		String ratioScoreValueString = split[indexesByHeaders.get(RATIOSCOREVALUE)];
		if (ratioScoreValueString.startsWith("'")) {
			ratioScoreValueString = ratioScoreValueString.substring(1);
		}
		if (!"".equals(ratioScoreValueString)) {
			ratioScoreValue = Double.valueOf(ratioScoreValueString);
		} else {
			ratioScoreValue = null;
		}
		numPSMs = Integer.valueOf(split[indexesByHeaders.get(NUMPSMS)]);
		numPeptides = Integer.valueOf(split[indexesByHeaders.get(NUMPEPTIDES)]);
		numMeasurements = Integer.valueOf(split[indexesByHeaders.get(NUMMEASUREMENTS)]);

		final String positionsInPeptideString = split[indexesByHeaders.get(POSITIONSINPEPTIDE)];
		positionsInPeptide.addAll(PositionInPeptide.parseStringToPositionInPeptide(positionsInPeptideString, "-"));
		proteins = split[indexesByHeaders.get(PROTEINS)];
		genes = split[indexesByHeaders.get(GENES)];
	}

	public List<PositionInPeptide> getPositionsInPeptide() {
		final List<PositionInPeptide> list = new ArrayList<PositionInPeptide>();
		list.addAll(positionsInPeptide);
		list.sort(new Comparator<PositionInPeptide>() {

			@Override
			public int compare(PositionInPeptide o1, PositionInPeptide o2) {
				return Integer.compare(o1.getPosition(), o2.getPosition());
			}
		});
		return list;
	}

	public String getNodeKey() {
		return nodeKey;
	}

	public String getSequence() {
		return sequence;
	}

	public Double getLog2Ratio() {
		return log2Ratio;
	}

	public Double getRatioScoreValue() {
		return ratioScoreValue;
	}

	public int getNumPSMs() {
		return numPSMs;
	}

	public int getNumPeptides() {
		return numPeptides;
	}

	public Double getLog2Ratio2() {
		return log2Ratio2;
	}

	public void setLog2Ratio2(Double log2Ratio2) {
		this.log2Ratio2 = log2Ratio2;
	}

	public Double getRatioScoreValue2() {
		return ratioScoreValue2;
	}

	public void setRatioScoreValue2(Double ratioScoreValue2) {
		this.ratioScoreValue2 = ratioScoreValue2;
	}

	public int getNumPSMs2() {
		return numPSMs2;
	}

	public void setNumPSMs2(int numPSMs2) {
		this.numPSMs2 = numPSMs2;
	}

	public int getNumPeptides2() {
		return numPeptides2;
	}

	public void setNumPeptides2(int numPeptides2) {
		this.numPeptides2 = numPeptides2;
	}

	public int getNumMeasurements() {
		return numMeasurements;
	}

	public int getNumMeasurements2() {
		return numMeasurements2;
	}

	public void setNumMeasurements2(int numMeasurements2) {
		this.numMeasurements2 = numMeasurements2;
	}

	public String getProteins() {
		return proteins;
	}

	public String getGenes() {
		return genes;
	}
}
