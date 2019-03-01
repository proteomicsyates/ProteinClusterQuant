package edu.scripps.yates.pcq.compare.model;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import edu.scripps.yates.utilities.sequence.PositionInPeptide;
import edu.scripps.yates.utilities.sequence.PositionInProtein;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.hash.THashSet;

/**
 * It represents a quantified site in a protein (or proteins) in a series of
 * experiments (that is why log ratios and other features are lists).<br>
 * The quantification of the site could come from the aggregation of multiple
 * individual measurements and that is why it has features such as
 * numMeasurements, numPeptides and numPSMs
 * 
 * @author salvador
 *
 */
public class QuantifiedSite {

	public static final String NODE_KEY = "Node key";
	public static final String SEQUENCE = "Sequence";
	public static final String LOG2RATIO = "Log2Ratio";
	public static final String RATIOSCOREVALUE = "Ratio Score Value";
	public static final String STDEV = "Ratio STDEV";
	public static final String NUMPSMS = "Num PSMs";
	public static final String NUMPEPTIDES = "Num Peptides";
	public static final String NUMMEASUREMENTS = "Num Measurements";
	public static final String PROTEINS = "Protein(s)";
	public static final String GENES = "Genes";
	public static final String QUANTPOSITIONSINPEPTIDE = "QuantSitePositionInPeptide";

	private final String nodeKey;
	private String sequence;
	private final TDoubleArrayList log2Ratio = new TDoubleArrayList();
	private final TDoubleArrayList ratioStdevs = new TDoubleArrayList();
	private final TIntArrayList numPSMs = new TIntArrayList();
	private final TIntArrayList numPeptides = new TIntArrayList();;
	private final TIntArrayList numMeasurements = new TIntArrayList();;
	private Set<PositionInPeptide> positionsInPeptide = new THashSet<PositionInPeptide>();
	private String proteins;
	private String genes;
	private final List<PositionInProtein> positionInProteinList;
	private TTestMatrix ttestMatrix;

	public QuantifiedSite(String nodeKey) {
		this.nodeKey = nodeKey;
		positionInProteinList = PositionInProtein.parseStringToPositionInProtein(nodeKey, "-");
		sequence = null;
		log2Ratio.add(Double.NaN);
		ratioStdevs.add(Double.NaN);
		numPSMs.add(0);
		numPeptides.add(0);
		numMeasurements.add(0);
		proteins = null;
		genes = null;
	}

	public QuantifiedSite(String[] split, TObjectIntHashMap<String> indexesByHeaders) {
		nodeKey = split[indexesByHeaders.get(NODE_KEY)];
		positionInProteinList = PositionInProtein.parseStringToPositionInProtein(nodeKey, "-");

		sequence = split[indexesByHeaders.get(SEQUENCE)];
		String numberString = split[indexesByHeaders.get(LOG2RATIO)];
		if (numberString.startsWith("'")) {
			numberString = numberString.substring(1);
		}
		if (!"".equals(numberString)) {
			final Double num = Double.valueOf(numberString);
			log2Ratio.add(num);
		} else {
			log2Ratio.add(Double.NaN);
		}

		// stdev
		String stdevString = null;
		if (indexesByHeaders.contains(STDEV)) {
			stdevString = split[indexesByHeaders.get(STDEV)];
		} else {
			stdevString = split[indexesByHeaders.get(RATIOSCOREVALUE)];
		}
		if (stdevString != null && stdevString.startsWith("'")) {
			stdevString = stdevString.substring(1);
		}
		if (!"".equals(stdevString)) {
			ratioStdevs.add(Double.valueOf(stdevString));
		} else {
			ratioStdevs.add(Double.NaN);
		}

		numPSMs.add(Integer.valueOf(split[indexesByHeaders.get(NUMPSMS)]));
		numPeptides.add(Integer.valueOf(split[indexesByHeaders.get(NUMPEPTIDES)]));
		if (indexesByHeaders.containsKey(NUMMEASUREMENTS)) {
			numMeasurements.add(Integer.valueOf(split[indexesByHeaders.get(NUMMEASUREMENTS)]));
		} else {
			numMeasurements.add(0);
		}
		String positionsInPeptideString = null;
		if (indexesByHeaders.contains(QUANTPOSITIONSINPEPTIDE)) {
			positionsInPeptideString = split[indexesByHeaders.get(QUANTPOSITIONSINPEPTIDE)];
		}
		if (positionsInPeptideString != null) {
			positionsInPeptide.addAll(PositionInPeptide.parseStringToPositionInPeptide(positionsInPeptideString, "-"));
		}
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

	public Double getLog2Ratio(int index) {
		return log2Ratio.get(index);
	}

	public Double getRatioStdev(int index) {
		return ratioStdevs.get(index);
	}

	public void addRatioStdev(Double stdev) {
		ratioStdevs.add(stdev);
	}

	public int getNumPSMs(int index) {
		return numPSMs.get(index);
	}

	public int getNumMeasurements(int index) {
		return numMeasurements.get(index);
	}

	public int getNumPeptides(int index) {
		return numPeptides.get(index);
	}

	public void addLog2Ratio(Double log2Ratio2) {
		log2Ratio.add(log2Ratio2);
	}

	public void addNumPSMs(int numPSMs2) {
		numPSMs.add(numPSMs2);
	}

	public void addNumPeptides(int numPeptides2) {
		numPeptides.add(numPeptides2);
	}

	public void addNumMeasurements(int numMeasurements2) {
		numMeasurements.add(numMeasurements2);
	}

	public String getProteins() {
		return proteins;
	}

	public String getGenes() {
		return genes;
	}

	public int getNumExperiments() {
		return log2Ratio.size();
	}

	public String getPositions() {
		final StringBuilder sb = new StringBuilder();
		for (final PositionInProtein positionInProtein : positionInProteinList) {
			if (!"".equals(sb.toString())) {
				sb.append("-");
			}
			sb.append(String.valueOf(positionInProtein.getAa()) + positionInProtein.getPosition());
		}
		return sb.toString();
	}

	public List<PositionInProtein> getPositionInProteinList() {
		return positionInProteinList;
	}

	public TTestMatrix getTtestMatrix() {
		return ttestMatrix;
	}

	public void setTtestMatrix(TTestMatrix ttestMatrix) {
		this.ttestMatrix = ttestMatrix;
	}

	public void setProteins(String proteins2) {
		this.proteins = proteins2;
	}

	public void setGenes(String genes2) {
		this.genes = genes2;
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

	public void setPositionsInPeptide(List<PositionInPeptide> positionsInPeptide) {
		this.positionsInPeptide = new THashSet<PositionInPeptide>();
		positionsInPeptide.addAll(positionsInPeptide);
	}
}
