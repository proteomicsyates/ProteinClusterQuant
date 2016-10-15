package edu.scripps.yates.pcq.xgmml.util;

import java.awt.Color;
import java.util.HashMap;
import java.util.Map;

import edu.scripps.yates.pcq.ProteinClusterQuantParameters;
import edu.scripps.yates.utilities.colors.ColorGenerator;

public class ColorManager {
	private final Map<String, Color> colorsByTaxonomies = new HashMap<String, Color>();
	// private final Map<Classification1Case, Color> colorsByCases1 = new
	// HashMap<Classification1Case, Color>();
	// private final Map<Classification2Case, Color> colorsByCases2 = new
	// HashMap<Classification2Case, Color>();
	private Color highlightColor = Color.red;
	// private Color sharedNodeLabelColor = Color.GRAY;
	private Color multiTaxonomyColor;
	private Color alignedPeptidesEdgeColor;
	private static ColorManager instance;

	private ColorManager() {

	}

	public static ColorManager getInstance() {
		if (instance == null) {
			instance = new ColorManager();
		}
		return instance;
	}

	public void addColorByTaxonomy(String taxonomy, Color color) {
		colorsByTaxonomies.put(taxonomy, color);
	}

	public Color getColorForProteinTaxonomy(String taxonomy) {
		for (String taxonomies : colorsByTaxonomies.keySet()) {
			if (taxonomies != null && taxonomy != null && (taxonomies.toLowerCase().contains(taxonomy.toLowerCase())
					|| (taxonomy != null && taxonomy.toLowerCase().contains(taxonomies.toLowerCase())))) {
				return colorsByTaxonomies.get(taxonomies);
			}
		}
		return Color.WHITE;
	}

	// public Color getColorByClassification1Case(Classification1Case case2) {
	// if (colorsByCases1.containsKey(case2)) {
	// return colorsByCases1.get(case2);
	// }
	// return Color.CYAN;
	// }
	//
	// public void addColorByClassification1Case(Classification1Case case2,
	// Color color) {
	// colorsByCases1.put(case2, color);
	// }
	//
	// public Color getColorByClassification2Case(Classification2Case case2) {
	// if (colorsByCases2.containsKey(case2)) {
	// return colorsByCases2.get(case2);
	// }
	// return Color.CYAN;
	// }
	//
	// public void addColorByClassification2Case(Classification2Case case2,
	// Color color) {
	// colorsByCases2.put(case2, color);
	// }

	public void addColorByTaxonomy(String tax, String colorString) {
		final Color hex2Rgb = ColorGenerator.hex2Rgb(colorString);
		if (hex2Rgb != null) {
			addColorByTaxonomy(tax, hex2Rgb);
		}

	}

	/**
	 * @return the highlightColor
	 */
	public Color getHighlightColor() {
		if (ProteinClusterQuantParameters.getInstance().isRemarkSignificantPeptides()) {
			return highlightColor;
		} else {
			return Color.black;
		}
	}

	/**
	 * @param highlightColor
	 *            the highlightColor to set
	 */
	public void setHighlightColor(Color highlightColor) {
		if (highlightColor != null)
			this.highlightColor = highlightColor;
	}

	public void setMultiTaxonomyColor(String multiTaxonomyColor) {
		this.multiTaxonomyColor = ColorGenerator.hex2Rgb(multiTaxonomyColor);
	}

	/**
	 * @return the multiTaxonomyColor
	 */
	public Color getMultiTaxonomyColor() {
		return multiTaxonomyColor;
	}

	/**
	 * @param multiTaxonomyColor
	 *            the multiTaxonomyColor to set
	 */
	public void setMultiTaxonomyColor(Color multiTaxonomyColor) {
		this.multiTaxonomyColor = multiTaxonomyColor;
	}

	// /**
	// * @return the sharedNodeLabelColor
	// */
	// public Color getSharedNodeLabelColor() {
	// return sharedNodeLabelColor;
	// }
	//
	// /**
	// * @param sharedNodeLabelColor
	// * the sharedNodeLabelColor to set
	// */
	// public void setSharedNodeLabelColor(Color sharedNodeLabelColor) {
	// this.sharedNodeLabelColor = sharedNodeLabelColor;
	// }
	//
	// public void setSharedNodeLabelColor(String sharedNodeLabelColor2) {
	// sharedNodeLabelColor = hex2Rgb(sharedNodeLabelColor2);
	// }

	public void setAlignedPeptidesEdgeColor(String alignedPeptideEdgecolor) {
		alignedPeptidesEdgeColor = ColorGenerator.hex2Rgb(alignedPeptideEdgecolor);

	}

	/**
	 * @return the alignedPeptidesEdgeColor
	 */
	public Color getAlignedPeptidesEdgeColor() {
		return alignedPeptidesEdgeColor;
	}

	/**
	 * @param alignedPeptidesEdgeColor
	 *            the alignedPeptidesEdgeColor to set
	 */
	public void setAlignedPeptidesEdgeColor(Color alignedPeptidesEdgeColor) {
		this.alignedPeptidesEdgeColor = alignedPeptidesEdgeColor;
	}

	public Color getFillColorForDiscardedNode() {
		return Color.LIGHT_GRAY;
	}

	public Color getLabelColorForDiscardedNode() {
		return Color.gray;
	}

}
