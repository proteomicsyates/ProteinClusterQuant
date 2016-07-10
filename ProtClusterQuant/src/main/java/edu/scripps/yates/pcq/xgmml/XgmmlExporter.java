package edu.scripps.yates.proteinclusters.xgmml;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.PropertyException;

import org.apache.derby.iapi.services.io.FileUtil;
import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.xml.Entry;
import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.Ion;
import edu.scripps.yates.census.read.model.IsobaricQuantifiedPeptide;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.proteinclusters.Classification1Case;
import edu.scripps.yates.proteinclusters.Classification2Case;
import edu.scripps.yates.proteinclusters.ProteinCluster;
import edu.scripps.yates.proteinclusters.ProteinPair;
import edu.scripps.yates.proteinclusters.Ratio;
import edu.scripps.yates.proteinclusters.util.ColorManager;
import edu.scripps.yates.proteinclusters.util.ProteinClusterQuantParameters;
import edu.scripps.yates.proteinclusters.util.ProteinLabel;
import edu.scripps.yates.proteinclusters.util.Shape;
import edu.scripps.yates.proteinclusters.util.Utils;
import edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph;
import edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Att;
import edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Edge;
import edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Graphics;
import edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Node;
import edu.scripps.yates.proteinclusters.xgmml.jaxb.ObjectFactory;
import edu.scripps.yates.utilities.alignment.nwalign.NWResult;
import edu.scripps.yates.utilities.colors.ColorGenerator;

public class XgmmlExporter {
	private final static Logger log = Logger.getLogger(XgmmlExporter.class);
	private static JAXBContext context;
	private QuantCondition cond1;
	private QuantCondition cond2;
	private int edgeCounter = 0;
	private final static String STRING = "string";
	private final static String BOOLEAN = "boolean";
	private final static String LIST = "list";
	private final static ObjectFactory factory = new ObjectFactory();
	private final static DecimalFormat formatter = new DecimalFormat("#.#");
	private ColorManager colorManager;
	private final Map<String, Node> visitedPeptideKeys = new HashMap<String, Node>();
	private final Map<QuantifiedProteinInterface, Node> visitedProteins = new HashMap<QuantifiedProteinInterface, Node>();
	private final Map<String, Edge> edgesIDs = new HashMap<String, Edge>();
	private Map<String, Entry> annotatedProteins;
	private static final String PCQ_ID = "PCQ_ID";

	private enum BORDER_TYPE {
		SOLID, DASHED
	};

	private static JAXBContext getJAXBContext() throws JAXBException {
		if (context == null) {
			context = JAXBContext.newInstance(Graph.class);
		}
		return context;
	}

	// public File exportToGmmlFromProteinPairs(File outputFile, String label,
	// Collection<ProteinPair> proteinPairs,
	// QuantCondition condition1, QuantCondition condition2, ColorManager
	// colorManager) throws JAXBException {
	// resetVisitedKeys();
	// cond1 = condition1;
	// cond2 = condition2;
	// this.colorManager = colorManager;
	// Graph ret = initializeGraph(label);
	// if (proteinPairs != null) {
	// createNodesAndEdgesFromProteinCluster(proteinPairs, ret);
	// }
	// final File file = createFile(ret, outputFile);
	// fixHeader(file);
	// System.out.println("File created: " + file.getAbsolutePath());
	// return file;
	// }

	public File exportToGmmlFromProteinClusters(File outputFile, String label, Collection<ProteinCluster> clusters,
			QuantCondition condition1, QuantCondition condition2, ColorManager colorManager) throws JAXBException {

		resetVisitedKeys();
		cond1 = condition1;
		cond2 = condition2;
		this.colorManager = colorManager;
		Graph ret = initializeGraph(label);
		if (clusters != null) {
			for (ProteinCluster proteinCluster : clusters) {
				createNodesAndEdgesFromProteinCluster(proteinCluster, ret);
			}
		}
		scaleColors(ret);

		final File file = createFile(ret, outputFile);
		fixHeader(file, label);
		return file;
	}

	private void scaleColors(Graph graph) {
		List<Node> nodes = graph.getNode();
		Double max = -Double.MAX_VALUE;
		Double min = Double.MAX_VALUE;
		double minimumRatioForColor = ProteinClusterQuantParameters.getInstance().getMinimumRatioForColor();
		double maximumRatioForColor = ProteinClusterQuantParameters.getInstance().getMaximumRatioForColor();
		for (Node node2 : nodes) {
			try {
				double ratio = 0.0;
				boolean valid = false;
				for (edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Node.Att att : node2.getAtt()) {
					if (att.getName().equals("countRatio")) {
						ratio = Double.valueOf(att.getValue());
						valid = true;
						break;
					}
				}
				if (valid && !Double.isInfinite(ratio) && Double.compare(ratio, Double.MAX_VALUE) != 0
						&& Double.compare(ratio, -Double.MAX_VALUE) != 0) {
					if (ratio < min) {
						min = ratio;
					}
					if (ratio > max) {
						max = ratio;
					}
				}
			} catch (NumberFormatException e) {

			}
		}
		if (Double.compare(max, -Double.MAX_VALUE) == 0 || Double.compare(min, Double.MAX_VALUE) == 0) {
			return;
		}
		if (max > maximumRatioForColor) {
			max = maximumRatioForColor;
		}
		if (min < minimumRatioForColor) {
			min = minimumRatioForColor;
		}
		for (Node node2 : nodes) {
			try {

				double ratio = 0.0;
				boolean valid = false;
				for (edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Node.Att att : node2.getAtt()) {
					if (att.getName().equals("countRatio")) {
						ratio = Double.valueOf(att.getValue());
						valid = true;
					}
				}
				if (valid) {
					if (Double.compare(Double.POSITIVE_INFINITY, ratio) == 0) {
						ratio = max;
					} else if (Double.compare(Double.NEGATIVE_INFINITY, ratio) == 0) {
						ratio = min;
					} else if (ratio < minimumRatioForColor) {
						ratio = min;
					} else if (ratio > maximumRatioForColor) {
						ratio = max;
					}
					// this is a peptide
					Color color = ColorGenerator.getColor(ratio, min, max,
							ProteinClusterQuantParameters.getInstance().getColorRatioMin(),
							ProteinClusterQuantParameters.getInstance().getColorRatioMax());

					node2.getGraphics().setFill(ColorManager.getHexString(color));
				}

			} catch (NumberFormatException e) {

			}
		}

	}

	private Graph initializeGraph(String label) {
		Graph ret = factory.createGraph();
		ret.setId(Math.abs(new Long(System.currentTimeMillis()).intValue()));
		ret.setLabel(label);
		ret.getAtt().add(createAtt(PCQ_ID, label, STRING));
		// ret.getAtt().add(createAtt("name", label, STRING));
		ret.getAtt().add(createAtt("selected", "1", BOOLEAN));
		ret.getAtt().add(createAtt("layoutAlgorithm", "Prefuse Force Directed Layout", STRING));
		ret.setGraphics(createGeneralGraphics(label));
		return ret;
	}

	private void resetVisitedKeys() {
		visitedPeptideKeys.clear();
		visitedProteins.clear();
		edgesIDs.clear();
	}

	private void fixHeader(File file, String label) {
		String line;
		BufferedReader br = null;
		File outFile = null;
		BufferedWriter bw = null;
		try {
			outFile = File.createTempFile("PCQ_cytoscape", "xml");
			FileInputStream fis = new FileInputStream(file);
			InputStreamReader isr = new InputStreamReader(fis, Charset.forName("UTF-8"));
			br = new BufferedReader(isr);

			OutputStream fos = new FileOutputStream(outFile);
			OutputStreamWriter osr = new OutputStreamWriter(fos, Charset.forName("UTF-8"));
			bw = new BufferedWriter(osr);
			boolean sustitute = true;
			while ((line = br.readLine()) != null) {
				String newLine = line;
				if (sustitute && line.trim().startsWith("<graph")) {
					newLine = "<graph id=\"1\" label=\"" + label
							+ "\" directed=\"1\" cy:documentVersion=\"3.0\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:cy=\"http://www.cytoscape.org\" xmlns=\"http://www.cs.rpi.edu/XGMML\">";
					sustitute = false;
				}
				bw.write(newLine + "\n");
			}

		} catch (IOException e) {
		} finally {
			if (br != null) {
				try {
					br.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			if (bw != null) {
				try {
					bw.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			// copy the temp file to the orginial one
			if (!FileUtil.copyFile(outFile, file)) {
				throw new IllegalArgumentException("Error copying files");
			}
			// remove the temp file
			outFile.deleteOnExit();
		}

	}

	private void createNodesAndEdgesFromProteinCluster(ProteinCluster proteinCluster, Graph graph) {
		if (!proteinCluster.getProteinPairs().isEmpty()) {
			createNodesAndEdgesFromProteinPairs(proteinCluster.getProteinPairs(), graph);
			// now the proteins and peptides without protein pair
		} else {
			createNodesAndEdgesFromProteins(proteinCluster.getProteinSet().iterator().next(), null, false, false, false,
					null, null, graph);
		}
		// draw an edge between the aligned peptides
		final Set<QuantifiedPeptideInterface> peptideSet = proteinCluster.getPeptideSet();
		for (QuantifiedPeptideInterface quantifiedPeptide : peptideSet) {
			final Set<QuantifiedPeptideInterface> alignedPeptides = proteinCluster
					.getAlignedPeptides(quantifiedPeptide);
			for (QuantifiedPeptideInterface quantifiedPeptide2 : alignedPeptides) {

				Set<QuantifiedPeptideInterface> peptidesToLink = new HashSet<QuantifiedPeptideInterface>();
				peptidesToLink.add(quantifiedPeptide);
				peptidesToLink.add(quantifiedPeptide2);
				String edgeID = getUniqueID(peptidesToLink);
				if (!edgesIDs.containsKey(edgeID)) {

					String edgeName3 = Utils.getPeptidesSequenceString(peptidesToLink);
					final NWResult alignmentResult = proteinCluster.getAlignmentResult(quantifiedPeptide,
							quantifiedPeptide2);
					Map<String, AttributeValueType> attributes3 = getAttributesForEdge(edgeName3, null, peptidesToLink,
							alignmentResult, null, null);
					String tooltip = getHtml(getTooltipFromAlignment(alignmentResult));
					Edge edge = createEdge(++edgeCounter, null, tooltip, attributes3, getUniqueID(quantifiedPeptide),
							getUniqueID(quantifiedPeptide2), colorManager.getAlignedPeptidesEdgeColor());
					graph.getEdge().add(edge);
					edgesIDs.put(edgeID, edge);
				}

			}
		}
	}

	private String getTooltipFromAlignment(NWResult alignmentResult) {
		StringBuilder sb = new StringBuilder();
		if (alignmentResult != null) {
			sb.append("<b>Alignment Score=</b>" + alignmentResult.getFinalAlignmentScore());
			sb.append("\n<b>peptide 1:</b>" + alignmentResult.getSeq1());
			sb.append("\n<b>peptide 2:</b>" + alignmentResult.getSeq2());
			sb.append("\n<b>Lenth of alignment=</b>" + alignmentResult.getAlignmentLength());
			sb.append("\n<b>Identical segment length=</b>" + alignmentResult.getIdenticalLength());
			sb.append("\n<b>Identity=</b>" + formatter.format(alignmentResult.getSequenceIdentity() * 100) + "%");
			sb.append("\n<b>Max consecutive identity=</b>" + alignmentResult.getMaxConsecutiveIdenticalAlignment());
			sb.append("\n<b>Alignment string=</b>\n" + alignmentResult.getAlignmentString());
		}
		return sb.toString();
	}

	private void createNodesAndEdgesFromProteinPairs(Collection<ProteinPair> proteinPairs, Graph graph) {
		for (ProteinPair proteinPair : proteinPairs) {
			createNodesAndEdgesFromProteins(proteinPair.getProtein1(), proteinPair.getProtein2(),
					proteinPair.isSharedPeptidesInconsistent(), proteinPair.isUniquePeptidesProt1Inconsistent(),
					proteinPair.isUniquePeptidesProt2Inconsistent(), proteinPair.getClassification1Case().values(),
					proteinPair.getClassification2Cases().values(), graph);
		}
	}

	private void createNodesAndEdgesFromProteins(QuantifiedProteinInterface protein1,
			QuantifiedProteinInterface protein2, boolean isSharedPeptidesInconsistent,
			boolean isUniquePeptidesProt1Inconsistent, boolean isUniquePeptidesProt2Inconsistent,
			Collection<Classification1Case> classification1Cases, Collection<Classification2Case> classification2Cases,
			Graph graph) {

		// COLORS
		Color edgeColorU1 = Color.black;
		Color edgeColorS1 = Color.black;
		Color edgeColorS2 = Color.black;
		Color edgeColorU2 = Color.black;
		if (isSharedPeptidesInconsistent) {
			edgeColorS2 = colorManager.getHighlightColor();
			edgeColorS1 = colorManager.getHighlightColor();
		}
		if (isUniquePeptidesProt1Inconsistent) {
			edgeColorU1 = colorManager.getHighlightColor();
		}
		if (isUniquePeptidesProt2Inconsistent) {
			edgeColorU2 = colorManager.getHighlightColor();
		}

		Color outlineColorU1 = Color.black;
		Color outlineColorP1 = Color.black;
		Color outlineColorS12 = Color.black;
		Color outlineColorP2 = Color.black;
		Color outlineColorU2 = Color.black;
		if (isUniquePeptidesProt1Inconsistent) {
			outlineColorU1 = colorManager.getHighlightColor();
		}
		if (isUniquePeptidesProt2Inconsistent) {
			outlineColorU2 = colorManager.getHighlightColor();
		}
		if (isSharedPeptidesInconsistent) {
			outlineColorS12 = colorManager.getHighlightColor();
		}

		/////////////////////
		// NODES:
		// U1 peptides unique to protein1
		// P1 protein1
		// S12 peptides shared by protein 1 and protein 2
		// P2 protein2
		// U2 peptides unique to protein 2
		//
		// EDGES:
		// U1
		// S1
		// S2
		// U2
		//////////////////////

		/////////////////////
		// U1 peptides unique to protein1
		if (protein1 != null) {
			List<List<QuantifiedPeptideInterface>> pep1DoubleList = new ArrayList<List<QuantifiedPeptideInterface>>();
			final List<QuantifiedPeptideInterface> peptidesU1 = Utils.getUniquePeptides(protein1, protein2, true);
			if (ProteinClusterQuantParameters.getInstance().isCollapseIndistinguishablePeptides()) {
				pep1DoubleList.add(peptidesU1);
			} else {
				for (QuantifiedPeptideInterface quantifiedPeptide : peptidesU1) {
					List<QuantifiedPeptideInterface> list = new ArrayList<QuantifiedPeptideInterface>();
					list.add(quantifiedPeptide);
					pep1DoubleList.add(list);
				}
			}
			for (List<QuantifiedPeptideInterface> uniquePeptides_U1 : pep1DoubleList) {
				final String peptidesSequenceString_U1 = Utils.getPeptidesSequenceString(uniquePeptides_U1);
				if (!uniquePeptides_U1.isEmpty()) {
					if (!visitedPeptideKeys.containsKey(peptidesSequenceString_U1)) {
						final String nodeID = getUniqueID(uniquePeptides_U1);
						final String label = formatNumber(Utils.getPepRatio(uniquePeptides_U1, cond1, cond2, null)
								.getLog2CountRatio(cond1, cond2));
						Node node = createNodeFromPeptides(nodeID, label, getPeptidesTooltip(label, uniquePeptides_U1),
								uniquePeptides_U1, outlineColorU1);
						graph.getNode().add(node);
						visitedPeptideKeys.put(peptidesSequenceString_U1, node);
					}
					// set highLight color anyway
					if (outlineColorU1.equals(colorManager.getHighlightColor())) {
						setNodeOutlineColor(visitedPeptideKeys.get(peptidesSequenceString_U1),
								colorManager.getHighlightColor());
					}
				}
			}
			// P1 protein1
			if (!visitedProteins.containsKey(protein1)) {
				Node node = createNodeFromProtein(protein1, classification1Cases, classification2Cases, outlineColorP1);
				graph.getNode().add(node);
				visitedProteins.put(protein1, node);
			}
			// set highLight color anyway
			if (outlineColorP1.equals(colorManager.getHighlightColor())) {
				setNodeOutlineColor(visitedProteins.get(protein1), colorManager.getHighlightColor());
			}
			// EDGES:
			// U1
			for (List<QuantifiedPeptideInterface> uniquePeptides_U1 : pep1DoubleList) {
				String edgeID = getUniqueID(uniquePeptides_U1) + getUniqueID(protein1);
				String edgeID2 = getUniqueID(protein1) + getUniqueID(uniquePeptides_U1);
				if (!edgesIDs.containsKey(edgeID) && !edgesIDs.containsKey(edgeID2)) {
					final Ratio uniquePepRatio3 = Utils.getPepRatio(uniquePeptides_U1, cond1, cond2, null);
					String edgeName3 = getUniqueID(protein1) + "-" + Utils.getPeptidesSequenceString(uniquePeptides_U1);
					Map<String, AttributeValueType> attributes3 = getAttributesForEdge(edgeName3, uniquePepRatio3,
							uniquePeptides_U1, null, classification1Cases, classification2Cases);
					Edge edge = createEdge(++edgeCounter, null, null, attributes3, getUniqueID(uniquePeptides_U1),
							getUniqueID(protein1), edgeColorU1);
					graph.getEdge().add(edge);
					edgesIDs.put(edgeID, edge);
					edgesIDs.put(edgeID2, edge);
				}
				updateCasesForEdge(edgesIDs.get(edgeID), edgeColorU1, classification1Cases, classification2Cases);
			}

		}
		// S12 peptides shared by protein 1 and protein 2
		final Map<String, Set<QuantifiedPeptideInterface>> sharedPeptidesMap_S12 = Utils.getSharedPeptidesMap(protein1,
				protein2, false);

		for (Set<QuantifiedPeptideInterface> peptidesS12 : sharedPeptidesMap_S12.values()) {
			List<Set<QuantifiedPeptideInterface>> pep12DoubleList = new ArrayList<Set<QuantifiedPeptideInterface>>();
			if (ProteinClusterQuantParameters.getInstance().isCollapseIndistinguishablePeptides()) {
				pep12DoubleList.add(peptidesS12);
			} else {
				for (QuantifiedPeptideInterface peptide : peptidesS12) {
					Set<QuantifiedPeptideInterface> set = new HashSet<QuantifiedPeptideInterface>();
					set.add(peptide);
					pep12DoubleList.add(set);
				}
			}
			for (Set<QuantifiedPeptideInterface> sharedPeptides_S12 : pep12DoubleList) {
				final String sharedSequenceString_S12 = Utils.getPeptidesSequenceString(sharedPeptides_S12);
				if (!sharedPeptides_S12.isEmpty()) {
					if (!visitedPeptideKeys.containsKey(sharedSequenceString_S12)) {
						final String nodeID = getUniqueID(sharedPeptides_S12);
						final String label = formatNumber(Utils.getPepRatio(sharedPeptides_S12, cond1, cond2, null)
								.getLog2CountRatio(cond1, cond2));
						Node node = createNodeFromPeptides(nodeID, label, getPeptidesTooltip(label, sharedPeptides_S12),
								sharedPeptides_S12, outlineColorS12);
						graph.getNode().add(node);
						visitedPeptideKeys.put(sharedSequenceString_S12, node);
					}
					// set highLight color anyway
					if (outlineColorS12.equals(colorManager.getHighlightColor())) {
						setNodeOutlineColor(visitedPeptideKeys.get(sharedSequenceString_S12),
								colorManager.getHighlightColor());
					}

					// P1S12
					String edgeID = getUniqueID(protein1) + getUniqueID(sharedPeptides_S12);
					String edgeID2 = getUniqueID(sharedPeptides_S12) + getUniqueID(protein1);
					final Ratio sharedPepRatio = Utils.getSharedPepRatio(protein1, protein2, cond1, cond2);
					if (!edgesIDs.containsKey(edgeID) && !edgesIDs.containsKey(edgeID2)) {
						String edgeName = getUniqueID(protein1) + "-"
								+ Utils.getPeptidesSequenceString(sharedPeptides_S12);
						Map<String, AttributeValueType> attributes = getAttributesForEdge(edgeName, sharedPepRatio,
								sharedPeptides_S12, null, classification1Cases, classification2Cases);

						Edge edge = createEdge(++edgeCounter, null, null, attributes, getUniqueID(protein1),
								getUniqueID(sharedPeptides_S12), edgeColorS1);
						graph.getEdge().add(edge);
						edgesIDs.put(edgeID, edge);
						edgesIDs.put(edgeID2, edge);
					}
					updateCasesForEdge(edgesIDs.get(edgeID), edgeColorS1, classification1Cases, classification2Cases);

					if (protein2 != null) {
						// S12P2
						edgeID = getUniqueID(protein2) + getUniqueID(sharedPeptides_S12);
						edgeID2 = getUniqueID(sharedPeptides_S12) + getUniqueID(protein2);
						if (!edgesIDs.containsKey(edgeID) && !edgesIDs.containsKey(edgeID2)) {
							String edgeName2 = getUniqueID(protein2) + "-"
									+ Utils.getPeptidesSequenceString(sharedPeptides_S12);
							Map<String, AttributeValueType> attributes2 = getAttributesForEdge(edgeName2,
									sharedPepRatio, sharedPeptides_S12, null, classification1Cases,
									classification2Cases);
							Edge edge = createEdge(++edgeCounter, null, null, attributes2, getUniqueID(protein2),
									getUniqueID(sharedPeptides_S12), edgeColorS2);
							graph.getEdge().add(edge);
							edgesIDs.put(edgeID, edge);
							edgesIDs.put(edgeID2, edge);
						}
						updateCasesForEdge(edgesIDs.get(edgeID), edgeColorS2, classification1Cases,
								classification2Cases);

					}
				}
			}
		}
		// P2 protein2
		if (protein2 != null) {
			if (!visitedProteins.containsKey(protein2)) {
				Node node = createNodeFromProtein(protein2, classification1Cases, classification2Cases, outlineColorP2);
				graph.getNode().add(node);
				visitedProteins.put(protein2, node);
			}
			// set highLight color anyway
			if (outlineColorP2.equals(colorManager.getHighlightColor())) {
				setNodeOutlineColor(visitedProteins.get(protein2), colorManager.getHighlightColor());
			}

			// U2 peptides unique to protein 2
			List<List<QuantifiedPeptideInterface>> pep2DoubleList = new ArrayList<List<QuantifiedPeptideInterface>>();
			final List<QuantifiedPeptideInterface> peptidesU2 = Utils.getUniquePeptides(protein2, protein1, true);
			if (ProteinClusterQuantParameters.getInstance().isCollapseIndistinguishablePeptides()) {
				pep2DoubleList.add(peptidesU2);
			} else {
				for (QuantifiedPeptideInterface quantifiedPeptide : peptidesU2) {
					List<QuantifiedPeptideInterface> list = new ArrayList<QuantifiedPeptideInterface>();
					list.add(quantifiedPeptide);
					pep2DoubleList.add(list);
				}
			}
			for (List<QuantifiedPeptideInterface> uniquePeptides_U2 : pep2DoubleList) {
				final String peptidesSequenceString_U2 = Utils.getPeptidesSequenceString(uniquePeptides_U2);
				if (!uniquePeptides_U2.isEmpty()) {
					if (!visitedPeptideKeys.containsKey(peptidesSequenceString_U2)) {
						final String nodeID = getUniqueID(uniquePeptides_U2);
						final String label = formatNumber(Utils.getPepRatio(uniquePeptides_U2, cond1, cond2, null)
								.getLog2CountRatio(cond1, cond2));
						Node node = createNodeFromPeptides(nodeID, label, getPeptidesTooltip(label, uniquePeptides_U2),
								uniquePeptides_U2, outlineColorU2);
						graph.getNode().add(node);
						visitedPeptideKeys.put(peptidesSequenceString_U2, node);
					}
					// set highLight color anyway
					if (outlineColorU2.equals(colorManager.getHighlightColor())) {
						setNodeOutlineColor(visitedPeptideKeys.get(peptidesSequenceString_U2),
								colorManager.getHighlightColor());
					}
				}
				// EDGES
				// U2
				if (!uniquePeptides_U2.isEmpty()) {
					String edgeID = getUniqueID(protein2) + getUniqueID(uniquePeptides_U2);
					String edgeID2 = getUniqueID(uniquePeptides_U2) + getUniqueID(protein2);
					if (!edgesIDs.containsKey(edgeID) && !edgesIDs.containsKey(edgeID2)) {
						final Ratio uniquePepRatio4 = Utils.getPepRatio(uniquePeptides_U2, cond1, cond2, null);
						String edgeName4 = getUniqueID(protein2) + "-"
								+ Utils.getPeptidesSequenceString(uniquePeptides_U2);
						Map<String, AttributeValueType> attributes4 = getAttributesForEdge(edgeName4, uniquePepRatio4,
								uniquePeptides_U2, null, classification1Cases, classification2Cases);
						Edge edge = createEdge(++edgeCounter, null, null, attributes4, getUniqueID(protein2),
								getUniqueID(uniquePeptides_U2), edgeColorU2);
						graph.getEdge().add(edge);
						edgesIDs.put(edgeID, edge);
						edgesIDs.put(edgeID2, edge);
					}
					updateCasesForEdge(edgesIDs.get(edgeID), edgeColorU2, classification1Cases, classification2Cases);
				}

			}
		}

	}

	private void updateCasesForEdge(Edge edge, Color fillColor, Collection<Classification1Case> classification1Cases,
			Collection<Classification2Case> classification2Cases) {
		if (fillColor != null && fillColor.equals(colorManager.getHighlightColor())) {
			edge.getGraphics().setFill(ColorManager.getHexString(fillColor));
		}
		final List<edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Edge.Graphics.Att> atts = edge.getGraphics()
				.getAtt();
		String tooltip = null;
		Set<Classification1Case> cases1 = new HashSet<Classification1Case>();
		Set<Classification2Case> cases2 = new HashSet<Classification2Case>();
		for (edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Edge.Graphics.Att att : atts) {
			if (att.getName().equals("EDGE_LABEL")) {
				if (att.getValue() != null) {
					final String oldValueCases1 = att.getValue().split("\n")[0];
					if (oldValueCases1 != null) {
						if (oldValueCases1.contains(",")) {
							final String[] split = oldValueCases1.split(",");
							for (String string : split) {
								cases1.add(Classification1Case.getByCaseID(Integer.valueOf(string)));
							}
						} else {
							if (!oldValueCases1.equals(" ")) {
								cases1.add(Classification1Case.getByCaseID(Integer.valueOf(oldValueCases1)));
							}
						}
					}
					final String oldValueCases2 = att.getValue().split("\n")[1];
					if (oldValueCases2 != null) {
						if (oldValueCases2.contains(",")) {
							final String[] split = oldValueCases2.split(",");
							for (String string : split) {
								cases2.add(Classification2Case.getByCaseID(Integer.valueOf(string)));
							}
						} else {
							if (!oldValueCases2.equals(" ")) {
								cases2.add(Classification2Case.getByCaseID(Integer.valueOf(oldValueCases2)));
							}
						}
					}
				}
				if (classification1Cases != null)
					cases1.addAll(classification1Cases);
				if (classification2Cases != null)
					cases2.addAll(classification2Cases);
				tooltip = "Classification 1: " + getClassificationCase1String(cases1, true) + "\nClassification 2: "
						+ getClassificationCase2String(cases2, true);
				String case1String = " ";
				String case2String = " ";
				final String classificationCase1String = getClassificationCase1String(cases1, false);
				if (classificationCase1String != null) {
					case1String = classificationCase1String;
				}
				final String classificationCase2String = getClassificationCase2String(cases2, false);
				if (classificationCase2String != null) {
					case2String = classificationCase2String;
				}
				String value = case1String + "\n" + case2String;

				if (ProteinClusterQuantParameters.getInstance().isShowCasesInEdges()) {
					att.setValue(value);
				}
				break;
			}
		}
		if (tooltip != null) {
			for (edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Edge.Graphics.Att att : atts) {
				if (att.getName().equals("EDGE_TOOLTIP")) {
					final String html = getHtml("Outlier unique peptide(s)\n" + tooltip);
					att.setValue(html);
				}
			}
		}
	}

	private String getClassificationCase1String(Collection<Classification1Case> classification1Cases,
			boolean addExplanationOfCases) {
		StringBuilder sb = new StringBuilder();
		int index = 0;
		List<Classification1Case> list = new ArrayList<Classification1Case>();
		for (Classification1Case classification1Case : classification1Cases) {
			if (classification1Case == Classification1Case.CASE4 || classification1Case == Classification1Case.CASE5
					|| classification1Case == Classification1Case.CASE3) {
				continue;
			}
			if (!list.contains(classification1Case))
				list.add(classification1Case);
		}

		Collections.sort(list, new Comparator<Classification1Case>() {
			public int compare(Classification1Case o1, Classification1Case o2) {
				return Integer.compare(o1.getCaseID(), o2.getCaseID());
			}
		});
		for (Classification1Case case1 : list) {
			if (index > 0)
				sb.append(",");
			sb.append(case1.getCaseID());
			if (addExplanationOfCases) {
				sb.append(" (" + case1.getExplanation() + ") ");
			}
			index++;
		}
		if ("".equals(sb.toString()))
			return null;
		return sb.toString();
	}

	private String getClassificationCase2String(Collection<Classification2Case> classification2Cases,
			boolean addExplanationOfCases) {
		StringBuilder sb = new StringBuilder();
		int index = 0;
		List<Classification2Case> list = new ArrayList<Classification2Case>();
		for (Classification2Case classification2Case : classification2Cases) {
			if (classification2Case == Classification2Case.unclassified
					|| classification2Case == Classification2Case.no_difference) {
				continue;
			}
			if (!list.contains(classification2Case))
				list.add(classification2Case);
		}

		Collections.sort(list, new Comparator<Classification2Case>() {

			public int compare(Classification2Case o1, Classification2Case o2) {
				return Integer.compare(o1.getCaseID(), o2.getCaseID());
			}
		});
		for (Classification2Case case2 : list) {
			if (index > 0)
				sb.append(",");
			sb.append(case2.getCaseID());
			if (addExplanationOfCases) {
				sb.append(" (" + case2.getExplanation() + ") ");
			}
			index++;
		}
		if ("".equals(sb.toString()))
			return null;
		return sb.toString();
	}

	private void setNodeOutlineColor(Node node, Color outlineColor) {
		node.getGraphics().setOutline(ColorManager.getHexString(outlineColor));
	}

	/**
	 * Gets the tooltip that is shown when a peptide node is hovered by the
	 * mouse.<br>
	 * It is a list of the peptide sequences, together with the number of PSMs
	 * and ions counts pero each one.
	 *
	 * @param label
	 *
	 * @param peptides
	 * @param cond22
	 * @param cond1
	 * @return
	 */
	private String getPeptidesTooltip(String prefix, Collection<QuantifiedPeptideInterface> peptides) {
		StringBuilder sb = new StringBuilder();
		StringBuilder totalNumerator = new StringBuilder();
		StringBuilder totalDenominator = new StringBuilder();
		for (QuantifiedPeptideInterface peptide : Utils.getSortedPeptidesBySequence(peptides)) {
			if (peptide instanceof IsobaricQuantifiedPeptide) {
				IsobaricQuantifiedPeptide isoPeptide = (IsobaricQuantifiedPeptide) peptide;
				if (!"".equals(totalNumerator.toString())) {
					totalNumerator.append(" + ");
				}

				if (isoPeptide.getIonsByCondition().containsKey(cond1)) {
					totalNumerator.append(isoPeptide.getIonsByCondition().get(cond1).size() + "/"
							+ peptide.getQuantifiedPSMs().size());
				} else {
					totalNumerator.append("0");
				}
				if (!"".equals(totalDenominator.toString())) {
					totalDenominator.append(" + ");
				}
				if (isoPeptide.getIonsByCondition().containsKey(cond2)) {
					totalDenominator.append(isoPeptide.getIonsByCondition().get(cond2).size() + "/"
							+ peptide.getQuantifiedPSMs().size());
				} else {
					totalDenominator.append("0");
				}
				if (!"".equals(sb.toString()))
					sb.append("\n");
				final Map<QuantificationLabel, Set<Ion>> ionsByLabel = isoPeptide.getIonsByLabel();
				StringBuilder ionNumberString = new StringBuilder();
				if (ionsByLabel != null && !ionsByLabel.isEmpty()) {
					for (QuantificationLabel label : QuantificationLabel.values()) {
						if (ionsByLabel.containsKey(label)) {
							if (!"".equals(ionNumberString.toString())) {
								ionNumberString.append(", ");
							}
							ionNumberString.append(label.name() + ":" + ionsByLabel.get(label).size());
						}
					}
				}

				sb.append(peptide.getSequence() + "\t" + peptide.getQuantifiedPSMs().size() + " PSMs\t"
						+ ionNumberString.toString() + "\t Shared by " + peptide.getQuantifiedProteins().size()
						+ " proteins");
			}
		}
		return "Log2Ratio:" + prefix + " = log2( (" + totalNumerator + ") / (" + totalDenominator + ") )\n"
				+ sb.toString();
	}

	private String formatNumber(Double number) {
		if (Double.isInfinite(number)) {
			if (Double.compare(Double.POSITIVE_INFINITY, number) == 0) {
				return "INF";
			}
			if (Double.compare(Double.NEGATIVE_INFINITY, number) == 0) {
				return "-INF";
			}
		} else if (Double.isNaN(number)) {
			return "N/A";
		}
		return formatter.format(number);

	}

	private String getUniqueID(Collection<QuantifiedPeptideInterface> peptides) {
		return Utils.getPeptidesSequenceString(peptides);
	}

	private String getUniqueID(QuantifiedPeptideInterface peptide) {
		Set<QuantifiedPeptideInterface> set = new HashSet<QuantifiedPeptideInterface>();
		set.add(peptide);
		return Utils.getPeptidesSequenceString(set);
	}

	private String getUniqueID(QuantifiedProteinInterface protein) {
		if (ProteinClusterQuantParameters.getInstance().getProteinLabel() == ProteinLabel.ACC) {
			return protein.getAccession();
		} else {
			return getProteinNameString(protein.getAccession());
		}
	}

	private Node createNodeFromPeptides(String nodeID, String label, String tooltip,
			Collection<QuantifiedPeptideInterface> peptides, Color outlineColor) {
		final String sequenceString = Utils.getPeptidesSequenceString(peptides);
		Map<String, AttributeValueType> attributes = new HashMap<String, AttributeValueType>();
		attributes.put("PeptideSequences", new AttributeValueType(sequenceString));
		attributes.put(PCQ_ID, new AttributeValueType(sequenceString));
		attributes.put("numProteinsInPeptide",
				new AttributeValueType(peptides.iterator().next().getQuantifiedProteins().size()));
		attributes.put("numPsmsInPeptide", new AttributeValueType(Utils.getPSMs(peptides).size()));
		final Ratio pepRatio = Utils.getPepRatio(peptides, cond1, cond2, null);
		attributes.put("countRatio", new AttributeValueType(pepRatio.getLog2CountRatio(cond1, cond2)));
		attributes.put("ionCount", new AttributeValueType(Utils.getIonCount(peptides)));
		attributes.put("isProtein", new AttributeValueType(0));
		attributes.put("normLightIons", new AttributeValueType(pepRatio.getIonCount(cond1)));
		attributes.put("normHeavyIons", new AttributeValueType(pepRatio.getIonCount(cond2)));
		attributes.put("lightIons", new AttributeValueType(Utils.getIonCount(peptides, cond1)));
		attributes.put("heavyIons", new AttributeValueType(Utils.getIonCount(peptides, cond2)));
		BORDER_TYPE borderType = BORDER_TYPE.SOLID;
		Set<QuantifiedProteinInterface> proteins = new HashSet<QuantifiedProteinInterface>();
		for (QuantifiedPeptideInterface quantifiedPeptide : peptides) {
			proteins.addAll(quantifiedPeptide.getQuantifiedProteins());
		}
		String label_sufix = "";
		Color labelColor = Color.black;
		if (proteins.size() > 1) {
			// borderType = BORDER_TYPE.DASHED;
			// labelColor = colorManager.getSharedNodeLabelColor();
			attributes.put("uniquePep", new AttributeValueType(0));
		} else {
			attributes.put("uniquePep", new AttributeValueType(1));
			label_sufix = " *";
		}
		Node node = createNode(nodeID, sequenceString, attributes);
		node.setGraphics(createGraphicsNode(label + label_sufix, getHtml(tooltip),
				ProteinClusterQuantParameters.getInstance().getPeptideNodeShape(),
				ProteinClusterQuantParameters.getInstance().getPeptideNodeHeight(),
				ProteinClusterQuantParameters.getInstance().getPeptideNodeWidth(), outlineColor, Color.cyan, labelColor,
				borderType));
		return node;
	}

	private Node createNodeFromProtein(QuantifiedProteinInterface protein,
			Collection<Classification1Case> classification1Cases, Collection<Classification2Case> classification2Cases,
			Color outlineColor) {
		Map<String, AttributeValueType> attributes = new HashMap<String, AttributeValueType>();
		attributes.put("UniprotKB", new AttributeValueType(protein.getAccession()));
		attributes.put("isProtein", new AttributeValueType(1));
		attributes.put("numPeptidesInProteins", new AttributeValueType(protein.getQuantifiedPeptides().size()));
		attributes.put("numPsmsInProtein", new AttributeValueType(protein.getQuantifiedPSMs().size()));

		StringBuilder classification1String = new StringBuilder();
		StringBuilder classification1StringNOHTML = new StringBuilder();
		if (classification1Cases != null) {
			for (Classification1Case classification1Case : classification1Cases) {
				final String value = "<b>" + String.valueOf(classification1Case.getCaseID()) + "</b>: "
						+ classification1Case.getExplanation();
				if (!"".equals(classification1String.toString())) {
					classification1String.append(", ");
				}
				classification1String.append(value);
				classification1StringNOHTML.append(classification1Case.getCaseID());

			}
			attributes.put("Classification1Case",
					new AttributeValueType(classification1StringNOHTML.toString(), AttType.string));
		}

		StringBuilder classification2String = new StringBuilder();
		StringBuilder classification2StringNOHTML = new StringBuilder();
		if (classification2Cases != null) {
			for (Classification2Case classification2Case : classification2Cases) {
				final String value = "<b>" + String.valueOf(classification2Case.getCaseID()) + "</b>: "
						+ classification2Case.getExplanation();
				if (!"".equals(classification2String.toString())) {
					classification2String.append(", ");
				}
				classification2String.append(value);
				classification2StringNOHTML.append(classification2Case.getCaseID());

			}
			attributes.put("Classification2Case",
					new AttributeValueType(classification2StringNOHTML.toString(), AttType.string));
		}

		attributes.put("ProteinDescription", new AttributeValueType(protein.getDescription()
				.replace(Utils.PROTEIN_DESCRIPTION_SEPARATOR, " " + Utils.PROTEIN_DESCRIPTION_SEPARATOR + " ")));
		attributes.put("Species", new AttributeValueType(protein.getTaxonomy()));
		attributes.put("GeneName", new AttributeValueType(getGeneString(protein)));
		attributes.put("ID", new AttributeValueType(getProteinNameString(protein.getAccession())));
		attributes.put("ACC", new AttributeValueType(protein.getAccession()));
		BORDER_TYPE borderType = BORDER_TYPE.SOLID;
		String label_sufix = "";
		Color labelColor = Color.black;
		if (Utils.getUniquePeptides(protein).isEmpty()) {
			// borderType = BORDER_TYPE.DASHED;
			// labelColor = colorManager.getSharedNodeLabelColor();
			attributes.put("conclusive", new AttributeValueType("0"));
		} else {
			attributes.put("conclusive", new AttributeValueType("1"));
			label_sufix = " *";
		}

		Node node = createNode(getUniqueID(protein), getUniqueID(protein) + label_sufix, attributes);
		String tooltipText = "<b>Protein ACC(s):</b>\n" + protein.getAccession() + "\n<b>Protein name(s):</b>\n "
				+ protein.getDescription().replace(Utils.PROTEIN_DESCRIPTION_SEPARATOR, "\n");

		if (classification1Cases != null) {
			tooltipText += "\n<b>Classification1 case:</b> " + classification1String.toString();
		}
		if (classification2Cases != null) {
			tooltipText += "\n<b>Classification2 case:</b> " + classification2String.toString();
		}

		tooltipText += "\n<b>TAX:</b> " + protein.getTaxonomy() + "\n<b>Gene name:</b> " + getGeneString(protein);
		final String tooltip = getHtml(tooltipText);
		node.setGraphics(createGraphicsNode(node.getLabel(), tooltip,
				ProteinClusterQuantParameters.getInstance().getProteinNodeShape(),
				ProteinClusterQuantParameters.getInstance().getProteinNodeHeight(),
				ProteinClusterQuantParameters.getInstance().getProteinNodeWidth(), outlineColor,
				getFillColorByTaxonomy(protein.getTaxonomy()), labelColor, borderType));
		return node;
	}

	private String getGeneString(QuantifiedProteinInterface protein) {
		final String geneNameString = Utils.getGeneNameString(getAnnotatedProtein(protein.getAccession()), protein,
				null, false);
		return geneNameString;
	}

	private String getProteinNameString(String accString) {
		List<String> list = new ArrayList<String>();
		if (accString.contains(Utils.PROTEIN_ACC_SEPARATOR)) {
			final String[] split = accString.split(Utils.PROTEIN_ACC_SEPARATOR);
			for (String acc : split) {
				list.add(getProteinName(acc));
			}
		} else {
			list.add(getProteinName(accString));
		}
		StringBuilder sb = new StringBuilder();
		for (String id : list) {
			if (!"".equals(sb.toString())) {
				sb.append(Utils.PROTEIN_ACC_SEPARATOR);
			}
			sb.append(id);
		}
		return sb.toString();
	}

	private String getProteinName(String acc) {

		final Map<String, Entry> annotatedProteins = getAnnotatedProtein(acc);
		if (annotatedProteins.containsKey(acc)) {
			final String proteinName = annotatedProteins.get(acc).getName().get(0);
			if (proteinName.contains("obsolete")) {
				return acc;
			}
			return proteinName;
		}
		return null;

	}

	private String getHtml(String tooltipText) {
		String text = tooltipText.replace("\n", "<br>");
		String ret = "<html>" + text + "</html>";
		return ret;
	}

	private Color getFillColorByTaxonomy(String taxonomy) {
		if (taxonomy.contains("-")) {
			return colorManager.getMultiTaxonomyColor();
		}
		return colorManager.getColorForProteinTaxonomy(taxonomy);
	}

	private Color getOutlineColorByCase1(Classification1Case classificationCase) {
		return colorManager.getColorByClassification1Case(classificationCase);
	}

	private Color getOutlineColorByCase2(Classification2Case classificationCase) {
		return colorManager.getColorByClassification2Case(classificationCase);
	}

	private edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Node.Graphics createGraphicsNode(String label,
			String tooltip, Shape shape, int heigth, int width, Color outlineColor, Color fillColor, Color labelColor,
			BORDER_TYPE borderType) {
		edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Node.Graphics ret = factory.createGraphNodeGraphics();
		ret.setType(shape.name());
		ret.setH(heigth);
		ret.setW(width);
		ret.setWidth(4); // border width
		ret.setOutline(ColorManager.getHexString(outlineColor));
		ret.setFill(ColorManager.getHexString(fillColor));
		ret.getAtt().add(createNodeGraphicAtt("NODE_TOOLTIP", tooltip, STRING));
		ret.getAtt().add(createNodeGraphicAtt("NODE_NESTED_NETWORK_IMAGE_VISIBLE", "true", STRING));
		ret.getAtt().add(createNodeGraphicAtt("NODE_BORDER_STROKE", borderType.name(), STRING));
		ret.getAtt().add(createNodeGraphicAtt("NODE_SELECTED", "false", STRING));
		ret.getAtt().add(createNodeGraphicAtt("NODE_TRANSPARENCY", "255", STRING));
		// TODO maybe change the label width?
		ret.getAtt().add(createNodeGraphicAtt("NODE_LABEL_WIDTH", "200", STRING));
		ret.getAtt().add(createNodeGraphicAtt("NODE_LABEL", label, STRING));
		ret.getAtt().add(createNodeGraphicAtt("NODE_LABEL_FONT_SIZE", "12", STRING));
		ret.getAtt().add(createNodeGraphicAtt("NODE_LABEL_TRANSPARENCY", "255", STRING));
		ret.getAtt().add(createNodeGraphicAtt("NODE_LABEL_COLOR", ColorManager.getHexString(labelColor), STRING));
		ret.getAtt().add(createNodeGraphicAtt("NODE_VISIBLE", "true", STRING));
		ret.getAtt().add(createNodeGraphicAtt("NODE_DEPTH", "0.0", STRING));
		ret.getAtt().add(createNodeGraphicAtt("NODE_BORDER_TRANSPARENCY", "255", STRING));
		ret.getAtt().add(createNodeGraphicAtt("NODE_LABEL_FONT_FACE", "Dialog,plain,12", STRING));
		ret.getAtt().add(createNodeGraphicAtt("NODE_LABEL_TRANSPARENCY", "255", STRING));

		return ret;
	}

	private Map<String, AttributeValueType> getAttributesForEdge(String edgeName, Ratio sharedPepRatio,
			Collection<QuantifiedPeptideInterface> sharedPeptides, NWResult alignmentResult,
			Collection<Classification1Case> classification1Cases,
			Collection<Classification2Case> classification2Cases) {
		Map<String, AttributeValueType> ret = new HashMap<String, AttributeValueType>();
		// ret.put("name", edgeName);
		ret.put(PCQ_ID, new AttributeValueType(edgeName));
		if (sharedPepRatio != null) {
			final Double log2CountRatio = sharedPepRatio.getLog2CountRatio(cond1, cond2);
			if (log2CountRatio != null) {
				ret.put("countRatio", new AttributeValueType(log2CountRatio));
			}
		}

		StringBuilder classification1String = new StringBuilder();
		if (classification1Cases != null) {
			for (Classification1Case classification1Case : classification1Cases) {
				final String value = String.valueOf(classification1Case.getCaseID());
				if (!"".equals(classification1String.toString())) {
					classification1String.append(", ");
				}
				classification1String.append(value);
			}
			ret.put("Classification1Case", new AttributeValueType(classification1String.toString(), AttType.string));
		}

		StringBuilder classification2String = new StringBuilder();
		if (classification2Cases != null) {
			for (Classification2Case classification2Case : classification2Cases) {
				final String value = String.valueOf(classification2Case.getCaseID());
				if (!"".equals(classification2String.toString())) {
					classification2String.append(", ");
				}
				classification2String.append(value);
			}
			ret.put("Classification2Case", new AttributeValueType(classification2String.toString(), AttType.string));
		}

		if (alignmentResult != null) {
			ret.put("Alignment score", new AttributeValueType(alignmentResult.getFinalAlignmentScore()));
			ret.put("Alignment length", new AttributeValueType(alignmentResult.getAlignmentLength()));
			ret.put("Alignment identity", new AttributeValueType(alignmentResult.getSequenceIdentity()));
			ret.put("Alignment identical length", new AttributeValueType(alignmentResult.getIdenticalLength()));
			ret.put("Alignment segment maximum length",
					new AttributeValueType(alignmentResult.getMaxConsecutiveIdenticalAlignment()));
			ret.put("Homology connection", new AttributeValueType("true"));
		}
		return ret;
	}

	private Edge createEdge(int edgeID, String label, String tooltip, Map<String, AttributeValueType> attributes3,
			String sourceID, String targetID, Color edgeColor) {
		final Edge edge = factory.createGraphEdge();
		edge.setId(edgeID);
		edge.setLabel(label);
		edge.setSource(sourceID);
		edge.setTarget(targetID);
		for (String attName : attributes3.keySet()) {
			final AttributeValueType attValue = attributes3.get(attName);
			edge.getAtt().add(createEdgeAttribute(attName, attValue.getValue().toString(), attValue.getType().name()));
		}
		edge.setGraphics(createGraphicsEdge(label, tooltip, edgeColor));
		return edge;
	}

	private edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Edge.Graphics createGraphicsEdge(String label,
			String tooltip, Color fillColor) {
		final edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Edge.Graphics edge = factory.createGraphEdgeGraphics();
		edge.setFill(ColorManager.getHexString(fillColor));
		edge.setWidth(3);
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_SELECTED", "false", STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_LABEL_COLOR", ColorManager.getHexString(Color.black), STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_TARGET_ARROW_SHAPE", "none", STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_SOURCE_ARROW_UNSELECTED_PAINT",
				ColorManager.getHexString(Color.black), STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_TARGET_ARROW_SELECTED_PAINT", "#FFFF00", STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_LABEL_TRANSPARENCY", "255", STRING));
		edge.getAtt().add(
				createEdgeGraphAttribute("EDGE_STROKE_SELECTED_PAINT", ColorManager.getHexString(Color.red), STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_LINE_TYPE", "SOLID", STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_TOOLTIP", tooltip, STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_CURVED", "true", STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_TRANSPARENCY", "255", STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_BEND", "", STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_LABEL_FONT_FACE", "Dialog,plain,10", STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_LABEL", label, STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_SOURCE_ARROW_SHAPE", "none", STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_TARGET_ARROW_UNSELECTED_PAINT",
				ColorManager.getHexString(Color.black), STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_SOURCE_ARROW_SELECTED_PAINT", "#FFFF00", STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_VISIBLE", "true", STRING));
		return edge;
	}

	private edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Edge.Att createEdgeAttribute(String name, String value,
			String type) {
		final edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Edge.Att ret = factory.createGraphEdgeAtt();
		ret.setName(name);
		ret.setValue(value);
		ret.setType(type);
		return ret;
	}

	private edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Edge.Graphics.Att createEdgeGraphAttribute(String name,
			String value, String type) {
		edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Edge.Graphics.Att ret = factory.createGraphEdgeGraphicsAtt();
		ret.setName(name);
		ret.setValue(value);
		ret.setType(type);
		return ret;
	}

	private Node createNode(String nodeID, String label, Map<String, AttributeValueType> attributes) {
		Node node = factory.createGraphNode();
		node.setId(nodeID);
		node.setLabel(label);
		for (String attName : attributes.keySet()) {
			AttributeValueType attValue = attributes.get(attName);

			node.getAtt().add(createNodeAtt(attName, attValue.getValue().toString(), attValue.getType().name()));

		}
		return node;
	}

	private edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Node.Att createNodeAtt(String name, String value,
			String type) {
		final edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Node.Att ret = factory.createGraphNodeAtt();
		ret.setName(name);
		ret.setValue(value);
		ret.setType(type);
		return ret;

	}

	private edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Node.Graphics.Att createNodeGraphicAtt(String name,
			String value, String type) {
		final edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Node.Graphics.Att ret = factory
				.createGraphNodeGraphicsAtt();
		ret.setName(name);
		ret.setValue(value);
		ret.setType(type);
		return ret;
	}

	private Graphics createGeneralGraphics(String label) {
		final Graphics ret = factory.createGraphGraphics();
		ret.getAtt().add(createGraphicAtt("NETWORK_NODE_SELECTION", "true", STRING));
		ret.getAtt().add(createGraphicAtt("NETWORK_HEIGHT", "381.0", STRING));
		ret.getAtt().add(createGraphicAtt("NETWORK_TITLE", label, STRING));
		ret.getAtt().add(createGraphicAtt("NETWORK_EDGE_SELECTION", "true", STRING));
		ret.getAtt().add(createGraphicAtt("NETWORK_SCALE_FACTOR", "0.23", STRING));
		ret.getAtt().add(createGraphicAtt("NETWORK_WIDTH", "946.0", STRING));
		ret.getAtt().add(createGraphicAtt("NETWORK_DEPTH", "0.0", STRING));
		ret.getAtt().add(createGraphicAtt("NETWORK_BACKGROUND_PAINT", "#FFFFFF", STRING));
		return ret;
	}

	private edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Graphics.Att createGraphicAtt(String name, String value,
			String type) {
		final edu.scripps.yates.proteinclusters.xgmml.jaxb.Graph.Graphics.Att ret = factory.createGraphGraphicsAtt();
		ret.setName(name);
		ret.setValue(value);
		ret.setType(type);
		return ret;
	}

	private Att createAtt(String name, String value, String type) {
		final Att ret = factory.createGraphAtt();
		ret.setName(name);
		ret.setValue(value);
		ret.setType(type);
		return ret;
	}

	private File createFile(Graph graph, File outputFile) throws JAXBException {
		final Marshaller marshaller = getJAXBContext().createMarshaller();

		marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, new Boolean(true));
		try {
			marshaller.setProperty("com.sun.xml.bind.indentString", "\t");
		} catch (PropertyException e) {
			marshaller.setProperty("com.sun.xml.internal.bind.indentString", "\t");
		}
		marshaller.marshal(graph, outputFile);
		log.debug(outputFile.getAbsolutePath() + " created");
		return outputFile;
	}

	private Map<String, Entry> getAnnotatedProtein(String accession) {
		Set<String> accs = new HashSet<String>();
		if (accession.contains(" ")) {
			final String[] split = accession.split("\\s+");
			for (String string : split) {
				accs.add(string);
			}
		} else {
			accs.add(accession);
		}
		Map<String, Entry> ret = new HashMap<String, Entry>();
		for (String string : accs) {
			if (annotatedProteins.containsKey(string)) {
				ret.put(string, annotatedProteins.get(string));
			}
		}

		return ret;
	}

	/**
	 * @return the annotatedProteins
	 */
	public Map<String, Entry> getAnnotatedProteins() {
		return annotatedProteins;
	}

	/**
	 * @param annotatedProteins
	 *            the annotatedProteins to set
	 */
	public void setAnnotatedProteins(Map<String, Entry> annotatedProteins) {
		this.annotatedProteins = annotatedProteins;
	}
}
