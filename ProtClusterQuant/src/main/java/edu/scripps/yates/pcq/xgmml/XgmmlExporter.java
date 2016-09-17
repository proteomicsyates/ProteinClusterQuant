package edu.scripps.yates.pcq.xgmml;

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

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.xml.Entry;
import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.IsobaricQuantifiedPeptide;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.pcq.ProteinClusterQuantParameters;
import edu.scripps.yates.pcq.cases.Classification1Case;
import edu.scripps.yates.pcq.cases.Classification2Case;
import edu.scripps.yates.pcq.model.IonCountRatio;
import edu.scripps.yates.pcq.model.PCQPeptideNode;
import edu.scripps.yates.pcq.model.PCQProteinNode;
import edu.scripps.yates.pcq.model.ProteinCluster;
import edu.scripps.yates.pcq.model.ProteinPair;
import edu.scripps.yates.pcq.util.PCQUtils;
import edu.scripps.yates.pcq.xgmml.jaxb.Graph;
import edu.scripps.yates.pcq.xgmml.jaxb.Graph.Att;
import edu.scripps.yates.pcq.xgmml.jaxb.Graph.Edge;
import edu.scripps.yates.pcq.xgmml.jaxb.Graph.Graphics;
import edu.scripps.yates.pcq.xgmml.jaxb.Graph.Node;
import edu.scripps.yates.pcq.xgmml.jaxb.ObjectFactory;
import edu.scripps.yates.pcq.xgmml.util.ColorManager;
import edu.scripps.yates.pcq.xgmml.util.ProteinNodeLabel;
import edu.scripps.yates.pcq.xgmml.util.Shape;
import edu.scripps.yates.utilities.alignment.nwalign.NWResult;
import edu.scripps.yates.utilities.colors.ColorGenerator;
import edu.scripps.yates.utilities.proteomicsmodel.Score;

public class XgmmlExporter {
	private final static Logger log = Logger.getLogger(XgmmlExporter.class);
	private static JAXBContext context;
	private QuantCondition cond1;
	private QuantCondition cond2;
	private int edgeCounter = 0;
	private final static String STRING = "string";
	private final static String BOOLEAN = "boolean";
	// private final static String LIST = "list";
	private final static ObjectFactory factory = new ObjectFactory();
	private final static DecimalFormat formatter = new DecimalFormat("#.#");
	private final static DecimalFormat formatter3Decimals = new DecimalFormat("#.###");

	private ColorManager colorManager;
	private final Map<String, Node> visitedPeptideKeys = new HashMap<String, Node>();
	private final Map<PCQProteinNode, Node> visitedProteinNodes = new HashMap<PCQProteinNode, Node>();
	private final Map<String, Edge> edgesIDs = new HashMap<String, Edge>();
	private Map<String, Entry> annotatedProteins;
	private static final String PCQ_ID = "PCQ_ID";
	private static final String COUNT_RATIO = "countRatio";
	private static final String FINAL_RATIO = "finalRatio";
	private static int numProteinNodes = 1;
	private static final Map<String, String> numProteinByProteinNodeID = new HashMap<String, String>();
	private static final Set<String> ratioAttributes = new HashSet<String>();
	private static final String WEIGHT = "Weight";
	private static final String VARIANCE = "Variance";
	private static final String SIGNIFICANTLY_REGULATED_ATTRIBUTE = "significant";
	private static final String IS_FILTERED = "isFiltered";
	static {
		ratioAttributes.add(COUNT_RATIO);
		ratioAttributes.add(FINAL_RATIO);
	}

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

	protected File exportToGmmlFromProteinClustersUsingNodes(File outputFile, String label,
			Collection<ProteinCluster> clusters, QuantCondition condition1, QuantCondition condition2,
			ColorManager colorManager) throws JAXBException {

		resetVisitedKeys();
		cond1 = condition1;
		cond2 = condition2;
		this.colorManager = colorManager;
		Graph ret = initializeGraph(label);
		if (clusters != null) {
			for (ProteinCluster proteinCluster : clusters) {
				createNodesAndEdgesFromProteinClusterUsingNodes(proteinCluster, ret);
			}
		}
		scaleColors(ret);

		final File file = createFile(ret, outputFile);
		fixHeader(file, label);
		return file;
	}

	/**
	 * Scale all colors of the peptide nodes according to input parameters
	 * settings. See parameters: minimumRatioForColor, maximumRatioForColor,
	 *
	 * @param graph
	 */
	private void scaleColors(Graph graph) {

		List<Node> nodes = graph.getNode();
		Double max = -Double.MAX_VALUE;
		Double min = Double.MAX_VALUE;
		final ProteinClusterQuantParameters params = ProteinClusterQuantParameters.getInstance();
		double minimumRatioForColor = params.getMinimumRatioForColor();
		double maximumRatioForColor = params.getMaximumRatioForColor();

		for (Node node2 : nodes) {
			try {
				double ratio = 0.0;
				boolean valid = false;

				for (edu.scripps.yates.pcq.xgmml.jaxb.Graph.Node.Att att : node2.getAtt()) {
					for (String ratioAttribute : ratioAttributes) {
						if (att.getName().equals(ratioAttribute)) {
							ratio = Double.valueOf(att.getValue());
							valid = true;
							break;
						}
					}
					if (valid) {
						break;
					}

				}
				if (valid && !Double.isNaN(ratio) && !Double.isInfinite(ratio)
						&& Double.compare(ratio, Double.MAX_VALUE) != 0
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
		if (min > max) {
			min = max;
		}
		for (Node node2 : nodes) {
			try {
				boolean isFiltered = false;
				double ratio = 0.0;
				boolean valid = false;
				boolean significantlyRegulated = false;
				for (edu.scripps.yates.pcq.xgmml.jaxb.Graph.Node.Att att : node2.getAtt()) {
					for (String ratioAttribute : ratioAttributes) {
						if (att.getName().equals(ratioAttribute)) {
							ratio = Double.valueOf(att.getValue());
							valid = true;
							break;
						}
					}
					if (att.getName().equals(SIGNIFICANTLY_REGULATED_ATTRIBUTE)) {
						if (att.getValue().equals("1")) {
							significantlyRegulated = true;
						}
					}
					if (att.getName().equals(IS_FILTERED)) {
						if (att.getValue().equals("1")) {
							isFiltered = true;
						}
					}

					if (valid) {
						break;
					}
				}
				if (valid) {
					if (isFiltered) {
						continue;
					}
					if (!significantlyRegulated && params.getColorNonRegulated() != null) {
						Color color = params.getColorNonRegulated();
						node2.getGraphics().setFill(ColorManager.getHexString(color));
					} else {
						if (Double.compare(Double.POSITIVE_INFINITY, ratio) == 0) {
							ratio = max;
						} else if (Double.compare(Double.NEGATIVE_INFINITY, ratio) == 0) {
							ratio = min;
						} else if (ratio < minimumRatioForColor) {
							ratio = min;
						} else if (ratio > maximumRatioForColor) {
							ratio = max;
						} else if (Double.isNaN(ratio)) {
							// skip it
							continue;
						}
						// this is a peptide
						Color color = ColorGenerator.getColor(ratio, min, max, params.getColorRatioMin(),
								params.getColorRatioMax());

						node2.getGraphics().setFill(ColorManager.getHexString(color));
					}
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
		visitedProteinNodes.clear();
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
			try {
				org.apache.commons.io.FileUtils.copyFile(outFile, file);
			} catch (IOException e) {
				throw new IllegalArgumentException(e);
			}

			// remove the temp file
			outFile.deleteOnExit();
		}

	}

	private void createNodesAndEdgesFromProteinClusterUsingNodes(ProteinCluster proteinCluster, Graph graph) {
		final Set<PCQProteinNode> proteinNodes = proteinCluster.getProteinNodes();
		if (!proteinCluster.getProteinPairs().isEmpty()) {
			createNodesAndEdgesFromProteinPairsUsingNodes(proteinCluster.getProteinPairs(), graph);
			// now the proteins and peptides without protein pair because they
			// were discarded probably
			for (PCQProteinNode pcqProteinNode : proteinNodes) {
				if (pcqProteinNode.getProteinPair() == null) {
					createNodesAndEdgesFromProteinNodes(pcqProteinNode, null, false, false, false, null, null, graph);
				}
			}

		} else {
			createNodesAndEdgesFromProteinNodes(proteinNodes.iterator().next(), null, false, false, false, null, null,
					graph);
		}
		// draw an edge between the aligned peptides
		final Set<PCQPeptideNode> peptideNodes = proteinCluster.getPeptideNodes();
		for (PCQPeptideNode peptideNode : peptideNodes) {
			final Set<PCQPeptideNode> alignedPeptides = proteinCluster.getAlignedPeptides(peptideNode);
			for (PCQPeptideNode peptideNode2 : alignedPeptides) {

				Set<PCQPeptideNode> peptidesToLink = new HashSet<PCQPeptideNode>();
				peptidesToLink.add(peptideNode);
				peptidesToLink.add(peptideNode2);
				String edgeID = getUniqueID(peptidesToLink);
				if (!edgesIDs.containsKey(edgeID)) {

					String edgeName3 = PCQUtils.getPeptideNodesSequenceString(peptidesToLink);
					final NWResult alignmentResult = proteinCluster.getAlignmentResult(peptideNode, peptideNode2);
					Map<String, AttributeValueType> attributes3 = getAttributesForEdge(edgeName3, null, alignmentResult,
							null, null);
					String tooltip = getHtml(getTooltipFromAlignment(alignmentResult));
					Edge edge = createEdge(++edgeCounter, null, tooltip, attributes3, peptideNode.getKey(),
							peptideNode2.getKey(), colorManager.getAlignedPeptidesEdgeColor());
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

	private void createNodesAndEdgesFromProteinPairsUsingNodes(Collection<ProteinPair> proteinPairs, Graph graph) {
		for (ProteinPair proteinPair : proteinPairs) {
			createNodesAndEdgesFromProteinNodes(proteinPair.getProteinNode1(), proteinPair.getProteinNode2(),
					proteinPair.isSharedPeptidesInconsistent(), proteinPair.isUniquePeptidesProt1Inconsistent(),
					proteinPair.isUniquePeptidesProt2Inconsistent(), proteinPair.getClassification1Case().values(),
					proteinPair.getClassification2Cases().values(), graph);
		}
	}

	private void createNodesAndEdgesFromProteinNodes(PCQProteinNode proteinNode1, PCQProteinNode proteinNode2,
			boolean isSharedPeptidesInconsistent, boolean isUniquePeptidesProt1Inconsistent,
			boolean isUniquePeptidesProt2Inconsistent, Collection<Classification1Case> classification1Cases,
			Collection<Classification2Case> classification2Cases, Graph graph) {

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
		if (proteinNode1 != null) {
			final Set<PCQPeptideNode> peptidesNodesU1 = PCQUtils.getUniquePeptideNodes(proteinNode1, proteinNode2, true,
					ProteinClusterQuantParameters.getInstance().isRemoveFilteredNodes());

			for (PCQPeptideNode uniquePeptideNode_U1 : peptidesNodesU1) {
				final String peptidesSequenceString_U1 = uniquePeptideNode_U1.getKey();

				if (!visitedPeptideKeys.containsKey(peptidesSequenceString_U1)) {
					final String nodeID = uniquePeptideNode_U1.getKey();
					final String label = formatNumber(
							PCQUtils.getRatioValue(uniquePeptideNode_U1.getConsensusRatio(cond1, cond2), cond1, cond2));
					Node node = createNodeFromPeptideNode(nodeID, label,
							getPeptideNodeTooltip(label, uniquePeptideNode_U1), uniquePeptideNode_U1, outlineColorU1);
					graph.getNode().add(node);
					visitedPeptideKeys.put(peptidesSequenceString_U1, node);
				}
				// set highLight color anyway
				if (outlineColorU1.equals(colorManager.getHighlightColor())) {
					setNodeOutlineColor(visitedPeptideKeys.get(peptidesSequenceString_U1),
							colorManager.getHighlightColor());
				}

			}
			// P1 protein1
			if (!visitedProteinNodes.containsKey(proteinNode1)) {
				Node node = createNodeFromProteinNode(proteinNode1, classification1Cases, classification2Cases,
						outlineColorP1);
				graph.getNode().add(node);
				visitedProteinNodes.put(proteinNode1, node);
			}
			// set highLight color anyway
			if (outlineColorP1.equals(colorManager.getHighlightColor())) {
				setNodeOutlineColor(visitedProteinNodes.get(proteinNode1), colorManager.getHighlightColor());
			}
			// EDGES:
			// U1
			for (PCQPeptideNode uniquePeptideNode_U1 : peptidesNodesU1) {
				String edgeID = uniquePeptideNode_U1.getKey() + getUniqueID(proteinNode1);
				String edgeID2 = getUniqueID(proteinNode1) + uniquePeptideNode_U1.getKey();
				if (!edgesIDs.containsKey(edgeID) && !edgesIDs.containsKey(edgeID2)) {
					final QuantRatio uniquePepRatio3 = uniquePeptideNode_U1.getConsensusRatio(cond1, cond2);
					String edgeName3 = getUniqueID(proteinNode1) + "-" + uniquePeptideNode_U1.getKey();
					Map<String, AttributeValueType> attributes3 = getAttributesForEdge(edgeName3, uniquePepRatio3, null,
							classification1Cases, classification2Cases);
					Edge edge = createEdge(++edgeCounter, null, null, attributes3, uniquePeptideNode_U1.getKey(),
							getUniqueID(proteinNode1), edgeColorU1);
					graph.getEdge().add(edge);
					edgesIDs.put(edgeID, edge);
					edgesIDs.put(edgeID2, edge);
				}
				updateCasesForEdge(edgesIDs.get(edgeID), edgeColorU1, classification1Cases, classification2Cases);
			}

		}
		// S12 peptides shared by protein 1 and protein 2
		final Map<String, Set<PCQPeptideNode>> sharedPeptideNodesMap_S12 = PCQUtils.getSharedPeptideNodesMap(
				proteinNode1, proteinNode2, false, ProteinClusterQuantParameters.getInstance().isRemoveFilteredNodes());

		for (Set<PCQPeptideNode> peptideNodesS12 : sharedPeptideNodesMap_S12.values()) {

			for (PCQPeptideNode sharedPeptideNode_S12 : peptideNodesS12) {
				final String sharedSequenceString_S12 = sharedPeptideNode_S12.getKey();

				if (!visitedPeptideKeys.containsKey(sharedSequenceString_S12)) {
					final String nodeID = sharedPeptideNode_S12.getKey();
					final String label = formatNumber(PCQUtils
							.getRatioValue(sharedPeptideNode_S12.getConsensusRatio(cond1, cond2), cond1, cond2));
					Node node = createNodeFromPeptideNode(nodeID, label,
							getPeptideNodeTooltip(label, sharedPeptideNode_S12), sharedPeptideNode_S12,
							outlineColorS12);
					graph.getNode().add(node);
					visitedPeptideKeys.put(sharedSequenceString_S12, node);
				}
				// set highLight color anyway
				if (outlineColorS12.equals(colorManager.getHighlightColor())) {
					setNodeOutlineColor(visitedPeptideKeys.get(sharedSequenceString_S12),
							colorManager.getHighlightColor());
				}

				// P1S12
				String edgeID = getUniqueID(proteinNode1) + sharedPeptideNode_S12.getKey();
				String edgeID2 = sharedPeptideNode_S12.getKey() + getUniqueID(proteinNode1);
				final QuantRatio sharedPepRatio = sharedPeptideNode_S12.getConsensusRatio(cond1, cond2);
				if (!edgesIDs.containsKey(edgeID) && !edgesIDs.containsKey(edgeID2)) {
					String edgeName = getUniqueID(proteinNode1) + "-" + sharedPeptideNode_S12.getKey();
					Map<String, AttributeValueType> attributes = getAttributesForEdge(edgeName, sharedPepRatio, null,
							classification1Cases, classification2Cases);

					Edge edge = createEdge(++edgeCounter, null, null, attributes, getUniqueID(proteinNode1),
							sharedPeptideNode_S12.getKey(), edgeColorS1);
					graph.getEdge().add(edge);
					edgesIDs.put(edgeID, edge);
					edgesIDs.put(edgeID2, edge);
				}
				updateCasesForEdge(edgesIDs.get(edgeID), edgeColorS1, classification1Cases, classification2Cases);

				if (proteinNode2 != null) {
					// S12P2
					edgeID = getUniqueID(proteinNode2) + sharedPeptideNode_S12.getKey();
					edgeID2 = sharedPeptideNode_S12.getKey() + getUniqueID(proteinNode2);
					if (!edgesIDs.containsKey(edgeID) && !edgesIDs.containsKey(edgeID2)) {
						String edgeName2 = getUniqueID(proteinNode2) + "-" + sharedPeptideNode_S12.getKey();
						Map<String, AttributeValueType> attributes2 = getAttributesForEdge(edgeName2, sharedPepRatio,
								null, classification1Cases, classification2Cases);
						Edge edge = createEdge(++edgeCounter, null, null, attributes2, getUniqueID(proteinNode2),
								sharedPeptideNode_S12.getKey(), edgeColorS2);
						graph.getEdge().add(edge);
						edgesIDs.put(edgeID, edge);
						edgesIDs.put(edgeID2, edge);
					}
					updateCasesForEdge(edgesIDs.get(edgeID), edgeColorS2, classification1Cases, classification2Cases);

				}

			}
		}
		// P2 protein2
		if (proteinNode2 != null) {
			if (!visitedProteinNodes.containsKey(proteinNode2)) {
				Node node = createNodeFromProteinNode(proteinNode2, classification1Cases, classification2Cases,
						outlineColorP2);
				graph.getNode().add(node);
				visitedProteinNodes.put(proteinNode2, node);
			}
			// set highLight color anyway
			if (outlineColorP2.equals(colorManager.getHighlightColor())) {
				setNodeOutlineColor(visitedProteinNodes.get(proteinNode2), colorManager.getHighlightColor());
			}

			// U2 peptides unique to protein 2
			final Set<PCQPeptideNode> peptidesU2 = PCQUtils.getUniquePeptideNodes(proteinNode2, proteinNode1, true,
					ProteinClusterQuantParameters.getInstance().isRemoveFilteredNodes());
			for (PCQPeptideNode uniquePeptides_U2 : peptidesU2) {
				final String peptidesSequenceString_U2 = uniquePeptides_U2.getKey();

				if (!visitedPeptideKeys.containsKey(peptidesSequenceString_U2)) {
					final String nodeID = uniquePeptides_U2.getKey();
					final String label = formatNumber(
							PCQUtils.getRatioValue(uniquePeptides_U2.getConsensusRatio(cond1, cond2), cond1, cond2));
					Node node = createNodeFromPeptideNode(nodeID, label,
							getPeptideNodeTooltip(label, uniquePeptides_U2), uniquePeptides_U2, outlineColorU2);
					graph.getNode().add(node);
					visitedPeptideKeys.put(peptidesSequenceString_U2, node);
				}
				// set highLight color anyway
				if (outlineColorU2.equals(colorManager.getHighlightColor())) {
					setNodeOutlineColor(visitedPeptideKeys.get(peptidesSequenceString_U2),
							colorManager.getHighlightColor());
				}

				// EDGES
				// U2

				String edgeID = getUniqueID(proteinNode2) + uniquePeptides_U2.getKey();
				String edgeID2 = uniquePeptides_U2.getKey() + getUniqueID(proteinNode2);
				if (!edgesIDs.containsKey(edgeID) && !edgesIDs.containsKey(edgeID2)) {
					final QuantRatio uniquePepRatio4 = uniquePeptides_U2.getConsensusRatio(cond1, cond2);
					String edgeName4 = getUniqueID(proteinNode2) + "-" + uniquePeptides_U2.getKey();
					Map<String, AttributeValueType> attributes4 = getAttributesForEdge(edgeName4, uniquePepRatio4, null,
							classification1Cases, classification2Cases);
					Edge edge = createEdge(++edgeCounter, null, null, attributes4, getUniqueID(proteinNode2),
							uniquePeptides_U2.getKey(), edgeColorU2);
					graph.getEdge().add(edge);
					edgesIDs.put(edgeID, edge);
					edgesIDs.put(edgeID2, edge);
				}
				updateCasesForEdge(edgesIDs.get(edgeID), edgeColorU2, classification1Cases, classification2Cases);

			}
		}

	}

	private void updateCasesForEdge(Edge edge, Color fillColor, Collection<Classification1Case> classification1Cases,
			Collection<Classification2Case> classification2Cases) {
		if (fillColor != null && fillColor.equals(colorManager.getHighlightColor())) {
			edge.getGraphics().setFill(ColorManager.getHexString(fillColor));
		}
		final List<edu.scripps.yates.pcq.xgmml.jaxb.Graph.Edge.Graphics.Att> atts = edge.getGraphics().getAtt();
		String tooltip = null;
		Set<Classification1Case> cases1 = new HashSet<Classification1Case>();
		Set<Classification2Case> cases2 = new HashSet<Classification2Case>();
		for (edu.scripps.yates.pcq.xgmml.jaxb.Graph.Edge.Graphics.Att att : atts) {
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
				if (classification1Cases != null) {
					cases1.addAll(classification1Cases);
				}
				if (classification2Cases != null) {
					cases2.addAll(classification2Cases);
				}
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
			for (edu.scripps.yates.pcq.xgmml.jaxb.Graph.Edge.Graphics.Att att : atts) {
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
					|| classification1Case == Classification1Case.CASE0) {
				continue;
			}
			if (!list.contains(classification1Case))
				list.add(classification1Case);
		}

		Collections.sort(list, new Comparator<Classification1Case>() {
			@Override
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
			if (classification2Case == Classification2Case.CASE7 || classification2Case == Classification2Case.CASE6) {
				continue;
			}
			if (!list.contains(classification2Case))
				list.add(classification2Case);
		}

		Collections.sort(list, new Comparator<Classification2Case>() {

			@Override
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
	 * and ions counts per each one.<br>
	 * It also includes the consensus ratio value and associated statistics.
	 *
	 * @param prefix
	 * @param peptideNode
	 * @return
	 */
	private String getPeptideNodeTooltip(String prefix, PCQPeptideNode peptideNode) {
		StringBuilder sb = new StringBuilder();
		if (peptideNode.isDiscarded()) {
			sb.append("<b>Peptide node discarded by applied filters</b>\n");
		}
		sb.append(peptideNode.getQuantifiedPeptides().size() + " Peptide sequences\n");

		sb.append(peptideNode.getSequence() + "\n" + peptideNode.getQuantifiedPSMs().size() + " PSMs\n" + "Shared by "
				+ peptideNode.getProteinNodes().size() + " protein Nodes\n" + "Shared by "
				+ PCQUtils.getIndividualProteinsMap(peptideNode).size() + " proteins\n" + "Detected in "
				+ peptideNode.getRawFileNames().size() + " different MS runs\n" + "Detected in "
				+ peptideNode.getFileNames().size() + " different Replicates\n");
		final Double confidenceValue = peptideNode.getConfidenceValue();
		if (confidenceValue != null) {
			sb.append(WEIGHT + ":\t" + formatNumberMoreDecimals(confidenceValue) + "\n");

			sb.append(VARIANCE + ": \t" + formatNumberMoreDecimals(1.0 / confidenceValue) + "\n");
		}
		final QuantRatio consensusRatio = peptideNode.getConsensusRatio(cond1, cond2);
		if (consensusRatio != null) {
			if (consensusRatio instanceof IonCountRatio) {
				final String ionCountRatioTooltip = getIonCountRatioTooltip((IonCountRatio) consensusRatio,
						peptideNode.getQuantifiedPeptides());
				if (ionCountRatioTooltip != null) {
					sb.append(ionCountRatioTooltip + "\n");
				}
			} else {
				sb.append("Ratio (log2):\t" + formatNumber(consensusRatio.getLog2Ratio(cond1, cond2)) + "\n");
				Set<IsobaricQuantifiedPeptide> isobaricQuantifiedPeptides = new HashSet<IsobaricQuantifiedPeptide>();
				for (QuantifiedPeptideInterface peptide : peptideNode.getQuantifiedPeptides()) {
					if (peptide instanceof IsobaricQuantifiedPeptide) {
						isobaricQuantifiedPeptides.add((IsobaricQuantifiedPeptide) peptide);
					}
				}
				if (!isobaricQuantifiedPeptides.isEmpty()) {
					final IonCountRatio consensusIonCountRatio = PCQUtils
							.getConsensusIonCountRatio(isobaricQuantifiedPeptides, cond1, cond2);
					final String ionCountRatioTooltip = getIonCountRatioTooltip(consensusIonCountRatio,
							peptideNode.getQuantifiedPeptides());
					if (ionCountRatioTooltip != null) {
						sb.append(ionCountRatioTooltip + "\n");
					}
				}
			}

			if (consensusRatio.getAssociatedConfidenceScore() != null) {
				sb.append(consensusRatio.getAssociatedConfidenceScore().getScoreName() + ":\t"
						+ consensusRatio.getAssociatedConfidenceScore().getValue() + "\n");
			}

		} else {
			sb.append("No ratio calculated");
		}

		return sb.toString();
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
	private String getIsobaricPeptidesTooltip(String log2ratioValue, Collection<QuantifiedPeptideInterface> peptides) {
		StringBuilder sb = new StringBuilder();
		StringBuilder totalNumerator = new StringBuilder();
		StringBuilder totalDenominator = new StringBuilder();
		for (QuantifiedPeptideInterface peptide : PCQUtils.getSortedPeptidesBySequence(peptides)) {
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
				// final Map<QuantificationLabel, Set<Ion>> ionsByLabel =
				// isoPeptide.getIonsByLabel();
				// StringBuilder ionNumberString = new StringBuilder();
				// if (ionsByLabel != null && !ionsByLabel.isEmpty()) {
				// for (QuantificationLabel label :
				// QuantificationLabel.values()) {
				// if (ionsByLabel.containsKey(label)) {
				// if (!"".equals(ionNumberString.toString())) {
				// ionNumberString.append(", ");
				// }
				// ionNumberString.append(label.name() + ":" +
				// ionsByLabel.get(label).size());
				// }
				// }
				// }

				sb.append(peptide.getSequence() + "\t" + peptide.getQuantifiedPSMs().size() + " PSMs\t"
				// + ionNumberString.toString()
						+ "\t Shared by " + peptide.getQuantifiedProteins().size() + " proteins");
			}
		}
		return "IonCountRatio (Log2):" + log2ratioValue + " = log2( (" + totalNumerator + ") / (" + totalDenominator
				+ ") )\n" + sb.toString();
	}

	private String getIonCountRatioTooltip(IonCountRatio ionCountRatio, Set<QuantifiedPeptideInterface> peptides) {
		final Double log2Ratio = ionCountRatio.getLog2Ratio(cond1, cond2);
		if (log2Ratio != null) {
			return getIsobaricPeptidesTooltip(formatNumber(log2Ratio), peptides);
		}
		return null;
	}

	private String formatNumber(Double number) {
		if (number == null) {
			return null;
		}
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

	private String formatNumberMoreDecimals(Double number) {
		if (number == null) {
			return null;
		}
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
		return formatter3Decimals.format(number);

	}

	private String getUniqueID(Collection<PCQPeptideNode> peptides) {
		return PCQUtils.getPeptideNodesSequenceString(peptides);
	}

	private String getUniqueID(PCQProteinNode protein) {
		if (ProteinClusterQuantParameters.getInstance().getProteinLabel() == ProteinNodeLabel.ACC) {
			return protein.getAccession();
		} else {
			return getProteinNameString(protein);
		}
	}

	/**
	 * Get the label for a protein node, depending on the 'getProteinLabel()'
	 * from the input parameters
	 *
	 * @param protein
	 * @return
	 */
	private String getProteinNodeLabel(PCQProteinNode protein) {
		if (ProteinClusterQuantParameters.getInstance().getProteinLabel() == ProteinNodeLabel.ACC) {
			return protein.getAccession();
		} else if (ProteinClusterQuantParameters.getInstance().getProteinLabel() == ProteinNodeLabel.ID) {
			return getProteinNameString(protein);
		} else {
			return getGeneString(protein);
		}
	}

	private static String controlProteinNodeIDLength(String id) {
		try {
			if (numProteinByProteinNodeID.containsKey(id)) {
				return numProteinByProteinNodeID.get(id);
			}

			if (id.length() > 100) {
				final String string2 = "Prot_" + numProteinNodes;
				numProteinByProteinNodeID.put(id, string2);
				return string2;
			} else {
				numProteinByProteinNodeID.put(id, id);
				return id;
			}
		} finally {
			numProteinNodes++;
		}
	}

	private Node createNodeFromPeptideNode(String nodeID, String label, String tooltip, PCQPeptideNode peptideNode,
			Color outlineColor) {
		Color fillColor = Color.cyan;
		final String sequenceString = peptideNode.getSequence();
		Map<String, AttributeValueType> attributes = new HashMap<String, AttributeValueType>();
		attributes.put("PeptideSequences", new AttributeValueType(sequenceString));
		attributes.put(PCQ_ID, new AttributeValueType(sequenceString));
		final int numIndividualProteins = PCQUtils.getIndividualProteinsMap(peptideNode).size();
		attributes.put("numProteins", new AttributeValueType(numIndividualProteins));
		attributes.put("numPsms", new AttributeValueType(peptideNode.getQuantifiedPSMs().size()));
		attributes.put("numMSRuns", new AttributeValueType(peptideNode.getRawFileNames().size()));
		attributes.put("numReplicates", new AttributeValueType(peptideNode.getFileNames().size()));
		attributes.put("numPeptideSequences", new AttributeValueType(peptideNode.getQuantifiedPeptides().size()));
		attributes.put("numConnectedProteinNodes", new AttributeValueType(peptideNode.getProteinNodes().size()));
		attributes.put("ionCount", new AttributeValueType(PCQUtils.getIonCount(peptideNode)));
		final boolean discarded = peptideNode.isDiscarded();
		attributes.put(IS_FILTERED, new AttributeValueType(getNumFromBoolean(discarded)));
		attributes.put("isProtein", new AttributeValueType(getNumFromBoolean(false)));

		final QuantRatio pepRatio = peptideNode.getConsensusRatio(cond1, cond2);
		final Double finalRatioValue = PCQUtils.getRatioValue(pepRatio, cond1, cond2);
		attributes.put(FINAL_RATIO, new AttributeValueType(finalRatioValue));
		int significant = 0;
		String label_sufix = "";
		final ProteinClusterQuantParameters params = ProteinClusterQuantParameters.getInstance();
		if (pepRatio != null) {
			if (pepRatio instanceof IonCountRatio) {
				attributes.put("normLightIons", new AttributeValueType(((IonCountRatio) pepRatio).getIonCount(cond1)));
				attributes.put("normHeavyIons", new AttributeValueType(((IonCountRatio) pepRatio).getIonCount(cond2)));
			}
			attributes.put("lightIons",
					new AttributeValueType(PCQUtils.getIonCountFromPeptideNode(peptideNode, cond1)));
			attributes.put("heavyIons",
					new AttributeValueType(PCQUtils.getIonCountFromPeptideNode(peptideNode, cond2)));
			final Score associatedConfidenceScore = pepRatio.getAssociatedConfidenceScore();
			if (associatedConfidenceScore != null) {
				attributes.put(associatedConfidenceScore.getScoreName(),
						new AttributeValueType(associatedConfidenceScore.getValue()));

				if (params.getSignificantFDRThreshold() != null && params.getSignificantFDRThreshold() >= Double
						.valueOf(associatedConfidenceScore.getValue())) {

					if (params.getColorNonRegulated() != null) {
						fillColor = params.getColorNonRegulated();
					}
					significant = 1;
				}
			}
			if (finalRatioValue != null && Double.isInfinite(finalRatioValue)) {
				significant = 1;
			}
		}
		if (significant == 1) {
			label_sufix = "*";
		}
		Color labelColor = Color.black;
		if (discarded) {
			fillColor = colorManager.getFillColorForDiscardedNode();
			labelColor = colorManager.getLabelColorForDiscardedNode();
		}
		attributes.put(SIGNIFICANTLY_REGULATED_ATTRIBUTE, new AttributeValueType(significant));
		if (peptideNode.getConfidenceValue() != null) {
			attributes.put(WEIGHT, new AttributeValueType(peptideNode.getConfidenceValue()));
			attributes.put(VARIANCE, new AttributeValueType(1.0 / peptideNode.getConfidenceValue()));
		}

		BORDER_TYPE borderType = BORDER_TYPE.SOLID;

		if (numIndividualProteins > 1) {
			attributes.put("uniquePeptide", new AttributeValueType(0));
		} else {
			attributes.put("uniquePeptide", new AttributeValueType(1));
		}
		if (peptideNode.getProteinNodes().size() > 1) {
			attributes.put("uniquePeptideNode", new AttributeValueType(0));
		} else {
			attributes.put("uniquePeptideNode", new AttributeValueType(1));
		}
		Node node = createNode(nodeID, sequenceString, attributes);

		node.setGraphics(createGraphicsNode(label + label_sufix, getHtml(tooltip), params.getPeptideNodeShape(),
				params.getPeptideNodeHeight(), params.getPeptideNodeWidth(), outlineColor, fillColor, labelColor,
				borderType));
		return node;
	}

	private int getNumFromBoolean(boolean b) {
		if (b) {
			return 1;
		} else {
			return 0;
		}
	}

	private Node createNodeFromProteinNode(PCQProteinNode proteinNode,
			Collection<Classification1Case> classification1Cases, Collection<Classification2Case> classification2Cases,
			Color outlineColor) {
		Map<String, AttributeValueType> attributes = new HashMap<String, AttributeValueType>();
		attributes.put("UniprotKB", new AttributeValueType(proteinNode.getAccession(), AttType.string));
		int numdifferentAccs = 1;
		if (proteinNode.getAccession().contains(PCQUtils.PROTEIN_ACC_SEPARATOR)) {
			numdifferentAccs = proteinNode.getAccession().split(PCQUtils.PROTEIN_ACC_SEPARATOR).length;
		}
		attributes.put(IS_FILTERED, new AttributeValueType(getNumFromBoolean(proteinNode.isDiscarded())));
		attributes.put("numProteins", new AttributeValueType(numdifferentAccs));
		attributes.put("isProtein", new AttributeValueType(getNumFromBoolean(true)));
		attributes.put("numPeptideSequencesInProteins",
				new AttributeValueType(proteinNode.getQuantifiedPeptides().size()));

		attributes.put("numPsmsInProtein", new AttributeValueType(proteinNode.getQuantifiedPSMs().size()));
		attributes.put("numConnectedPeptideNodes", new AttributeValueType(proteinNode.getQuantifiedPeptides().size()));

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

		attributes.put("ProteinDescription",
				new AttributeValueType(proteinNode.getDescription().replace(PCQUtils.PROTEIN_DESCRIPTION_SEPARATOR,
						" " + PCQUtils.PROTEIN_DESCRIPTION_SEPARATOR + " "), AttType.string));
		if (proteinNode.getTaxonomies() != null) {
			attributes.put("Species", new AttributeValueType(proteinNode.getTaxonomies(), AttType.string));
		}
		final String geneString = getGeneString(proteinNode);
		if (geneString != null && !"".equals(geneString)) {
			attributes.put("GeneName", new AttributeValueType(geneString, AttType.string));
		}
		attributes.put("ID", new AttributeValueType(getProteinNameString(proteinNode), AttType.string));
		attributes.put("ACC", new AttributeValueType(proteinNode.getAccession(), AttType.string));
		BORDER_TYPE borderType = BORDER_TYPE.SOLID;

		Color labelColor = Color.black;
		if (PCQUtils.getUniquePeptideNodes(proteinNode).isEmpty()) {
			attributes.put("conclusiveProteinNode", new AttributeValueType("0"));
		} else {
			attributes.put("conclusiveProteinNode", new AttributeValueType("1"));
		}
		if (proteinNode.getQuantifiedPeptides().size() == 1) {
			attributes.put("conclusiveProtein", new AttributeValueType("1"));
		} else {
			attributes.put("conclusiveProtein", new AttributeValueType("0"));
		}

		final String label = controlProteinNodeIDLength(getProteinNodeLabel(proteinNode));
		Node node = createNode(getUniqueID(proteinNode), label, attributes);

		String tooltipText = getProteinNodeTooltip(proteinNode, geneString, classification1Cases, classification1String,
				classification2Cases, classification2String);

		final String tooltip = getHtml(tooltipText);
		node.setGraphics(createGraphicsNode(node.getLabel(), tooltip,
				ProteinClusterQuantParameters.getInstance().getProteinNodeShape(),
				ProteinClusterQuantParameters.getInstance().getProteinNodeHeight(),
				ProteinClusterQuantParameters.getInstance().getProteinNodeWidth(), outlineColor,
				getFillColorByTaxonomy(proteinNode), labelColor, borderType));
		return node;
	}

	private String getProteinNodeTooltip(PCQProteinNode proteinNode, String geneString,
			Collection<Classification1Case> classification1Cases, StringBuilder classification1String,
			Collection<Classification2Case> classification2Cases, StringBuilder classification2String) {
		String tooltipText = "<b>Protein ACC(s):</b>\n" + proteinNode.getAccession() + "\n<b>Protein name(s):</b>\n "
				+ proteinNode.getDescription().replace(PCQUtils.PROTEIN_DESCRIPTION_SEPARATOR, "\n");

		if (classification1Cases != null) {
			tooltipText += "\n<b>Classification1 case:</b> " + classification1String.toString();
		}
		if (classification2Cases != null) {
			tooltipText += "\n<b>Classification2 case:</b> " + classification2String.toString();
		}

		tooltipText += "\n<b>TAX:</b> " + proteinNode.getTaxonomies() + "\n<b>Gene name:</b> " + geneString;
		return tooltipText;
	}

	private String getGeneString(PCQProteinNode proteinNode) {
		final String geneNameString = PCQUtils.getGeneNameString(getAnnotatedProtein(proteinNode.getAccession()),
				proteinNode, null, false);
		return geneNameString;
	}

	private String getProteinNameString(PCQProteinNode proteinNode) {
		String accString = proteinNode.getAccession();
		List<String> list = new ArrayList<String>();
		if (accString.contains(PCQUtils.PROTEIN_ACC_SEPARATOR)) {
			final String[] split = accString.split(PCQUtils.PROTEIN_ACC_SEPARATOR);
			for (String acc : split) {
				final String proteinName = getProteinNameFromUniprot(acc);
				if (proteinName != null) {
					list.add(proteinName);
				} else {
					list.add(acc);
				}
			}
		} else {
			final String proteinName = getProteinNameFromUniprot(accString);
			if (proteinName != null) {
				list.add(proteinName);
			} else {
				list.add(accString);
			}
		}
		StringBuilder sb = new StringBuilder();
		for (String id : list) {
			if (id == null) {
				continue;
			}
			if (!"".equals(sb.toString())) {
				sb.append(PCQUtils.PROTEIN_ACC_SEPARATOR);
			}
			sb.append(id);
		}
		return sb.toString();
	}

	/**
	 * Get protein name from uniprot entry, that is the Uniprot ID, like
	 * ALDOA_HUMAN
	 *
	 * @param acc
	 * @return
	 */
	private String getProteinNameFromUniprot(String acc) {

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

	private Color getFillColorByTaxonomy(PCQProteinNode proteinNode) {
		if (proteinNode.isDiscarded()) {
			return colorManager.getFillColorForDiscardedNode();
		}
		Set<String> taxonomies = proteinNode.getTaxonomies();
		if (taxonomies != null) {
			if (taxonomies.size() > 1) {
				Color ret = colorManager.getMultiTaxonomyColor();
				if (ret != null) {
					return ret;
				}
			} else {
				return colorManager.getColorForProteinTaxonomy(taxonomies.iterator().next());
			}
		}
		return Color.WHITE;

	}

	// private Color getOutlineColorByCase1(Classification1Case
	// classificationCase) {
	// return colorManager.getColorByClassification1Case(classificationCase);
	// }
	//
	// private Color getOutlineColorByCase2(Classification2Case
	// classificationCase) {
	// return colorManager.getColorByClassification2Case(classificationCase);
	// }

	private edu.scripps.yates.pcq.xgmml.jaxb.Graph.Node.Graphics createGraphicsNode(String label, String tooltip,
			Shape shape, int heigth, int width, Color outlineColor, Color fillColor, Color labelColor,
			BORDER_TYPE borderType) {
		edu.scripps.yates.pcq.xgmml.jaxb.Graph.Node.Graphics ret = factory.createGraphNodeGraphics();
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

	private Map<String, AttributeValueType> getAttributesForEdge(String edgeName, QuantRatio sharedPepRatio,
			NWResult alignmentResult, Collection<Classification1Case> classification1Cases,
			Collection<Classification2Case> classification2Cases) {
		Map<String, AttributeValueType> ret = new HashMap<String, AttributeValueType>();
		// ret.put("name", edgeName);
		ret.put(PCQ_ID, new AttributeValueType(edgeName));
		if (sharedPepRatio != null) {
			final Double log2CountRatio = sharedPepRatio.getLog2Ratio(cond1, cond2);
			if (log2CountRatio != null) {
				ret.put(COUNT_RATIO, new AttributeValueType(log2CountRatio));
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

	private edu.scripps.yates.pcq.xgmml.jaxb.Graph.Edge.Graphics createGraphicsEdge(String label, String tooltip,
			Color fillColor) {
		final edu.scripps.yates.pcq.xgmml.jaxb.Graph.Edge.Graphics edge = factory.createGraphEdgeGraphics();
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

	private edu.scripps.yates.pcq.xgmml.jaxb.Graph.Edge.Att createEdgeAttribute(String name, String value,
			String type) {
		final edu.scripps.yates.pcq.xgmml.jaxb.Graph.Edge.Att ret = factory.createGraphEdgeAtt();
		ret.setName(name);
		ret.setValue(value);
		ret.setType(type);
		return ret;
	}

	private edu.scripps.yates.pcq.xgmml.jaxb.Graph.Edge.Graphics.Att createEdgeGraphAttribute(String name, String value,
			String type) {
		edu.scripps.yates.pcq.xgmml.jaxb.Graph.Edge.Graphics.Att ret = factory.createGraphEdgeGraphicsAtt();
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

	private edu.scripps.yates.pcq.xgmml.jaxb.Graph.Node.Att createNodeAtt(String name, String value, String type) {
		final edu.scripps.yates.pcq.xgmml.jaxb.Graph.Node.Att ret = factory.createGraphNodeAtt();
		ret.setName(name);
		ret.setValue(value);
		ret.setType(type);
		return ret;

	}

	private edu.scripps.yates.pcq.xgmml.jaxb.Graph.Node.Graphics.Att createNodeGraphicAtt(String name, String value,
			String type) {
		final edu.scripps.yates.pcq.xgmml.jaxb.Graph.Node.Graphics.Att ret = factory.createGraphNodeGraphicsAtt();
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

	private edu.scripps.yates.pcq.xgmml.jaxb.Graph.Graphics.Att createGraphicAtt(String name, String value,
			String type) {
		final edu.scripps.yates.pcq.xgmml.jaxb.Graph.Graphics.Att ret = factory.createGraphGraphicsAtt();
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
			if (annotatedProteins != null && annotatedProteins.containsKey(string)) {
				ret.put(string, annotatedProteins.get(string));
			}
		}

		return ret;
	}

	/**
	 * @param annotatedProteins
	 *            the annotatedProteins to set
	 */
	public void setAnnotatedProteins(Map<String, Entry> annotatedProteins) {
		this.annotatedProteins = annotatedProteins;
	}

	public void exportToXGMMLUsingNodes(Collection<ProteinCluster> clusterCollection,
			Map<String, Entry> annotatedProteins, QuantCondition condition1, QuantCondition condition2) {
		setAnnotatedProteins(annotatedProteins);

		final ProteinClusterQuantParameters params = ProteinClusterQuantParameters.getInstance();
		ColorManager colorManager = params.getColorManager();
		File outputFileFolder = params.getTemporalOutputFolder();
		String outputPrefix = params.getOutputPrefix();
		String outputSuffix = params.getOutputSuffix();

		log.info("Creating XGMML files for Cytoscape...");
		try {
			log.info("Creating XGMML for the entire network...");
			// export the total network
			File xgmmlOutPutFile = new File(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix
					+ "_cytoscape_ALL_" + outputSuffix + ".xgmml");
			exportToGmmlFromProteinClustersUsingNodes(xgmmlOutPutFile, outputPrefix + "_" + outputSuffix,
					clusterCollection, condition1, condition2, colorManager);

			// fdr
			Double fdrThreshold = params.getSignificantFDRThreshold();
			log.info("Creating XGMML for the cluster containing a peptide node that is significantly changing...");
			Set<ProteinCluster> significantlyRegulatedProteinClusters = getSignificantlyRegulatedProteinClusters(
					clusterCollection, fdrThreshold);
			String fdrText = "";
			if (params.isPerformRatioIntegration() && params.getSignificantFDRThreshold() != null) {
				fdrText = params.getSignificantFDRThreshold() + "_";
			}
			final String fileName = outputFileFolder.getAbsolutePath() + File.separator + outputPrefix
					+ "_cytoscape_Significants_" + fdrText + outputSuffix + ".xgmml";

			File xgmmlOutPutFileFDR = new File(fileName);
			exportToGmmlFromProteinClustersUsingNodes(xgmmlOutPutFileFDR,
					outputPrefix + "_FDR" + fdrThreshold + "_" + outputSuffix, significantlyRegulatedProteinClusters,
					condition1, condition2, colorManager);

			// classification 1
			final Map<Classification1Case, Set<ProteinCluster>> proteinPairsByClassification1 = getProteinClustersByClassification1(
					clusterCollection);
			final Classification1Case[] cases1 = Classification1Case.values();
			for (Classification1Case case1 : cases1) {
				if (proteinPairsByClassification1.containsKey(case1)) {
					log.info("Creating XGMML for case " + case1.getCaseID() + "(" + case1.getExplanation() + ")...");

					File xgmmlOutPutFile2 = new File(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix
							+ "_cytoscape_" + case1.getCaseID() + "-" + case1.name() + "_" + outputSuffix + ".xgmml");
					exportToGmmlFromProteinClustersUsingNodes(xgmmlOutPutFile2,
							outputPrefix + "_" + case1.getCaseID() + "-" + case1.name() + "_" + outputSuffix,
							proteinPairsByClassification1.get(case1), condition1, condition2, colorManager);
				}
			}
			// classification 2
			final Map<Classification2Case, Set<ProteinCluster>> proteinPairsByClassification2 = getProteinClustersByClassification2(
					clusterCollection);
			final Classification2Case[] cases2 = Classification2Case.values();
			for (Classification2Case case2 : cases2) {
				if (proteinPairsByClassification2.containsKey(case2)) {
					log.info("Creating XGMML for case " + case2.getCaseID() + "(" + case2.getExplanation() + ")...");
					File xgmmlOutPutFile2 = new File(outputFileFolder.getAbsolutePath() + File.separator + outputPrefix
							+ "_cytoscape_" + case2.getCaseID() + "-" + case2.name() + "_" + outputSuffix + ".xgmml");
					exportToGmmlFromProteinClustersUsingNodes(xgmmlOutPutFile2,
							outputPrefix + "_" + case2.getCaseID() + "-" + case2.name() + "_" + outputSuffix,
							proteinPairsByClassification2.get(case2), condition1, condition2, colorManager);
				}
			}

		} catch (JAXBException e) {
			e.printStackTrace();
			log.error(e.getMessage());
		}

	}

	/**
	 * Get significantly regulated protein clusters, that is, the ones having at
	 * least one peptide node with FDR less or equals to the input parameter
	 * FDRThredhold, or having an INFINITY value
	 *
	 * @param clusters
	 * @param fdrThreshold
	 * @return
	 */
	private Set<ProteinCluster> getSignificantlyRegulatedProteinClusters(Collection<ProteinCluster> clusters,
			Double fdrThreshold) {
		Set<ProteinCluster> ret = new HashSet<ProteinCluster>();
		for (ProteinCluster proteinCluster : clusters) {
			final Set<PCQPeptideNode> peptideNodes = proteinCluster.getPeptideNodes();
			for (PCQPeptideNode peptideNode : peptideNodes) {
				final QuantRatio consensusRatio = peptideNode.getConsensusRatio(cond1, cond2);
				if (consensusRatio != null) {
					final Score score = consensusRatio.getAssociatedConfidenceScore();
					if (score != null) {
						try {
							if (score.getValue() != null) {
								double fdr = Double.valueOf(score.getValue());
								if (fdrThreshold != null && fdr <= fdrThreshold) {
									ret.add(proteinCluster);
									break;
								}
							}

						} catch (NumberFormatException e) {

						}
					}
					if (Double.isInfinite(consensusRatio.getLog2Ratio(cond1, cond2))) {
						ret.add(proteinCluster);
						break;
					}
				}
			}
		}
		return ret;
	}

	private Map<Classification1Case, Set<ProteinCluster>> getProteinClustersByClassification1(
			Collection<ProteinCluster> clusters) {
		Map<Classification1Case, Set<ProteinCluster>> ret = new HashMap<Classification1Case, Set<ProteinCluster>>();
		for (ProteinCluster proteinCluster : clusters) {
			final Set<ProteinPair> proteinPairs = proteinCluster.getProteinPairs();
			for (ProteinPair proteinPair : proteinPairs) {
				final Collection<Classification1Case> classification1PairCases = proteinPair.getClassification1Case()
						.values();
				for (Classification1Case classification1PairCase : classification1PairCases) {
					if (ret.containsKey(classification1PairCase)) {
						ret.get(classification1PairCase).add(proteinCluster);
					} else {
						Set<ProteinCluster> set = new HashSet<ProteinCluster>();
						set.add(proteinCluster);
						ret.put(classification1PairCase, set);
					}
				}
			}
		}
		return ret;
	}

	private Map<Classification2Case, Set<ProteinCluster>> getProteinClustersByClassification2(
			Collection<ProteinCluster> clusters) {
		Map<Classification2Case, Set<ProteinCluster>> ret = new HashMap<Classification2Case, Set<ProteinCluster>>();
		for (ProteinCluster proteinCluster : clusters) {
			final Set<ProteinPair> proteinPairs = proteinCluster.getProteinPairs();
			for (ProteinPair proteinPair : proteinPairs) {
				final Collection<Classification2Case> classification2PairCases = proteinPair.getClassification2Cases()
						.values();
				for (Classification2Case classification2PairCase : classification2PairCases) {
					if (ret.containsKey(classification2PairCase)) {
						ret.get(classification2PairCase).add(proteinCluster);
					} else {
						Set<ProteinCluster> set = new HashSet<ProteinCluster>();
						set.add(proteinCluster);
						ret.put(classification2PairCase, set);
					}
				}
			}
		}
		return ret;
	}

}
