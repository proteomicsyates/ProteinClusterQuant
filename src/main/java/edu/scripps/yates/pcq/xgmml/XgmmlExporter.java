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
import edu.scripps.yates.census.read.model.IonCountRatio;
import edu.scripps.yates.census.read.model.IsobaricQuantifiedPeptide;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.pcq.cases.Classification1Case;
import edu.scripps.yates.pcq.cases.Classification2Case;
import edu.scripps.yates.pcq.model.PCQPeptideNode;
import edu.scripps.yates.pcq.model.PCQProteinNode;
import edu.scripps.yates.pcq.model.ProteinCluster;
import edu.scripps.yates.pcq.model.ProteinPair;
import edu.scripps.yates.pcq.params.ProteinClusterQuantParameters;
import edu.scripps.yates.pcq.util.AnalysisInputType;
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
import edu.scripps.yates.pcq.xgmml.util.UniprotAnnotationColumn;
import edu.scripps.yates.utilities.alignment.nwalign.NWResult;
import edu.scripps.yates.utilities.colors.ColorGenerator;
import edu.scripps.yates.utilities.proteomicsmodel.Score;
import edu.scripps.yates.utilities.sequence.PositionInPeptide;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

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
	private final Map<String, Node> visitedPeptideKeys = new THashMap<String, Node>();
	private final Map<PCQProteinNode, Node> visitedProteinNodes = new THashMap<PCQProteinNode, Node>();
	private final Map<String, Edge> edgesIDs = new THashMap<String, Edge>();
	private Map<String, Entry> annotatedProteins;
	private static final String PCQ_ID = "PCQ_ID";
	private static final String COUNT_RATIO = "countRatio";
	private static final String FINAL_RATIO = "finalRatio";
	private static int numProteinNodes = 1;
	private static final Map<String, String> numProteinByProteinNodeID = new THashMap<String, String>();
	private static final Set<String> ratioAttributes = new THashSet<String>();
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
						node2.getGraphics().setFill(ColorGenerator.getHexString(color));
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

						node2.getGraphics().setFill(ColorGenerator.getHexString(color));
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
		if (!proteinCluster.getAllProteinPairs().isEmpty()) {
			createNodesAndEdgesFromProteinPairsUsingNodes(proteinCluster.getProteinPairs(), graph);
			// now the proteins and peptides without protein pair because they
			// were discarded probably
			for (PCQProteinNode pcqProteinNode : proteinNodes) {
				if (pcqProteinNode.getProteinPair() == null) {
					createNodesAndEdgesFromProteinNodes(pcqProteinNode, null, false, false, false, null, graph);
				}
			}

		} else {
			createNodesAndEdgesFromProteinNodes(proteinNodes.iterator().next(), null, false, false, false, null, graph);
		}
		// draw an edge between the aligned peptides
		final Set<PCQPeptideNode> peptideNodes = proteinCluster.getPeptideNodes();
		for (PCQPeptideNode peptideNode : peptideNodes) {
			final Set<PCQPeptideNode> alignedPeptides = proteinCluster.getAlignedPeptideNodes(peptideNode);

			for (PCQPeptideNode peptideNode2 : alignedPeptides) {

				Set<PCQPeptideNode> peptidesToLink = new THashSet<PCQPeptideNode>();
				peptidesToLink.add(peptideNode);
				peptidesToLink.add(peptideNode2);
				String edgeID = getUniqueID(peptidesToLink);
				if (!edgesIDs.containsKey(edgeID)) {

					String edgeName3 = PCQUtils.getPeptideNodesSequenceString(peptidesToLink);
					final NWResult alignmentResult = proteinCluster.getAlignmentResult(peptideNode, peptideNode2);
					Map<String, AttributeValueType> attributes3 = getAttributesForEdge(edgeName3, null, alignmentResult,
							null);
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
					proteinPair.isUniquePeptidesProt2Inconsistent(), proteinPair.getClassification2Cases().values(),
					graph);
		}
	}

	private void createNodesAndEdgesFromProteinNodes(PCQProteinNode proteinNode1, PCQProteinNode proteinNode2,
			boolean isSharedPeptidesInconsistent, boolean isUniquePeptidesProt1Inconsistent,
			boolean isUniquePeptidesProt2Inconsistent, Collection<Classification2Case> classification2Cases,
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
		if (proteinNode1 != null) {
			final Set<PCQPeptideNode> peptidesNodesU1 = PCQUtils.getUniquePeptideNodes(proteinNode1, proteinNode2, true,
					ProteinClusterQuantParameters.getInstance().isRemoveFilteredNodes());

			for (PCQPeptideNode uniquePeptideNode_U1 : peptidesNodesU1) {
				final String peptidesSequenceString_U1 = uniquePeptideNode_U1.getKey();
				// System.out.println(uniquePeptideNode_U1);
				if (!visitedPeptideKeys.containsKey(peptidesSequenceString_U1)) {
					final String nodeID = uniquePeptideNode_U1.getKey();
					final String label = formatNumber(
							PCQUtils.getLog2RatioValue(PCQUtils.getRepresentativeRatioForPeptideNode(
									uniquePeptideNode_U1, cond1, cond2, true, null), cond1, cond2));
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
				Node node = createNodeFromProteinNode(proteinNode1, classification2Cases, outlineColorP1);
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
							classification2Cases);
					Edge edge = createEdge(++edgeCounter, null, null, attributes3, uniquePeptideNode_U1.getKey(),
							getUniqueID(proteinNode1), edgeColorU1);
					graph.getEdge().add(edge);
					edgesIDs.put(edgeID, edge);
					edgesIDs.put(edgeID2, edge);
				}
				updateCasesForEdge(edgesIDs.get(edgeID), edgeColorU1, classification2Cases);
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
					final String label = formatNumber(
							PCQUtils.getLog2RatioValue(PCQUtils.getRepresentativeRatioForPeptideNode(
									sharedPeptideNode_S12, cond1, cond2, true, null), cond1, cond2));
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
							classification2Cases);

					Edge edge = createEdge(++edgeCounter, null, null, attributes, getUniqueID(proteinNode1),
							sharedPeptideNode_S12.getKey(), edgeColorS1);
					graph.getEdge().add(edge);
					edgesIDs.put(edgeID, edge);
					edgesIDs.put(edgeID2, edge);
				}
				updateCasesForEdge(edgesIDs.get(edgeID), edgeColorS1, classification2Cases);

				if (proteinNode2 != null) {
					// S12P2
					edgeID = getUniqueID(proteinNode2) + sharedPeptideNode_S12.getKey();
					edgeID2 = sharedPeptideNode_S12.getKey() + getUniqueID(proteinNode2);
					if (!edgesIDs.containsKey(edgeID) && !edgesIDs.containsKey(edgeID2)) {
						String edgeName2 = getUniqueID(proteinNode2) + "-" + sharedPeptideNode_S12.getKey();
						Map<String, AttributeValueType> attributes2 = getAttributesForEdge(edgeName2, sharedPepRatio,
								null, classification2Cases);
						Edge edge = createEdge(++edgeCounter, null, null, attributes2, getUniqueID(proteinNode2),
								sharedPeptideNode_S12.getKey(), edgeColorS2);
						graph.getEdge().add(edge);
						edgesIDs.put(edgeID, edge);
						edgesIDs.put(edgeID2, edge);
					}
					updateCasesForEdge(edgesIDs.get(edgeID), edgeColorS2, classification2Cases);

				}

			}
		}
		// P2 protein2
		if (proteinNode2 != null) {
			if (!visitedProteinNodes.containsKey(proteinNode2)) {
				Node node = createNodeFromProteinNode(proteinNode2, classification2Cases, outlineColorP2);
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
					final String label = formatNumber(PCQUtils.getLog2RatioValue(
							PCQUtils.getRepresentativeRatioForPeptideNode(uniquePeptides_U2, cond1, cond2, true, null),
							cond1, cond2));
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
							classification2Cases);
					Edge edge = createEdge(++edgeCounter, null, null, attributes4, getUniqueID(proteinNode2),
							uniquePeptides_U2.getKey(), edgeColorU2);
					graph.getEdge().add(edge);
					edgesIDs.put(edgeID, edge);
					edgesIDs.put(edgeID2, edge);
				}
				updateCasesForEdge(edgesIDs.get(edgeID), edgeColorU2, classification2Cases);

			}
		}

	}

	private void updateCasesForEdge(Edge edge, Color fillColor, Collection<Classification2Case> classification2Cases) {
		if (fillColor != null && fillColor.equals(colorManager.getHighlightColor())) {
			edge.getGraphics().setFill(ColorGenerator.getHexString(fillColor));
		}
		final List<edu.scripps.yates.pcq.xgmml.jaxb.Graph.Edge.Graphics.Att> atts = edge.getGraphics().getAtt();
		String tooltip = "";
		// Set<Classification1Case> cases1 = new
		// THashSet<Classification1Case>();
		Set<Classification2Case> cases2 = new THashSet<Classification2Case>();
		for (edu.scripps.yates.pcq.xgmml.jaxb.Graph.Edge.Graphics.Att att : atts) {
			if (att.getName().equals("EDGE_LABEL")) {
				if (att.getValue() != null) {
					// final String oldValueCases1 =
					// att.getValue().split("\n")[0];
					// if (oldValueCases1 != null &&
					// !"null".equals(oldValueCases1)) {
					// if (oldValueCases1.contains(",")) {
					// final String[] split = oldValueCases1.split(",");
					// for (String string : split) {
					// cases1.add(Classification1Case.getByCaseID(Integer.valueOf(string)));
					// }
					// } else {
					// if (!oldValueCases1.equals(" ")) {
					// cases1.add(Classification1Case.getByCaseID(Integer.valueOf(oldValueCases1)));
					// }
					// }
					// }
					final String oldValueCases2 = att.getValue();
					if (oldValueCases2 != null && !"null".equals(oldValueCases2)) {
						if (oldValueCases2.contains("\n")) {
							final String[] split = oldValueCases2.split("\n");
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
				// if (classification1Cases != null) {
				// cases1.addAll(classification1Cases);
				// }
				if (classification2Cases != null) {
					cases2.addAll(classification2Cases);
				}
				// final String classificationCase1WithExplanationString =
				// getClassificationCase1String(cases1, true);
				// if (classificationCase1WithExplanationString != null) {
				// tooltip = "Classification 1: " +
				// classificationCase1WithExplanationString;
				// }

				final String classificationCase2WithExplanationString = getClassificationCase2String(cases2, true,
						false);
				if (classificationCase2WithExplanationString != null) {
					if (!"".equals(tooltip)) {
						tooltip += "\n";
					}
					tooltip += classificationCase2WithExplanationString;
				}
				// String case1String = null;
				// String case2String = null;
				// final String classificationCase1String =
				// getClassificationCase1String(cases1, false);
				// if (classificationCase1String != null) {
				// case1String = classificationCase1String;
				// }

				if (ProteinClusterQuantParameters.getInstance().isShowCasesInEdges()) {
					// only if significant
					final String classificationCase2String = getClassificationCase2String(cases2, false, true);
					att.setValue(classificationCase2String);
				}
				break;
			}
		}
		if (tooltip != null) {
			for (edu.scripps.yates.pcq.xgmml.jaxb.Graph.Edge.Graphics.Att att : atts) {
				if (att.getName().equals("EDGE_TOOLTIP")) {
					final String html = getHtml("This protein pair has been classified as:\n" + tooltip);
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
			if (!classification1Case.isInconsistence()) {
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
			boolean addExplanationOfCases, boolean onlySignificantOnes) {
		StringBuilder sb = new StringBuilder();
		int index = 0;
		List<Classification2Case> list = new ArrayList<Classification2Case>();
		for (Classification2Case classification2Case : classification2Cases) {
			if (onlySignificantOnes && !classification2Case.isInconsistence()) {
				continue;
			}
			if (!list.contains(classification2Case)) {
				list.add(classification2Case);
			}
		}

		Collections.sort(list, new Comparator<Classification2Case>() {

			@Override
			public int compare(Classification2Case o1, Classification2Case o2) {
				return Integer.compare(o1.getCaseID(), o2.getCaseID());
			}
		});
		for (Classification2Case case2 : list) {
			if (index > 0)
				sb.append("\n");
			if (addExplanationOfCases && case2.isInconsistence()) {
				sb.append("<b>");
			}

			sb.append(case2.getCaseID());
			if (addExplanationOfCases) {
				sb.append(" (" + case2.getExplanation() + ") ");
			}
			if (addExplanationOfCases && case2.isInconsistence()) {
				sb.append("</b>");
			}
			index++;
		}
		if ("".equals(sb.toString()))
			return null;
		return sb.toString();
	}

	private void setNodeOutlineColor(Node node, Color outlineColor) {
		node.getGraphics().setOutline(ColorGenerator.getHexString(outlineColor));
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
		final Set<QuantifiedPeptideInterface> quantifiedPeptides = peptideNode.getQuantifiedPeptides();
		if (!peptideNode.getKey().contentEquals(peptideNode.getFullSequence())) {
			sb.append(peptideNode.getKey() + "\n");
		}
		sb.append(getSequenceAnnotated(peptideNode) + "\n");

		sb.append(quantifiedPeptides.size() + " Peptide sequences\n" + peptideNode.getQuantifiedPSMs().size()
				+ " PSMs\n" + "Shared by " + peptideNode.getProteinNodes().size() + " protein Nodes\n" + "Shared by "
				+ PCQUtils.getIndividualProteinsMap(peptideNode).size() + " proteins\n" + "Detected in "
				+ peptideNode.getRawFileNames().size() + " MS runs\n" + "Detected in "
				+ peptideNode.getFileNames().size() + " Replicates\n");
		final Double confidenceValue = peptideNode.getConfidenceValue();
		if (confidenceValue != null) {
			sb.append(WEIGHT + " = " + formatNumberMoreDecimals(confidenceValue) + "\n");

			sb.append(VARIANCE + " = " + formatNumberMoreDecimals(1.0 / confidenceValue) + "\n");
		}
		QuantRatio finalRatio = PCQUtils.getRepresentativeRatioForPeptideNode(peptideNode, cond1, cond2, true, null);

		if (finalRatio != null) {
			if (finalRatio instanceof IonCountRatio) {
				final String ionCountRatioTooltip = getIonCountRatioTooltip((IonCountRatio) finalRatio,
						quantifiedPeptides);
				if (ionCountRatioTooltip != null) {
					sb.append(ionCountRatioTooltip + "\n");
				}
			} else {
				sb.append(finalRatio.getDescription() + " = " + formatNumber(finalRatio.getLog2Ratio(cond1, cond2))
						+ "\n");

				final IonCountRatio ionCountRatio = PCQUtils.getNormalizedIonCountRatioForPeptideNode(peptideNode,
						cond1, cond2, null);
				if (ionCountRatio != null) {
					final String ionCountRatioTooltip = getIonCountRatioTooltip(ionCountRatio, quantifiedPeptides);
					if (ionCountRatioTooltip != null) {
						sb.append(ionCountRatioTooltip + "\n");
					}
				}
			}

			if (finalRatio.getAssociatedConfidenceScore() != null) {
				sb.append(finalRatio.getAssociatedConfidenceScore().getScoreName() + " = "
						+ finalRatio.getAssociatedConfidenceScore().getValue() + "\n");
			}

		} else {
			sb.append("No ratio calculated\n");
		}
		sb.append("Individual peptides in the node:\n");
		for (QuantifiedPeptideInterface peptide : PCQUtils.getSortedPeptidesBySequence(quantifiedPeptides)) {
			final QuantRatio individualPeptideRatio = peptide.getConsensusRatio(cond1, cond2);
			sb.append(peptide.getFullSequence() + ", " + peptide.getQuantifiedPSMs().size() + " PSMs, " + "Shared by "
					+ peptide.getQuantifiedProteins().size() + " proteins, " + "Detected in "
					+ peptide.getRawFileNames().size() + " MS Runs, " + "Detected in " + peptide.getFileNames().size()
					+ " replicates, " + individualPeptideRatio.getDescription() + " = "
					+ individualPeptideRatio.getLog2Ratio(cond1, cond2) + "\n");
		}
		if (peptideNode.getTaxonomies() != null && !peptideNode.getTaxonomies().isEmpty()) {
			sb.append("\n<b>TAX:</b> " + PCQUtils.getSpeciesString(peptideNode.getTaxonomies()));
		}
		return sb.toString();
	}

	private String getSequenceAnnotated(PCQPeptideNode peptideNode) {
		if (ProteinClusterQuantParameters.getInstance().isCollapseBySites()) {
			StringBuilder sb = new StringBuilder();
			final List<Pair<IsobaricQuantifiedPeptide, PositionInPeptide>> peptidesWithPositionsInPeptide = peptideNode
					.getPeptidesWithPositionsInPeptide();
			for (Pair<IsobaricQuantifiedPeptide, PositionInPeptide> pair : peptidesWithPositionsInPeptide) {
				if (!"".equals(sb.toString())) {
					sb.append("-");
				}
				final String fullSequence = pair.getFirstelement().getFullSequence();
				int position = pair.getSecondElement().getPosition();
				int currentposition = 0;
				boolean isPTM = false;
				for (int i = 0; i < fullSequence.length(); i++) {

					final char charAt = fullSequence.charAt(i);
					if (charAt == '(' || charAt == '[') {
						isPTM = true;
						continue;
					}
					if (charAt == ')' || charAt == ']') {
						isPTM = false;
						continue;
					}
					if (!isPTM) {
						currentposition++;
						if (currentposition == position) {
							sb.append("<b>" + charAt + "</b>");
						} else {
							sb.append(charAt);
						}
					}
				}
			}
			return sb.toString();
		} else {
			return peptideNode.getFullSequence();
		}
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
		// StringBuilder sb = new StringBuilder();
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

			}
		}
		return "Ion count ratio Rc = " + log2ratioValue + " = log2( (" + totalNumerator + ") / (" + totalDenominator
				+ ") )";
		// + sb.toString();
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
				return "Infinity";
			}
			if (Double.compare(Double.NEGATIVE_INFINITY, number) == 0) {
				return "-Infinity";
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
			return protein.getKey();
		} else {
			return getProteinNameString(protein);
		}
	}

	/**
	 * Get the label for a protein node, depending on the 'getProteinLabel()'
	 * from the input parameters
	 *
	 * @param proteinNode
	 * @return
	 */
	private String getProteinNodeLabel(PCQProteinNode proteinNode) {
		if (ProteinClusterQuantParameters.getInstance().getProteinLabel() == ProteinNodeLabel.ACC) {
			return proteinNode.getKey();
		} else if (ProteinClusterQuantParameters.getInstance().getProteinLabel() == ProteinNodeLabel.ID) {
			return getProteinNameString(proteinNode);
		} else {
			return getGeneString(proteinNode);
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
		final String sequenceString = peptideNode.getFullSequence();
		Map<String, AttributeValueType> attributes = new THashMap<String, AttributeValueType>();
		attributes.put("PeptideSequences", new AttributeValueType(sequenceString));
		attributes.put(PCQ_ID, new AttributeValueType(nodeID));
		final int numIndividualProteins = PCQUtils.getIndividualProteinsMap(peptideNode).size();
		attributes.put("numProteins", new AttributeValueType(numIndividualProteins));
		attributes.put("numPsms", new AttributeValueType(peptideNode.getQuantifiedPSMs().size()));
		attributes.put("numMSRuns", new AttributeValueType(peptideNode.getRawFileNames().size()));
		attributes.put("numReplicates", new AttributeValueType(peptideNode.getFileNames().size()));
		attributes.put("numPeptideSequences", new AttributeValueType(peptideNode.getQuantifiedPeptides().size()));
		attributes.put("numConnectedProteinNodes", new AttributeValueType(peptideNode.getProteinNodes().size()));
		if (peptideNode.getTaxonomies() != null) {
			attributes.put("Species",
					new AttributeValueType(PCQUtils.getSpeciesString(peptideNode.getTaxonomies()), AttType.string));
		}
		attributes.put("ionCount", new AttributeValueType(PCQUtils.getIonCount(peptideNode)));
		final boolean discarded = peptideNode.isDiscarded();
		attributes.put(IS_FILTERED, new AttributeValueType(getNumFromBoolean(discarded)));
		attributes.put("isProtein", new AttributeValueType(getNumFromBoolean(false)));
		attributes.put("containsPTMs", new AttributeValueType(getNumFromBoolean(peptideNode.containsPTMs())));
		final QuantRatio pepRatio = PCQUtils.getRepresentativeRatioForPeptideNode(peptideNode, cond1, cond2, true,
				null);
		final Double finalRatioValue = PCQUtils.getLog2RatioValue(pepRatio, cond1, cond2);
		attributes.put(FINAL_RATIO, new AttributeValueType(finalRatioValue));

		// in case of isotopologues
		if (ProteinClusterQuantParameters.getInstance().getInputType() == AnalysisInputType.CENSUS_CHRO) {
			// distinguish between Ri and Rc
			final Double log2RiRatio = peptideNode.getConsensusRatio(cond1, cond2).getLog2Ratio(cond1, cond2);
			attributes.put("Ri", new AttributeValueType(log2RiRatio));
			attributes.put("Rc",
					new AttributeValueType(
							PCQUtils.getNormalizedIonCountRatioForPeptideNode(peptideNode, cond1, cond2, null)
									.getLog2Ratio(cond1, cond2)));
		}

		int significant = 0;
		String label_sufix = "";
		final ProteinClusterQuantParameters params = ProteinClusterQuantParameters.getInstance();
		if (pepRatio != null) {
			if (pepRatio instanceof IonCountRatio && !pepRatio.getLog2Ratio(cond1, cond2).isNaN()) {
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

				if (params.getSignificantFDRThreshold() != null
						&& associatedConfidenceScore.getScoreName().equals(PCQUtils.FDR_CONFIDENCE_SCORE_NAME) && params
								.getSignificantFDRThreshold() >= Double.valueOf(associatedConfidenceScore.getValue())) {

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
		Node node = createNode(nodeID, nodeID, attributes);

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
			Collection<Classification2Case> classification2Cases, Color outlineColor) {
		Map<String, AttributeValueType> attributes = new THashMap<String, AttributeValueType>();
		attributes.put("UniprotKB", new AttributeValueType(proteinNode.getKey(), AttType.string));
		int numdifferentAccs = 1;
		if (proteinNode.getKey().contains(PCQUtils.PROTEIN_ACC_SEPARATOR)) {
			numdifferentAccs = proteinNode.getKey().split(PCQUtils.PROTEIN_ACC_SEPARATOR).length;
		}
		attributes.put(IS_FILTERED, new AttributeValueType(getNumFromBoolean(proteinNode.isDiscarded())));
		attributes.put("numProteins", new AttributeValueType(numdifferentAccs));
		attributes.put("isProtein", new AttributeValueType(getNumFromBoolean(true)));
		attributes.put("containsPTMs", new AttributeValueType(getNumFromBoolean(proteinNode.containsPTMs())));
		attributes.put("numPeptideSequencesInProteins",
				new AttributeValueType(proteinNode.getQuantifiedPeptides().size()));

		attributes.put("numPsmsInProtein", new AttributeValueType(proteinNode.getQuantifiedPSMs().size()));
		attributes.put("numConnectedPeptideNodes", new AttributeValueType(proteinNode.getQuantifiedPeptides().size()));

		// chack uniprot columsn
		List<UniprotAnnotationColumn> uniprotAnnotationColumns = ProteinClusterQuantParameters.getInstance()
				.getUniprotAnnotationColumns();
		if (!uniprotAnnotationColumns.isEmpty()) {
			for (UniprotAnnotationColumn uniprotAnnotationColumn : uniprotAnnotationColumns) {
				StringBuilder sb = new StringBuilder();
				for (QuantifiedProteinInterface protein : proteinNode.getQuantifiedProteins()) {

					List<String> valuesForProtein = uniprotAnnotationColumn.getValuesForProtein(protein.getAccession(),
							annotatedProteins);
					for (String value : valuesForProtein) {
						if (!"".equals(sb.toString())) {
							sb.append(" ");
						}
						sb.append(value);
					}
				}
				if (!"".equals(sb.toString())) {
					attributes.put(uniprotAnnotationColumn.getColumnName(),
							new AttributeValueType(sb.toString(), AttType.string));
				}
			}
		}

		// StringBuilder classification1String = new StringBuilder();
		// StringBuilder classification1StringNOHTML = new StringBuilder();
		// if (classification1Cases != null) {
		// for (Classification1Case classification1Case : classification1Cases)
		// {
		// final String value = "<b>" +
		// String.valueOf(classification1Case.getCaseID()) + "</b>: "
		// + classification1Case.getExplanation();
		// if (!"".equals(classification1String.toString())) {
		// classification1String.append(", ");
		// }
		// classification1String.append(value);
		// classification1StringNOHTML.append(classification1Case.getCaseID());
		//
		// }
		// attributes.put("Classification1Case",
		// new AttributeValueType(classification1StringNOHTML.toString(),
		// AttType.string));
		// }

		StringBuilder classification2StringTooltip = new StringBuilder();
		StringBuilder classification2StringNOHTML = new StringBuilder();
		Set<Classification2Case> cases2 = new THashSet<Classification2Case>();
		if (classification2Cases != null) {
			List<Classification2Case> sortedList = new ArrayList<Classification2Case>();
			sortedList.addAll(classification2Cases);
			Collections.sort(sortedList, new Comparator<Classification2Case>() {

				@Override
				public int compare(Classification2Case o1, Classification2Case o2) {
					return Integer.compare(o1.getCaseID(), o2.getCaseID());
				}
			});
			for (Classification2Case classification2Case : sortedList) {
				if (cases2.contains(classification2Case)) {
					continue;
				}
				cases2.add(classification2Case);
				final String value = "<b>" + String.valueOf(classification2Case.getCaseID()) + "</b>: "
						+ classification2Case.getExplanation();
				if (!"".equals(classification2StringTooltip.toString())) {
					classification2StringTooltip.append("\n");
					classification2StringNOHTML.append(",");
				}
				classification2StringTooltip.append(value);
				classification2StringNOHTML.append(classification2Case.getCaseID());

			}
			attributes.put("Classification2Case",
					new AttributeValueType(classification2StringNOHTML.toString(), AttType.string));
		}

		attributes.put("ProteinDescription",
				new AttributeValueType(proteinNode.getDescription().replace(PCQUtils.PROTEIN_DESCRIPTION_SEPARATOR,
						" " + PCQUtils.PROTEIN_DESCRIPTION_SEPARATOR + " "), AttType.string));
		if (proteinNode.getTaxonomies() != null) {
			attributes.put("Species",
					new AttributeValueType(PCQUtils.getSpeciesString(proteinNode.getTaxonomies()), AttType.string));
		}
		final String geneString = getGeneString(proteinNode);
		if (geneString != null && !"".equals(geneString)) {
			attributes.put("GeneName", new AttributeValueType(geneString, AttType.string));
		}
		attributes.put("ID", new AttributeValueType(getProteinNameString(proteinNode), AttType.string));
		attributes.put("ACC", new AttributeValueType(proteinNode.getKey(), AttType.string));
		BORDER_TYPE borderType = BORDER_TYPE.SOLID;

		Color labelColor = Color.black;
		if (PCQUtils.getUniquePeptideNodes(proteinNode, true).isEmpty()) {
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

		String tooltipText = getProteinNodeTooltip(proteinNode, geneString, classification2Cases,
				classification2StringTooltip);

		final String tooltip = getHtml(tooltipText);
		node.setGraphics(createGraphicsNode(node.getLabel(), tooltip,
				ProteinClusterQuantParameters.getInstance().getProteinNodeShape(),
				ProteinClusterQuantParameters.getInstance().getProteinNodeHeight(),
				ProteinClusterQuantParameters.getInstance().getProteinNodeWidth(), outlineColor,
				getFillColorByTaxonomy(proteinNode), labelColor, borderType));
		return node;
	}

	private String getProteinNodeTooltip(PCQProteinNode proteinNode, String geneString,
			Collection<Classification2Case> classification2Cases, StringBuilder classification2String) {
		StringBuilder tooltipText = new StringBuilder("<b>Protein ACC(s):</b>\n").append(proteinNode.getKey())
				.append("\n<b>Protein name(s):</b>\n ")
				.append(proteinNode.getDescription().replace(PCQUtils.PROTEIN_DESCRIPTION_SEPARATOR, "\n"));

		if (classification2Cases != null) {
			tooltipText.append("\n<b>Classification case(s):</b>\n").append(classification2String.toString());
		}
		if (proteinNode.getTaxonomies() != null && !proteinNode.getTaxonomies().isEmpty()) {
			tooltipText.append("\n<b>TAX:</b> ").append(PCQUtils.getSpeciesString(proteinNode.getTaxonomies()))
					.append("\n<b>Gene name:</b> ").append(geneString);
		}
		return tooltipText.toString();
	}

	private String getGeneString(PCQProteinNode proteinNode) {
		final String geneNameString = PCQUtils.getGeneNameString(getAnnotatedProtein(proteinNode.getKey()), proteinNode,
				null, false, true);
		return geneNameString;
	}

	private String getProteinNameString(PCQProteinNode proteinNode) {
		String accString = proteinNode.getKey();
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
			Entry entry = annotatedProteins.get(acc);
			if (entry != null && entry.getName() != null && !entry.getName().isEmpty()) {
				final String proteinName = entry.getName().get(0);
				if (proteinName.contains("obsolete")) {
					return acc;
				}
				return proteinName;
			}
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
				if (!taxonomies.isEmpty()) {
					return colorManager.getColorForProteinTaxonomy(taxonomies.iterator().next());
				}
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
		ret.setOutline(ColorGenerator.getHexString(outlineColor));
		ret.setFill(ColorGenerator.getHexString(fillColor));
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
		ret.getAtt().add(createNodeGraphicAtt("NODE_LABEL_COLOR", ColorGenerator.getHexString(labelColor), STRING));
		ret.getAtt().add(createNodeGraphicAtt("NODE_VISIBLE", "true", STRING));
		ret.getAtt().add(createNodeGraphicAtt("NODE_DEPTH", "0.0", STRING));
		ret.getAtt().add(createNodeGraphicAtt("NODE_BORDER_TRANSPARENCY", "255", STRING));
		ret.getAtt().add(createNodeGraphicAtt("NODE_LABEL_FONT_FACE", "Dialog,plain,12", STRING));
		ret.getAtt().add(createNodeGraphicAtt("NODE_LABEL_TRANSPARENCY", "255", STRING));

		return ret;
	}

	private Map<String, AttributeValueType> getAttributesForEdge(String edgeName, QuantRatio sharedPepRatio,
			NWResult alignmentResult, Collection<Classification2Case> classification2Cases) {
		Map<String, AttributeValueType> ret = new THashMap<String, AttributeValueType>();
		// ret.put("name", edgeName);
		ret.put(PCQ_ID, new AttributeValueType(edgeName));
		if (sharedPepRatio != null) {
			final Double log2CountRatio = sharedPepRatio.getLog2Ratio(cond1, cond2);
			if (log2CountRatio != null) {
				ret.put(COUNT_RATIO, new AttributeValueType(log2CountRatio));
			}
		}

		// StringBuilder classification1String = new StringBuilder();
		// if (classification1Cases != null) {
		// for (Classification1Case classification1Case : classification1Cases)
		// {
		// final String value = String.valueOf(classification1Case.getCaseID());
		// if (!"".equals(classification1String.toString())) {
		// classification1String.append(", ");
		// }
		// classification1String.append(value);
		// }
		// ret.put("Classification1Case", new
		// AttributeValueType(classification1String.toString(),
		// AttType.string));
		// }

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
		edge.setFill(ColorGenerator.getHexString(fillColor));
		edge.setWidth(3);
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_SELECTED", "false", STRING));
		edge.getAtt()
				.add(createEdgeGraphAttribute("EDGE_LABEL_COLOR", ColorGenerator.getHexString(Color.black), STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_TARGET_ARROW_SHAPE", "none", STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_SOURCE_ARROW_UNSELECTED_PAINT",
				ColorGenerator.getHexString(Color.black), STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_TARGET_ARROW_SELECTED_PAINT", "#FFFF00", STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_LABEL_TRANSPARENCY", "255", STRING));
		edge.getAtt().add(
				createEdgeGraphAttribute("EDGE_STROKE_SELECTED_PAINT", ColorGenerator.getHexString(Color.red), STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_LINE_TYPE", "SOLID", STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_TOOLTIP", tooltip, STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_CURVED", "true", STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_TRANSPARENCY", "255", STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_BEND", "", STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_LABEL_FONT_FACE", "Dialog,plain,10", STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_LABEL", label, STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_SOURCE_ARROW_SHAPE", "none", STRING));
		edge.getAtt().add(createEdgeGraphAttribute("EDGE_TARGET_ARROW_UNSELECTED_PAINT",
				ColorGenerator.getHexString(Color.black), STRING));
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
		Set<String> accs = new THashSet<String>();
		if (accession.contains(" ")) {
			final String[] split = accession.split("\\s+");
			for (String string : split) {
				accs.add(string);
			}
		} else {
			accs.add(accession);
		}
		Map<String, Entry> ret = new THashMap<String, Entry>();
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
			if (!significantlyRegulatedProteinClusters.isEmpty()) {
				String fdrText = "";
				if (params.isPerformRatioIntegration() && params.getSignificantFDRThreshold() != null) {
					fdrText = params.getSignificantFDRThreshold() + "_";
				}
				final String fileName = outputFileFolder.getAbsolutePath() + File.separator + outputPrefix
						+ "_cytoscape_Significants_" + fdrText + outputSuffix + ".xgmml";

				File xgmmlOutPutFileFDR = new File(fileName);
				exportToGmmlFromProteinClustersUsingNodes(xgmmlOutPutFileFDR,
						outputPrefix + "_FDR" + fdrThreshold + "_" + outputSuffix,
						significantlyRegulatedProteinClusters, condition1, condition2, colorManager);
			}

			if (params.isApplyClassificationsByProteinPair()) {
				if (params.isCollapseIndistinguishablePeptides() && params.isCollapseIndistinguishableProteins()) {
					// classification 2
					final Map<Classification2Case, Set<ProteinCluster>> proteinPairsByClassification2 = getProteinClustersByClassification2(
							clusterCollection);
					final Classification2Case[] cases2 = Classification2Case.values();
					for (Classification2Case case2 : cases2) {
						if (case2 == Classification2Case.CASE6) {
							// skip
							continue;
						}
						if (proteinPairsByClassification2.containsKey(case2)) {
							log.info("Creating XGMML for case " + case2.getCaseID() + "(" + case2.getExplanation()
									+ ")...");
							File xgmmlOutPutFile2 = new File(
									outputFileFolder.getAbsolutePath() + File.separator + outputPrefix + "_cytoscape_"
											+ case2.getCaseID() + "-" + case2.name() + "_" + outputSuffix + ".xgmml");
							exportToGmmlFromProteinClustersUsingNodes(xgmmlOutPutFile2,
									outputPrefix + "_" + case2.getCaseID() + "-" + case2.name() + "_" + outputSuffix,
									proteinPairsByClassification2.get(case2), condition1, condition2, colorManager);
						}
					}
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
		Set<ProteinCluster> ret = new THashSet<ProteinCluster>();
		for (ProteinCluster proteinCluster : clusters) {
			final Set<PCQPeptideNode> peptideNodes = proteinCluster.getPeptideNodes();
			for (PCQPeptideNode peptideNode : peptideNodes) {
				final QuantRatio consensusRatio = peptideNode.getConsensusRatio(cond1, cond2);
				if (consensusRatio != null) {
					final Score score = consensusRatio.getAssociatedConfidenceScore();
					if (score != null) {
						try {
							if (score.getScoreName().equals(PCQUtils.FDR_CONFIDENCE_SCORE_NAME)
									&& score.getValue() != null) {
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
		Map<Classification1Case, Set<ProteinCluster>> ret = new THashMap<Classification1Case, Set<ProteinCluster>>();
		for (ProteinCluster proteinCluster : clusters) {
			final Set<ProteinPair> proteinPairs = proteinCluster.getProteinPairs();
			for (ProteinPair proteinPair : proteinPairs) {
				final Collection<Classification1Case> classification1PairCases = proteinPair.getClassification1Case()
						.values();
				for (Classification1Case classification1PairCase : classification1PairCases) {
					if (ret.containsKey(classification1PairCase)) {
						ret.get(classification1PairCase).add(proteinCluster);
					} else {
						Set<ProteinCluster> set = new THashSet<ProteinCluster>();
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
		Map<Classification2Case, Set<ProteinCluster>> ret = new THashMap<Classification2Case, Set<ProteinCluster>>();
		for (ProteinCluster proteinCluster : clusters) {
			final Set<ProteinPair> proteinPairs = proteinCluster.getProteinPairs();
			for (ProteinPair proteinPair : proteinPairs) {
				final Collection<Classification2Case> classification2PairCases = proteinPair.getClassification2Cases()
						.values();
				for (Classification2Case classification2PairCase : classification2PairCases) {
					if (ret.containsKey(classification2PairCase)) {
						ret.get(classification2PairCase).add(proteinCluster);
					} else {
						Set<ProteinCluster> set = new THashSet<ProteinCluster>();
						set.add(proteinCluster);
						ret.put(classification2PairCase, set);
					}
				}
			}
		}
		return ret;
	}

}
