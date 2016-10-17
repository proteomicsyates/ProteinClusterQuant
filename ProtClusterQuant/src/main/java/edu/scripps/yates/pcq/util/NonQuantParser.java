package edu.scripps.yates.pcq.util;

import java.io.IOException;
import java.util.HashMap;
import java.util.Set;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.census.analysis.util.KeyUtils;
import edu.scripps.yates.census.read.AbstractQuantParser;
import edu.scripps.yates.census.read.model.QuantStaticMaps;
import edu.scripps.yates.census.read.model.QuantifiedPeptide;
import edu.scripps.yates.census.read.model.QuantifiedProteinFromDBIndexEntry;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.dbindex.IndexedProtein;
import edu.scripps.yates.dbindex.util.PeptideNotFoundInDBIndexException;
import edu.scripps.yates.dtaselectparser.DTASelectParser;
import edu.scripps.yates.dtaselectparser.util.DTASelectPSM;
import edu.scripps.yates.dtaselectparser.util.DTASelectProtein;
import edu.scripps.yates.pcq.model.NonQuantifiedPSM;
import edu.scripps.yates.pcq.model.NonQuantifiedProtein;
import edu.scripps.yates.utilities.fasta.FastaParser;

public class NonQuantParser extends AbstractQuantParser {
	private final static Logger log = Logger.getLogger(NonQuantParser.class);
	private final DTASelectParser dtaSelectParser;

	public NonQuantParser(DTASelectParser dtaSelectParser) {
		this.dtaSelectParser = dtaSelectParser;
	}

	@Override
	protected void process() {
		try {
			final HashMap<String, DTASelectPSM> psms = dtaSelectParser.getDTASelectPSMsByPSMID();
			for (DTASelectPSM psm : psms.values()) {
				processPSM(psm);

			}
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private void processPSM(DTASelectPSM psm) throws IOException {

		final String inputFileName = psm.getMsRunId();
		String experimentKey = FilenameUtils.getBaseName(inputFileName);

		QuantifiedPSMInterface quantifiedPSM = new NonQuantifiedPSM(psm, chargeStateSensible);
		if (quantifiedPSM.getSequence().equals("AGALQVSAMSMVNGFFSFSAALELKERHAK")) {
			log.info(quantifiedPSM);
		}
		quantifiedPSM.addFileName(psm.getRunID());
		final String psmKey = KeyUtils.getSpectrumKey(quantifiedPSM, chargeStateSensible);
		// in case of TMT, the psm may have been created before
		if (QuantStaticMaps.psmMap.containsKey(psmKey)) {
			quantifiedPSM = QuantStaticMaps.psmMap.getItem(psmKey);
		}
		QuantStaticMaps.psmMap.addItem(quantifiedPSM);

		// psms.add(quantifiedPSM);
		// add to map
		if (!localPsmMap.containsKey(quantifiedPSM.getKey())) {
			localPsmMap.put(quantifiedPSM.getKey(), quantifiedPSM);
		}

		// create the peptide
		QuantifiedPeptideInterface quantifiedPeptide = null;
		final String peptideKey = KeyUtils.getSequenceKey(quantifiedPSM, true);
		if (peptideKey.equals("EHALAQAELLK")) {
			log.info(quantifiedPSM);
		}
		if (QuantStaticMaps.peptideMap.containsKey(peptideKey)) {
			quantifiedPeptide = QuantStaticMaps.peptideMap.getItem(peptideKey);
		} else {
			if (quantifiedPSM.getFullSequence().equals("PNS(114.042927)VPQE(14.01565)LAATTEKTEPNSQEDKNDGGK")) {
				log.info(quantifiedPSM);
			}
			quantifiedPeptide = new QuantifiedPeptide(quantifiedPSM);
			quantifiedPeptide.addFileName(inputFileName);
		}
		QuantStaticMaps.peptideMap.addItem(quantifiedPeptide);

		quantifiedPSM.setQuantifiedPeptide(quantifiedPeptide, true);
		// add peptide to map
		if (!localPeptideMap.containsKey(peptideKey)) {
			localPeptideMap.put(peptideKey, quantifiedPeptide);
		}

		if (dbIndex != null) {
			String cleanSeq = quantifiedPSM.getSequence();
			final Set<IndexedProtein> indexedProteins = dbIndex.getProteins(cleanSeq);
			if (indexedProteins.isEmpty()) {
				if (!super.ignoreNotFoundPeptidesInDB) {
					throw new PeptideNotFoundInDBIndexException("The peptide " + cleanSeq
							+ " is not found in Fasta DB.\nReview the default indexing parameters such as the number of allowed misscleavages.");
				}
				// log.warn("The peptide " + cleanSeq +
				// " is not found in Fasta DB.");
				// continue;
			}
			// create a new Quantified Protein for each
			// indexedProtein
			for (IndexedProtein indexedProtein : indexedProteins) {
				String proteinKey = KeyUtils.getProteinKey(indexedProtein);
				if (proteinKey.equals("Q03181")) {
					log.info(indexedProtein);
				}
				QuantifiedProteinInterface quantifiedProtein = null;
				if (QuantStaticMaps.proteinMap.containsKey(proteinKey)) {
					quantifiedProtein = QuantStaticMaps.proteinMap.getItem(proteinKey);
				} else {
					quantifiedProtein = new QuantifiedProteinFromDBIndexEntry(indexedProtein);
					quantifiedProtein.addFileName(inputFileName);
				}
				// add psm to the proteins
				quantifiedProtein.addPSM(quantifiedPSM, true);
				// add protein to the psm
				quantifiedPSM.addQuantifiedProtein(quantifiedProtein, true);
				// add peptide to the protein
				quantifiedProtein.addPeptide(quantifiedPeptide, true);
				// add to the map (if it was already there
				// is not a problem, it will be only once)
				addToMap(proteinKey, proteinToPeptidesMap, KeyUtils.getSequenceKey(quantifiedPSM, true));
				// add protein to protein map
				localProteinMap.put(proteinKey, quantifiedProtein);
				// add to protein-experiment map
				addToMap(experimentKey, experimentToProteinsMap, proteinKey);

			}
		}
		final Set<DTASelectProtein> proteins = psm.getProteins();
		for (DTASelectProtein dtaSelectProtein : proteins) {
			String proteinKey = FastaParser.getACC(dtaSelectProtein.getLocus()).getFirstelement();
			QuantifiedProteinInterface quantifiedProtein = null;
			if (QuantStaticMaps.proteinMap.containsKey(proteinKey)) {
				quantifiedProtein = QuantStaticMaps.proteinMap.getItem(proteinKey);
			} else {
				quantifiedProtein = new NonQuantifiedProtein(dtaSelectProtein);
				quantifiedProtein.addFileName(inputFileName);
			}
			// add psm to the proteins
			quantifiedProtein.addPSM(quantifiedPSM, true);
			// add protein to the psm
			quantifiedPSM.addQuantifiedProtein(quantifiedProtein, true);
			// add peptide to the protein
			quantifiedProtein.addPeptide(quantifiedPeptide, true);
			// add to the map (if it was already there
			// is not a problem, it will be only once)
			addToMap(proteinKey, proteinToPeptidesMap, KeyUtils.getSequenceKey(quantifiedPSM, true));
			// add protein to protein map
			localProteinMap.put(proteinKey, quantifiedProtein);
			// add to protein-experiment map
			addToMap(experimentKey, experimentToProteinsMap, proteinKey);

		}
	}
}
