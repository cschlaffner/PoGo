#include "MapEntry.h"

MapEntry::MapEntry(GeneEntry* geneentry_p) :
	m_p_gene_entry(geneentry_p),
	m_peptide_entries(std::map<std::string, PeptideEntry*>()),
	m_transcripts(std::set<std::string>()) {}
MapEntry::~MapEntry(void) {}

void MapEntry::add_transcript_id(const std::string& transcriptID) {
	m_transcripts.insert(transcriptID);
}

bool MapEntry::operator<(const MapEntry& rhs) const {
	return (*m_p_gene_entry) < (*(rhs.m_p_gene_entry));
}

std::ostream& MapEntry::to_gtf(const std::string& source, std::ostream& os, bool chrincluded) {
	if (m_peptide_entries.size() > 0) {
		m_p_gene_entry->to_gtf(source, os) << "\n";

		std::set<PeptideEntry*, peptideentry_p_compare> peptide_entries_set = std::set<PeptideEntry*, peptideentry_p_compare>();
		for (std::map<std::string, PeptideEntry*>::iterator it = m_peptide_entries.begin(); it != m_peptide_entries.end(); ++it) {
			peptide_entries_set.insert(it->second);
		}

		for (std::set<PeptideEntry*, peptideentry_p_compare>::iterator it = peptide_entries_set.begin(); it != peptide_entries_set.end(); ++it) {
			(*it)->to_gtf(source, os) << "\n";
		}
	}
	return os;
}

std::ostream& MapEntry::to_bed(std::ostream& os, bool chrincluded) {
	if (m_peptide_entries.size() > 0) {
		std::set<PeptideEntry*, peptideentry_p_compare> peptide_entries_set = std::set<PeptideEntry*, peptideentry_p_compare>();
		for (std::map<std::string, PeptideEntry*>::iterator it = m_peptide_entries.begin(); it != m_peptide_entries.end(); ++it) {
			peptide_entries_set.insert(it->second);
		}

		for (std::set<PeptideEntry*, peptideentry_p_compare>::iterator it = peptide_entries_set.begin(); it != peptide_entries_set.end(); ++it) {
			(*it)->to_bed(os);
		}
	}
	return os;
}

std::ostream& MapEntry::to_gct(std::vector<std::string> const& tissuevector, std::ostream& os, bool chrincluded) {
	if (m_peptide_entries.size() > 0) {
		std::set<PeptideEntry*, peptideentry_p_compare> peptide_entries_set = std::set<PeptideEntry*, peptideentry_p_compare>();
		for (std::map<std::string, PeptideEntry*>::iterator it = m_peptide_entries.begin(); it != m_peptide_entries.end(); ++it) {
			peptide_entries_set.insert(it->second);
		}

		for (std::set<PeptideEntry*, peptideentry_p_compare>::iterator it = peptide_entries_set.begin(); it != peptide_entries_set.end(); ++it) {
			(*it)->to_gct(m_p_gene_entry->get_id(), tissuevector, os);
		}
	}
	return os;
}

std::ostream& MapEntry::to_ptmbed(std::ostream& os, std::ostream& os2, bool chrincluded) {
	if (m_peptide_entries.size() > 0) {
		std::set<PeptideEntry*, peptideentry_p_compare> peptide_entries_set = std::set<PeptideEntry*, peptideentry_p_compare>();
		for (std::map<std::string, PeptideEntry*>::iterator it = m_peptide_entries.begin(); it != m_peptide_entries.end(); ++it) {
			peptide_entries_set.insert(it->second);
		}

		for (std::set<PeptideEntry*, peptideentry_p_compare>::iterator it = peptide_entries_set.begin(); it != peptide_entries_set.end(); ++it) {
			(*it)->to_ptmbed(os);
			if ((*it)->noPTM()) {
				(*it)->to_bed(os2, true);
			}
		}
	}
	return os;
}

void MapEntry::remove_peptides() {
	m_peptide_entries.clear();
}

unsigned int MapEntry::add_peptide(CoordinateWrapper& coordwrapper, const std::string& sequence, const std::string& tag, unsigned int sigPSMs, unsigned int genes, std::ofstream& ofstream, double quant, gene_id_map_t::value_type const& transcriptsEntry) {
	unsigned int newly_added = 0;
	std::string sequence_wo_ptm = remove_ptms(sequence);

	if (transcriptsEntry.second.m_entries.size() > 0) {
		if (m_peptide_entries.count(sequence_wo_ptm) != 1) {
			PeptideEntry* new_peptide = new PeptideEntry(m_p_gene_entry);
			m_peptide_entries.insert(std::pair<std::string, PeptideEntry*>(sequence_wo_ptm, new_peptide));
			coordwrapper.add_to_existing_peptides(sequence_wo_ptm, new_peptide);
			newly_added = 1;
		}
		m_peptide_entries[sequence_wo_ptm]->add_peptide(coordwrapper, sequence_wo_ptm, sequence, tag, sigPSMs, transcriptsEntry.second, genes, ofstream, quant);
	}
	return newly_added;
}

bool mapentry_p_compare::operator()(const MapEntry* lhs, const MapEntry* rhs) const {
	return (*lhs) < (*rhs);
}