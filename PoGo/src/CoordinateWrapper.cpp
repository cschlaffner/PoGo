#include "CoordinateWrapper.h"
#include "KmerMap.h"
#include "ExistingPeptides.h"

CoordinateWrapper::CoordinateWrapper(void) :
	m_map(std::map<std::string, ProteinEntry>()),
	m_existing_peptides(new ExistingPeptides()) {}
CoordinateWrapper::~CoordinateWrapper(void) {}

int CoordinateWrapper::size() const {
	return m_map.size();
}

ProteinEntry& CoordinateWrapper::lookup_entry(std::string transcriptId) {
	return m_map[transcriptId];
}

void CoordinateWrapper::add(ProteinEntry entry) {
	m_map[entry.get_transcript_id()] = entry;
}

void CoordinateWrapper::read_fasta_file(std::string file) {
	if (!FastaParser::get_instance()->open(file)) {
		throw FastaParser__file_not_found_exception("");
	}

	FastaEntry fasta_entry;
	while (!(fasta_entry = FastaParser::get_instance()->next_entry()).is_empty()) {
		add(ProteinEntry(fasta_entry));
	}

	FastaParser::get_instance()->close();
}

void CoordinateWrapper::add_all_proteins_to_kmer_map(KmerMap& kmerMap) {
	for (std::map<std::string, ProteinEntry>::iterator it = m_map.begin(); it != m_map.end(); ++it) {
		kmerMap.add_protein(it->second);
	}
}

void CoordinateWrapper::add_to_existing_peptides(std::string const& peptideSequence, PeptideEntry* peptideEntry) const {
	m_existing_peptides->add(peptideSequence, peptideEntry);
}

std::vector<PeptideEntry*>& CoordinateWrapper::get_existing_peptides_at(std::string const& peptideSequence) const {
	return (*m_existing_peptides)[peptideSequence];
}

bool CoordinateWrapper::peptide_already_exists(std::string const& peptideSequence) const {
	return m_existing_peptides->contains(peptideSequence);
}