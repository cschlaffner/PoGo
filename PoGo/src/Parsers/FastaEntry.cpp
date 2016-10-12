//#include "FastaEntry.h"

FastaEntry::FastaEntry(std::string const& header, std::string const&AAsequence)
	: m_header(header), m_aa_sequence(AAsequence) {}
FastaEntry::FastaEntry() : m_header(), m_aa_sequence() {}
FastaEntry::~FastaEntry(void) {}

std::string FastaEntry::get_header() const {
	return m_header;
}

std::string FastaEntry::get_sequence() const {
	return m_aa_sequence;
}

bool FastaEntry::is_empty() const {
	return m_header.empty() && m_aa_sequence.empty();
}
