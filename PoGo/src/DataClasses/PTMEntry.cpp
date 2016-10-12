//#include "PTMEntry.h"

PTMEntry::PTMEntry() :
	m_peptide_start_coord(std::string::npos), m_peptide_end_coord(std::string::npos) {}
PTMEntry::PTMEntry(std::string const& name, size_t start, size_t end)
	: m_ptm_psi_name(name), m_peptide_start_coord(start), m_peptide_end_coord(end) {}
PTMEntry::PTMEntry(PTMEntry const& rhs) :
	m_ptm_psi_name(rhs.m_ptm_psi_name),
	m_peptide_start_coord(rhs.m_peptide_start_coord),
	m_peptide_end_coord(rhs.m_peptide_end_coord),
	m_ptm_coord(rhs.m_ptm_coord) {}
PTMEntry::~PTMEntry() {}

bool PTMEntry::operator==(PTMEntry const& rhs) const {
	return m_ptm_psi_name.compare(rhs.m_ptm_psi_name) == 0 && m_peptide_start_coord == rhs.m_peptide_start_coord && m_peptide_end_coord == rhs.m_peptide_end_coord;
}

std::pair<size_t, size_t> PTMEntry::get_range() const {
	return std::pair<size_t, size_t>(m_peptide_start_coord, m_peptide_end_coord);
}

void PTMEntry::add_coord(size_t coord) {
	if (m_peptide_start_coord == std::string::npos && m_peptide_end_coord == std::string::npos) {
		m_peptide_start_coord = coord;
		m_peptide_end_coord = coord;
	} else if (coord <= m_peptide_start_coord) {
		m_peptide_start_coord = coord;
	} else if (coord >= m_peptide_end_coord) {
		m_peptide_end_coord = coord;
	}
}

std::vector<std::pair<PeptideCoordinates*, GenomeCoordinates>> PTMEntry::get_genome_coordinates() {
	std::vector<std::pair<PeptideCoordinates*, GenomeCoordinates>> ptm_coord;
	for (std::set<std::pair<PeptideCoordinates*, GenomeCoordinates>>::iterator it = m_ptm_coord.begin(); it != m_ptm_coord.end(); ++it) {
		ptm_coord.push_back(*it);
	}
	return ptm_coord;
}

void PTMEntry::add_genome_coordinates(PeptideCoordinates* coords, GenomeCoordinates ptmcoords) {
	m_ptm_coord.insert(std::pair<PeptideCoordinates*, GenomeCoordinates>(coords, ptmcoords));
}
