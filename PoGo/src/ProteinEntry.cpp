#include "ProteinEntry.h"
#include "PeptideMapperUtils.h"

ProteinEntry::ProteinEntry(void) :
	m_fasta_header(),
	m_transcript_id(),
	m_gene_id(),
	m_aa_sequence(),
	m_coordinates_map(CoordinateMapType()),
	m_cds_annotation_correct(0) {}
ProteinEntry::ProteinEntry(std::string const&fastaHeader, std::string const&AAsequence) {
	init(fastaHeader, AAsequence);
}
ProteinEntry::ProteinEntry(FastaEntry const&fastaEntry) {
	init(fastaEntry);
}
ProteinEntry::~ProteinEntry(void) {}

void ProteinEntry::init(std::string fastaHeader, std::string AAsequence) {
	if (fastaHeader.substr(0, 1) == ">") {
		m_fasta_header = fastaHeader;
		m_transcript_id = extract_transcript_id_fasta(fastaHeader);
		m_gene_id = extract_gene_id_fasta(fastaHeader);
		m_aa_sequence = AAsequence;
		m_coordinates_map = CoordinateMapType();
		m_cds_annotation_correct = 0;
	}
}

std::string ProteinEntry::extract_transcript_id_fasta(std::string str) {
	size_t index = str.find(GENOME_MAPPER_GLOBALS::ID::FASTA_TRANSCRIPT_ID) + GENOME_MAPPER_GLOBALS::ID::FASTA_TRANSCRIPT_ID.length();
	std::string value("");

	if (index != std::string::npos) {
		while (str[index] != ' ') {
			value = value + str[index];
			index += 1;
		}
	}

	return value;
}

std::string ProteinEntry::extract_gene_id_fasta(std::string str) {
	size_t index = str.find(GENOME_MAPPER_GLOBALS::ID::FASTA_GENE_ID) + GENOME_MAPPER_GLOBALS::ID::FASTA_GENE_ID.length();
	std::string value("");
	if (index != std::string::npos) {
		while (str[index] != ' ') {
			value = value + str[index];
			index += 1;
		}
	} 
	return value;
}

void ProteinEntry::set_coordinate_map(CoordinateMapType coordinatesMap) {
	m_coordinates_map = coordinatesMap;
}

void ProteinEntry::init(FastaEntry const&fastaEntry) {
	init(fastaEntry.get_header(), fastaEntry.get_sequence());
}

std::string const& ProteinEntry::get_transcript_id() const {
	return m_transcript_id;
}

std::string const& ProteinEntry::get_gene_id() const {
	return m_gene_id;
}

std::string const& ProteinEntry::get_sequence() const {
	return m_aa_sequence;
}

unsigned int ProteinEntry::is_cds_annotation_correct() const {
	return m_cds_annotation_correct;
}

std::vector<std::vector<GenomeCoordinates>> ProteinEntry::find_coordinates(int peptideseqSize, std::vector<position_mismatch_t> const& positions) {
	std::vector<std::vector<GenomeCoordinates> > found_coordinates;
	Coordinates peptide_coordinates;
	peptide_coordinates.Cterm = off3;
	peptide_coordinates.Nterm = off3;

	//iterate all found positions
	for (std::vector<position_mismatch_t>::const_iterator it = positions.begin(); it != positions.end(); ++it) {
		peptide_coordinates.start = (*it).position_in_protein();
		peptide_coordinates.end = (*it).position_in_protein() + (peptideseqSize - 1);
		//equal range will return two iterators: one to the first occurence of peptide_coordinates, and one to the last occurence.
		//m_coordinates map is a multimap, so elements can appear more than one (hence the range) and will always be in sorted order (which they have to be for equal range to work)
		std::pair<CoordinateMapType::iterator, CoordinateMapType::iterator> search = m_coordinates_map.equal_range(peptide_coordinates);
		std::vector<GenomeCoordinates> single;
		for (CoordinateMapType::iterator it_inner = search.first; it_inner != search.second; ++it_inner) {
			//the get_coordinates function will generate a CoordinateMapType entry that holds the proteinCoordinates and the GenomeCoordinates for a specific peptide.
			CoordinateMapType::value_type coordinatesPartial = get_coordinates(it_inner->first, it_inner->second, peptide_coordinates);
			//this has to be done several times to find all genomic & protein coordinates of a peptide
			single.push_back(coordinatesPartial.second);
		}
		//and has to be done several times to find all peptides.
		found_coordinates.push_back(single);
	}
	return found_coordinates;
}




