#include "PeptideCoordinates.h"

PeptideCoordinates::PeptideCoordinates(void) : m_coordinates(), m_transcript_coordinates(), m_cds_annotation_correct(0) {}

PeptideCoordinates::PeptideCoordinates(CoordinateMapType const&coordinates, unsigned int CDSannotationcorrect) :
	m_coordinates(coordinates), m_cds_annotation_correct(CDSannotationcorrect) {
	get_exon_coordinates();
	m_transcript_coordinates = get_transcript_coordinates();
}
PeptideCoordinates::~PeptideCoordinates(void) {}

std::vector<GenomeCoordinates> PeptideCoordinates::find_coordinates(std::pair<unsigned int, unsigned int> peptideproteincoords) {
	Coordinates ptm_coordinates;
	ptm_coordinates.start = peptideproteincoords.first;
	ptm_coordinates.end = peptideproteincoords.second;
	ptm_coordinates.Cterm = off3;
	ptm_coordinates.Nterm = off3;

	//equal range returns 2 iterators: the first and the last occurence of ptm_coordinates in the sorted multimap (m_coordinates)
	std::pair<CoordinateMapType::iterator, CoordinateMapType::iterator> ptm_search_it = m_coordinates.equal_range(ptm_coordinates);

	std::vector<GenomeCoordinates> ptm;

	for (CoordinateMapType::iterator it = ptm_search_it.first; it != ptm_search_it.second; ++it) {
		//for each found ptm_coordinates entry the genomic coordinates of that PTM are computed
		CoordinateMapType::value_type coordinates_partial = get_coordinates(it->first, it->second, ptm_coordinates);
		ptm.push_back(coordinates_partial.second);
	}
	//sorted
	std::sort(ptm.begin(), ptm.end(), compare_coordinates_ascending);
	return ptm;
}

std::vector<GenomeCoordinates> PeptideCoordinates::get_exon_coordinates() {
	//used to generate the exon information in the output
	//this is mainly important to find reverse exons.
	if (m_coordinate_vector.empty()) {
		for (CoordinateMapType::iterator it = m_coordinates.begin(); it != m_coordinates.end(); ++it) {
			GenomeCoordinates coord = it->second;
			if (coord.strand == rev) {
				unsigned int start = coord.end;
				unsigned int end = coord.start;
				coord.start = start;
				coord.end = end;
			}
			m_coordinate_vector.push_back(coord);
		}
		std::sort(m_coordinate_vector.begin(), m_coordinate_vector.end(), compare_coordinates_ascending);
	}
	return m_coordinate_vector;
}

GenomeCoordinates PeptideCoordinates::get_transcript_coordinates() {
	GenomeCoordinates transcript_coordinates;
	for (size_t i(0); i < m_coordinate_vector.size(); ++i) {
		if (i == 0) {
			transcript_coordinates = m_coordinate_vector.at(i);
		} else {
			if (transcript_coordinates.start > m_coordinate_vector.at(i).start) {
				transcript_coordinates.start = m_coordinate_vector.at(i).start;
			}
			if (transcript_coordinates.end < m_coordinate_vector.at(i).end) {
				transcript_coordinates.end = m_coordinate_vector.at(i).end;
			}
		}
	}

	if (m_cds_annotation_correct != 0) {
		if (transcript_coordinates.strand == fwd) {
			transcript_coordinates.start = transcript_coordinates.start - m_cds_annotation_correct;
			transcript_coordinates.end = transcript_coordinates.end - m_cds_annotation_correct;
		} else if (transcript_coordinates.strand == rev) {
			transcript_coordinates.end = transcript_coordinates.end + m_cds_annotation_correct;
			transcript_coordinates.start = transcript_coordinates.start + m_cds_annotation_correct;
		}
	}
	return transcript_coordinates;
}

bool PeptideCoordinates::operator<(const PeptideCoordinates& rhs) const {
	if (compare_coordinates_ascending_whole(m_transcript_coordinates, rhs.m_transcript_coordinates)) {
		return true;
	} else if (compare_coordinates_ascending_whole(rhs.m_transcript_coordinates, m_transcript_coordinates)) {
		return false;
	}
	return compare_genome_coordinate_sets_ascending(m_coordinate_vector, rhs.m_coordinate_vector);

}

bool peptidecoords_p_compare::operator()(const PeptideCoordinates* lhs, const PeptideCoordinates* rhs) const {
	return (*lhs) < (*rhs);
}

bool peptidecoords_pair_p_compare::operator()(const std::pair<PeptideCoordinates*, GenomeCoordinates>& lhs, const std::pair<PeptideCoordinates*, GenomeCoordinates>& rhs) const {
	return (*(lhs.first)) < (*(rhs.first));
}
