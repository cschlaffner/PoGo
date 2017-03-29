#ifndef PEPTIDECOORDINATES_H
#define PEPTIDECOORDINATES_H

#include "CoordinateWrapper.h"


//in this class all the information will be used to find the 
//genomic location of a peptide.
class PeptideCoordinates {
public:
	//ctr / dtr
	PeptideCoordinates(void);
	PeptideCoordinates(CoordinateMapType const& coordinates, unsigned int CDSannotationcorrect);
	~PeptideCoordinates(void);

	//generates exon coordinates from the coordinate map type.
	std::vector<GenomeCoordinates>	get_exon_coordinates();
	//generates GenomeCoordinates for the associated transcript.
	GenomeCoordinates				get_transcript_coordinates();
	//this function is the actual search function that will find genomic coordinates for a peptide.
	std::vector<GenomeCoordinates> find_coordinates(std::pair<unsigned int, unsigned int> peptideproteincoords);

	//lesser than operator. returns compare_coordinates_ascending_whole (this, other)
	//otherwise returns false if compare_coordinates_ascending_whole (other, this)
	//otherwise returns compare_genome_coordinate_sets_ascending(this, other)
	bool operator<(const PeptideCoordinates& rhs) const;

	std::set<std::string> get_trasncript_ids();
	std::set<std::string> get_exon_ids();

private:
	//holds the found coordinates.
	std::vector<GenomeCoordinates>	m_coordinate_vector;
	//passed to the constuctor. holds the peptides coordinates and the associated genomic coordinates. 
	//peptide coordinates sorted; used for computing ptm genomic coordinates.
	CoordinateMapType				m_coordinates;
	//the coordinates of the current transcript.
	GenomeCoordinates				m_transcript_coordinates;
	//tests if the annotation of the coding sequence is dividable by 3bp and not offset 
	//due to incomplete transcript annotation.
	unsigned int					m_cds_annotation_correct;

	std::set<std::string> m_transcriptids;
	std::set<std::string> m_exonids;

	void add_ids();
};

//allows for the comparision of two peptide coordinate pointers.
//it dereferences them and calls the operator<
struct peptidecoords_p_compare {
	bool operator()(const PeptideCoordinates* lhs, const PeptideCoordinates* rhs) const;
};

//allows for the comparision of std::pair<PeptideCoordinates*, GenomeCoordinates> pointers. 
//dereferences and compares the std::pair::first (PeptideCoordinates::operator<)
struct peptidecoords_pair_p_compare {
	bool operator()(const std::pair<PeptideCoordinates*, GenomeCoordinates>& lhs, const std::pair<PeptideCoordinates*, GenomeCoordinates>& rhs) const;
};

#endif
