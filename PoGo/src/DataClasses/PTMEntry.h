#ifndef PTMENTRY_H
#define PTMENTRY_H

//#include "GeneEntry.h"

//the PTMEntry contains information about post translational modifications like methlation.
//and at which position they appear.
class PTMEntry {
public:
	//ctr /dtr
	PTMEntry();
	PTMEntry(std::string const& name, size_t start, size_t end);
	PTMEntry(PTMEntry const& rhs);
	~PTMEntry();
	//ctr /dtr 

	//returns true if the ptm_psi_name, start and end coords are the same, 
	//otherwise returns false.
	bool operator==(PTMEntry const& rhs) const;
	//returns a pair of start and end coord for the current PTM.
	std::pair<size_t, size_t> get_range() const;
	//decreases m_peptide_start coord if coord is smaller than that 
	//or increases m_peptide start coor if coord is larger than that.
	//sets both to coord if they havent been set before.
	void add_coord(size_t coord);
	//returns a vector of all found genome coordinates.
	std::vector<std::pair<PeptideCoordinates*, GenomeCoordinates>> get_genome_coordinates();
	//adds genome coordinates.
	void add_genome_coordinates(PeptideCoordinates* coords, GenomeCoordinates ptmcoords);

private:
	//holds the name of the modification
	std::string m_ptm_psi_name;
	//holds the lowest start coord
	size_t m_peptide_start_coord;
	//holds the highest end coord.
	size_t m_peptide_end_coord;
	//holds all added ptms and their genomic locations.
	std::set<std::pair<PeptideCoordinates*, GenomeCoordinates>, peptidecoords_pair_p_compare> m_ptm_coord;


};

#endif 
