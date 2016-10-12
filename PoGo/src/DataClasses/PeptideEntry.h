#ifndef PEPTIDEENTRY_H
#define PEPTIDEENTRY_H

//#include "PTMEntry.h"

struct transcripts_t;

//a peptide entry contains information about a single peptide:
//all its occurences in the genome, 
//the associated gene, if it is unique to transcript or gene, 
//the tissues it appears in, 
//as well as the first and last appearance in the genome.
class PeptideEntry {
public:
	//ctr &dtr
	explicit PeptideEntry(GeneEntry* associatedgene);
	~PeptideEntry(void);
	//end ctr&dtr
	
	//adds a peptide. this funciton is used if this specific peptide has not yet been found. it will also find the peptides genomic coordinates.
	void add_peptide(CoordinateWrapper& coordwrapper, const std::string& sequence, const std::string& ptmsequence, const std::string& tag, unsigned int sigPSMs, transcripts_t const& transcripts, unsigned int genes, std::ofstream& ofstream, double quant);
	//adds a peptide. this function is used if a peptide has appeared before.
	void add_peptide(std::string const& ptmSequence, std::string const& tag, unsigned int sigPSMs, double quant);

	//comparator to check if one PeptideEntry is smaller than another one.
	//returns true if the startcoordinate are smaller than rhs' startcoordinate
	//otherwise returns true if endcoordinate is smaller than rhs' endcoordinate
	//otherwise returns true if lhs' sequence is lexicographically lesser than rhs'
	//otherwise returns false.
	bool operator<(const PeptideEntry& rhs) const;
	//generates a string in the gtf format and writes it to the specified ostream.
	std::ostream& to_gtf(const std::string& source, std::ostream& os = std::cout);
	//generates a bed line and writes it to the specified ostream.
	std::ostream& to_bed(std::ostream& os = std::cout);
	//generates a gct line and writes it to the specified ostream.
	std::ostream& to_gct(const std::string& geneID, std::vector<std::string> const& tissuevector, std::ostream& os = std::cout);
	//generates a bed line with ptms and writes it to the specified ostream.
	std::ostream& to_ptmbed(std::ostream& os = std::cout);

	//adds new tissuetags if they havent existed before.
	void add_tags(std::string const& tag, unsigned int sigPSMs, double quant);
	//adds ptms if a sequence matches another with ptms.
	void add_ptm(std::string const& ptmSequence);

private:
	//holds the peptide sequence.
	std::string		m_sequence;
	//number of found transcripts.
	unsigned int	m_numtranscripts;
	//is the peptide unique to one gene?
	bool			m_geneunique;
	//is the peptide unique to one transcript?
	bool			m_transcriptunique;
	//set of coordinates.
	std::set<PeptideCoordinates*, peptidecoords_p_compare>								m_pepcoordinates;
	//set of tissues.
	std::map<std::string, std::pair<std::vector<unsigned int>, std::vector<double>>>	m_tissuetags;
	//lowest startcoord of peptides that share the same sequence.
	unsigned int	m_startcoord;
	//highest endcoord of peptides that share the same sequence.
	unsigned int	m_endcoord;
	//pointer to the associated gene
	GeneEntry*		m_associatedgene;
	//this map holds the different PTMs for a peptide.
	std::map<std::string, std::map<std::string, PTMEntry>> m_pepforms;

	//this function is used to generate a date for a gct line.
	std::string tissue_quant_to_string(std::vector<std::string> const& tissueVector);
	//generates a mep that holds all PTMs in the given sequence.
	static std::map<std::string, PTMEntry> ptm_set(std::string sequence);
	//this function creates a coordinate_map_type and works similar to a CoordinateMapTypeConstructor.
	static CoordinateMapType create_coordinate_map_type(std::vector<GenomeCoordinates> const& genomeCoords);
};

//this comparator lets you compare PeptideEntry pointers.
//it just calls the operator < for the dereferenced pointers.
struct peptideentry_p_compare {
	bool operator()(const PeptideEntry* lhs, const PeptideEntry* rhs) const;
};

#endif
