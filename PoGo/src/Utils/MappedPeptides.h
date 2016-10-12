#ifndef MAPPEDPEPTIDES_H
#define MAPPEDPEPTIDES_H


//#include "MapEntry.h"

//
//the mapped peptides class connects the different classes
//and offers functions to wite all found peptides.
class MappedPeptides {
public:
	//ctr / dtr
	MappedPeptides(void);
	~MappedPeptides(void);

	//adds a new gene from a gtfline.
	void add_gene_from_gtf(std::string const& gtfGeneLine);
	//maps a transcript id to a gene id
	void add_transcript_id_to_gene(std::string const& gtfTranscriptLine);
	//passes the peptide insertion and lookup to the MapEntry.
	void add_peptide(CoordinateWrapper& coordwrapper, const std::string& sequence, const std::string& tag, unsigned int sigPSMs, unsigned int genes, std::ofstream& ofstream, double quant, gene_id_map_t::value_type const& transcriptsEntry);

	//print functions
	//converts all peptides to gtf lines
	void to_gtf(const std::string& filename, const std::string& source);
	void to_gtf(std::string source, std::ostream& os = std::cout);
	//converts all peptides to bed lines
	void to_bed(std::string filename);
	void to_bed(std::ostream& os = std::cout);
	//converts all peptides to gct lines
	void to_gct(std::string filename);
	void to_gct(std::ostream& os = std::cout);
	//converts all peptides to _ptm.bed lines
	void to_ptmbed(std::string filename);
	void to_ptmbed(std::ostream& os = std::cout);
	//removes all peptides from the MappedPeptides.
	void remove_all_peptides();
private:
	typedef std::map<std::string, MapEntry*> map_t; 
	//maps the peptide sequence (without PTM) to the corresponding MapEntry.
	map_t		 m_mapping;
	//counts the number of found peptides. similar to m_mapping.size().
	unsigned int m_count_peptides;

	//maps a tissue to an index.
	std::map<std::string, unsigned int> m_tissuemap;
	//the index, to be incremented durig the program.
	unsigned int m_tissueindex;

	//converts the above map to a string that contains all found tissues, delimited by '\t' and sorted.
	std::string tissuemap_to_sorted_string(const std::string& sep);

};

#endif
