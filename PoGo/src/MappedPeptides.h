#ifndef MAPPEDPEPTIDES_H
#define MAPPEDPEPTIDES_H


#include "MapEntry.h"

//
//the mapped peptides class connects the different classes
//and offers functions to wite all found peptides.
class MappedPeptides {
public:
	//ctr / dtr
	MappedPeptides(void);
	~MappedPeptides(void);

	//adds a new gene from a gtfline.
	assembly add_gene_from_gtf(std::string const& gtfGeneLine);
	//maps a transcript id to a gene id
	void add_transcript_id_to_gene(std::string const& gtfTranscriptLine);
	//passes the peptide insertion and lookup to the MapEntry.
	void add_peptide(CoordinateWrapper& coordwrapper, const std::string& sequence, const std::string& tag, unsigned int sigPSMs, unsigned int genes, std::ofstream& ofstream, double quant, gene_id_map_t::value_type const& transcriptsEntry);

	//print functions
	//converts all peptides to gtf lines
	void to_gtf(const std::string& filename, const std::string& source, assembly assem = primary, bool chrincluded = true);
	void to_gtf(assembly assem, std::string source, std::ostream& os = std::cout, bool chrincluded = true);
	//converts all peptides to bed lines
	void to_bed(std::string filename, assembly assem = primary, bool chrincluded = true);
	void to_bed(assembly assem, std::ostream& os = std::cout, bool chrincluded = true);
	//converts all peptides to gct lines
	void to_gct(std::string filename, assembly assem = primary, bool chrincluded = true);
	void to_gct(assembly assem, std::ostream& os = std::cout, bool chrincluded = true);
	//converts all peptides to _ptm.bed lines
	void to_ptmbed(std::string filename, std::string filename2, assembly assem = primary, bool chrincluded = true);
	void to_ptmbed(assembly assem, std::ostream& os = std::cout, std::ostream& os2 = std::cout, bool chrincluded = true);
	//removes all peptides from the MappedPeptides.
	void remove_all_peptides();
private:
	typedef std::map<std::string, MapEntry*> map_t; 
	//maps the peptide sequence (without PTM) to the corresponding MapEntry (chromosomes).
	map_t		 m_mapping;
	//maps the peptide sequence (without PTM) to the corresponding MapEntry (patches, haplotypes, scaffolds).
	map_t		m_mapping_phs;
	//counts the number of found peptides. similar to m_mapping.size().
	unsigned int m_count_peptides;
	unsigned int m_count_peptides_phs;

	//maps a tissue to an index.
	std::map<std::string, unsigned int> m_tissuemap;
	//the index, to be incremented durig the program.
	unsigned int m_tissueindex;

	//converts the above map to a string that contains all found tissues, delimited by '\t' and sorted.
	std::string tissuemap_to_sorted_string(const std::string& sep);

};

#endif
