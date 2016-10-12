#ifndef GENEENTRY_H
#define GENEENTRY_H

//#include "Utils\PeptideCoordinates.h"

class GeneEntry {
public:
	//ctr/dtr
	GeneEntry(void);
	explicit GeneEntry(std::string const& gtfGeneLine);
	~GeneEntry(void);
	// end ctr/dtr

	//returns the gene ID
	std::string get_id() const;
	//returns the type
	std::string get_type() const;
	//returns the genes status
	std::string get_status() const;
	//returns the gene name
	std::string get_name() const;

	//comapres two genes. 
	//returns true if the chromosome number is smaller than rhs' chromosome number 
	//otherwise returns true if the startposition in the chromosome is smaller 
	//otherwise returns true if the endposition in the chromosome is smaller 
	//otherwise returns false.	
	bool operator<(const GeneEntry& rhs) const;
	//converts a gene into a gtf line and prints it to the given outputstream.
	std::ostream& to_gtf(const std::string& source, std::ostream& os = std::cout);

	//looks for the text specified in GENOME_MAPPER_GLOBALS::ID::GENE_ID and returns the ID (including the trailing number of length GENOME_MAPPER_GLOBALS::ID::LENGTH - GENOME_MAPPER_GLOBALS::ID::GENE_ID.length()).
	std::string static extract_gene_id(std::string gtfgeneline);
	//looks for the text specified in GENOME_MAPPER_GLOBALS::ID::TRANSCRIPT_ID and returns the ID (including the trailing number of length GENOME_MAPPER_GLOBALS::ID::LENGTH - GENOME_MAPPER_GLOBALS::ID::TRANSCRIPT_ID.length()).
	std::string static extract_transcript_id(std::string gtfgeneline);

private:
	//the genomic coordinates for a gene are stored here.
	GenomeCoordinates m_coord;
	//gene id as extracted by extract_gene_id();
	std::string m_id;
	//type (protein_coding,...)
	std::string m_type;
	//status (KNOWN,...)
	std::string m_status;
	//gene name (gene symbol)
	std::string m_gene_name;
	//tags. (ncRNA_host,...)
	std::vector<std::string> m_tags;

	//private member function called in constructor to initialize a gene entry
	void init(const std::string& gtfGeneLine);
	//private member function called in constructor to initialize a gene entry
	void init(std::string ID, GenomeCoordinates coordinates, std::string type, std::string status, std::string geneName, std::vector<std::string> tags = std::vector<std::string>());
	
	//after tokenizing the geneline these functions can be used to extract information.
	
	//extracts the gene type (proein_coding,...) 
	static std::string extract_type(const std::vector<std::string> & tokens);
	//extracts the gene status (KNOWN,...)
	static std::string extract_status(const std::vector<std::string> & tokens);
	//extracts the gene symbol
	static std::string extract_gene_name(const std::vector<std::string> & tokens);
	//extracts all possible tags to  be used with extract_by_tag.
	static std::vector<std::string> extract_tags(const std::vector<std::string> & tokens);
	//tag takes the name of the tag to extract, tagList contains all tags separated by \. 
	//possible values for tag could be gene_name, gene_type,...
	static std::vector<std::string> extract_by_tag(const std::string& tag, const std::string& tagList);
};

#endif
