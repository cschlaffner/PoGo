#ifndef MAPENTRY_H
#define MAPENTRY_H

#include "KmerMap.h"

//a map entry is used to connect a peptide string to all its variations 
//it maps a peptide to a gene, and the transcripts.
class MapEntry {
public:
	//ctr /dtr 
	explicit MapEntry(GeneEntry* geneentry_p);
	~MapEntry(void);

	//delegates to peptide_entry after checking if that peptide already 
	//exists and creating it if it doesn't
	unsigned int add_peptide(CoordinateWrapper& coordwrapper, const std::string& sequence, const std::string& tag, unsigned int sigPSMs, unsigned int genes, std::ofstream& ofstream, double quant, gene_id_map_t::value_type const& transcriptsEntry);
	//adds a transcript id to the mapping.
	void add_transcript_id(const std::string& transcriptID);
	//removes all peptides that are associated with a specific sequence.
	void remove_peptides();

	//compares two MapEntry objects. returns true if the lhs' GeneEntry is lesser than rhs' 
	bool operator<(const MapEntry& rhs) const;
	//calls the PeptideEntry::to_gtf metod for every peptide.
	std::ostream& to_gtf(const std::string& source, std::ostream& os = std::cout, bool chrincluded = true);
	//calls the PeptideEntr::to_bed metod for every peptide.
	std::ostream& to_bed(std::ostream& os = std::cout, bool chrincluded = true);
	//calls the PeptideEntry::to_gct metod for every peptide.
	std::ostream& to_gct(std::vector<std::string> const& tissuevector, std::ostream& os = std::cout, bool chrincluded = true);
	//calls the PeptideEntry::to_ptmbed metod for every peptide.
	std::ostream& to_ptmbed(std::ostream& os = std::cout, std::ostream&os2 = std::cout, bool chrincluded = true);
private:
	//pointer to the associated GeneEntry
	GeneEntry* m_p_gene_entry;
	//peptideentries, maps sequence without ptms to the corresponding PeptideEntry
	std::map<std::string, PeptideEntry*> m_peptide_entries;
	//holds all transcripts for this MapEntry
	std::set<std::string> m_transcripts;
};

//comparator for MapEntry pointers. dereferences the pointers and calls the operator< method.
struct mapentry_p_compare {
	bool operator()(const MapEntry* lhs, const MapEntry* rhs) const;
};

#endif
