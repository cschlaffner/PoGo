#if !defined	PEPTIDEMAPPERUTILS_H
#define			PEPTIDEMAPPERUTILS_H

//#include <map>
#include <vector>
#include <string>

//this is made for a maximum of 2 mismatches but can easily be extended 
//(for example just replace the ints with a vaector of int.)
//care: this change has to be applied to KmereMap::insert_into_gene_id_map

//this struct holds information about the position of a peptide in a protein.
//it also knows how many and where mismatches occured in the matching process.
struct position_mismatch_t {
private:
	//position of the peptide in the proteinsequence (zerobased)
	int m_position_in_protein;
	//first mismatch position in the proteinsequence (zerobased)
	int m_first_mismatch_positon;
	//second mismatch position in the proteinsequence (zerobased)
	int m_second_mismatch_positon;
public:
	//ctr
	position_mismatch_t(int posInProtein, int firstMismatch, int secondMismatch);
	//returns the positon of the peptide in the protein or -1 if the peptide does not match.
	int position_in_protein() const;
	//returns the position of the first mismatch or -1 if there was no mismatch.
	int first() const;
	//returns the position of the second mismatch or -1 if there was no mismatch.
	int second() const;
};

struct transcripts_t {
	//				 transcript id, mismatch positions
	typedef std::map<std::string, std::vector<position_mismatch_t>> transcript_id_map_t;
	transcript_id_map_t m_entries;
};

//				    gene id,  transcript ids
typedef std::map<std::string, transcripts_t> gene_id_map_t;


#endif
