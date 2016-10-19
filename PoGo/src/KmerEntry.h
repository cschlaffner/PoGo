#if !defined	KMER_ENTRY
#define			KMER_ENTRY

#include "PeptideEntry.h"

//a kmer Entry holds information on a kmer.
//it has pointers to the proteinsequence and associated protein to which it belongs, 
//and knows its position therein
struct KmerEntry {
	//the pointer to the first letter of the key in the protein sequence
	const char* m_p_key;
	//the pointer to the protein. 
	ProteinEntry const* m_p_protein;
	//the (0 based) index of the first letter of the kmer in the protein string. 
	int m_pos_in_protein;

	//ctr
	KmerEntry(const char* key, ProteinEntry const* protein, unsigned int pos_in_protein)
		: m_p_key(key), m_p_protein(protein), m_pos_in_protein(pos_in_protein) {}

	KmerEntry() : m_p_key(nullptr), m_p_protein(nullptr), m_pos_in_protein(0) {}
};
#endif
