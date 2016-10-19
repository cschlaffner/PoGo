#if !defined	KMER_MAP
#define			KMER_MAP

#include "KmerEntry.h"
#include "PossibleKeyGenerator.h"
#include <unordered_map>

class KmerMap {
public:
	KmerMap();
	~KmerMap();
	//digests and adds a protein to the map.
	void add_protein(ProteinEntry const& protein);
	//searches for all matches (imperfect matching, set via PEPTIDE_MAPPER) and returns a map that contains all finds.
	gene_id_map_t const& find_peptide(std::string const& peptideString);
	// returns true if a kmer (key) is in the digested proteins
	bool contains(const std::string& key) const;
	//returns the number of fragments that were created during digestion
	size_t size() const;

private:
	//kmerMap
	//in the gencode fasta file, data analysis showed that the average size of the vector is 18 elements (for approx. 93 000 proteins)
	//and the map size is at approx 1.8m elements.
	typedef std::unordered_map<std::string, std::vector<KmerEntry>> kmer_map_t;
	kmer_map_t m_kmers;

	//gene_id map
	gene_id_map_t m_gene_id_map;
	//inserts a found peptide into the current gene id map.
	void insert_into_gene_id_map(KmerEntry const& entry, std::vector<unsigned int> const& mismatches, int offset = 0);

	//the keygenerator will generate keys for one and two mismatch matching.
	PossibleKeyGenerator m_key_gen;

	//these function pointers allow for changing the "one in five mode". 
	//it also makes it easy to implement other matching algorithms (just change to the one 
	// you want to use in the kmermap constructor.)
	typedef bool (*match_protein_fun_t)				(char const*, KmerEntry const&, std::vector<unsigned int>&, size_t);
	typedef bool (*match_protein_backwards_fun_t)	(char const*, KmerEntry const&, std::vector<unsigned int>&, size_t, int offset);

	static match_protein_fun_t				m_match_protein;
	static match_protein_backwards_fun_t	m_match_protein_backwards;

	//the basic forward matching algorithm
	static bool match_protein(char const* peptideString, KmerEntry const& kmerEntry, std::vector<unsigned int>& mismatches, size_t peptideLength);
	//forward matching with one in five stop criterion
	static bool match_protein_one_in_five_mode(char const* peptideString, KmerEntry const& kmerEntry, std::vector<unsigned int>& mismatches, size_t peptideLength);
	//backwards matching functionality. this is possible by using cstrings.
	static bool match_backwards(char const* peptideString, KmerEntry const& kmerEntry, std::vector<unsigned int>& mismatches, size_t peptideLength, int offset);
	//backwards matching functionality with one in five stop criterion.
	static bool match_backwards_one_in_five_mode(char const* peptideString, KmerEntry const& kmerEntry, std::vector<unsigned int>& mismatches, size_t peptideLength, int offset);
};


#endif
