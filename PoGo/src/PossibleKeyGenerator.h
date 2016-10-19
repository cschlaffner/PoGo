#if !defined	POSSIBLE_KEY_GENERATOR
#define			POSSIBLE_KEY_GENERATOR
#include "PeptideMapperUtils.h"
#include "Globals.h"
#include <set>

class KmerMap;

class PossibleKeyGenerator {
public:
	//ctr
	explicit PossibleKeyGenerator(KmerMap const* k);

	//returns
	//-1: key is shorter than KMER_LENGTH 
	// 0: key is shorter than (ALLOWED MISMATCHES + 1) * KMER_LENGTH 
	// 1: key is >= KMER_LENGTH
	int set_original_key(std::string const& key);
	
	//assigns the next key to the passed string. 
	//returns true if there is another key otherwise returns false.
	//the first call might be a bit slower as the keys are generated in a lazy way.
	bool get_next_key(std::string& key);

private:
	//string used to generate keys.
	std::string m_key;
	//are there any generated keys left? 
	bool m_keys_generated;

	// kmermap, to check if a certain key exists in there, 
	//otherwise the key will not be added to the generated keys.
	KmerMap const* m_kmers;
	
	//current set of generated keys.
	std::vector<std::string> m_keys;
	//iterator pointing to the current element in m_keys.
	std::vector<std::string>::iterator m_curr_element;

	//used to set the kmerlength long key.
	void set_short_original_key(std::string const& key);
	
	//calls the other generator functions.
	void generate_keys();
	//generates keys for one mismatch matching.
	void generate_keys_one_mismatch();
	//generates keys for two mismatch matching.
	void generate_keys_two_mismatches();

};

#endif
