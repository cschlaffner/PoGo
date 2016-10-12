#include "PossibleKeyGenerator.h"
#include "KmereMap.h"

PossibleKeyGenerator::PossibleKeyGenerator(KmereMap const* k)
: m_keys_generated(false), m_kmeres(k) {}

int PossibleKeyGenerator::set_original_key(std::string const& key) {
	if (key.length() >= ((GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_MISMATCHES + 1) * GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::KMERE_LENGTH)) {
		//if key.length >= (allowed_mismatches+1)*kmere_lenght
		//in this case there is no way there are more than allowed_mismatches in every kmere.
		//one of those has to be matching perfectly:
		//protein:  TESTTESTTEST allowed mismatches = 2 
		//key:	    TEXTTEXTTEST  length = 15, 2 mismatches on the protein 
		//mismatches +1 * kmerelength = 15 therefore this if will be used
		//splitting the key into TEXT, TEXT, and TEST. the TEST does match the protein perfectly.
		//therefore the peptide will be found. this uses the match_backwards function.
		//if the key is longer this approach will also work. in the 2 mismatches case it generates 3 keys,
		//as opposed to the brute force key generation wich will produce 3331 unique keys.
		for (size_t i = 0; i <= GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_MISMATCHES; i++) {
			m_keys.push_back(key.substr((i * GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::KMERE_LENGTH), GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::KMERE_LENGTH));
		}
		m_keys_generated = true;
		m_curr_element = m_keys.begin();
		return 1;
	}
	if (key.length() >= GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::KMERE_LENGTH) {
		set_short_original_key(key.substr(0, GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::KMERE_LENGTH));
		return 0;
	}
	return -1;
}

void PossibleKeyGenerator::set_short_original_key(std::string const& key) {
	if (GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_MISMATCHES > 0) {
		m_key = key;
		m_keys_generated = false;
	} else {
		m_keys.push_back(key);
		m_keys_generated = true;
		m_curr_element = m_keys.begin();
	}
}

bool PossibleKeyGenerator::get_next_key(std::string& key) {
	if (!m_keys_generated) {
		//if the keys are not yet generated
		if (m_key.length() > 0 && GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_AMINO_ACIDS.size() > 0) {
			//generate them and return the first one (lazy generation)
			generate_keys();
			m_keys_generated = true;
			m_curr_element = m_keys.begin();
			return get_next_key(key);
		}
	} else {
		//otherwise just return the next key
		if (m_curr_element != m_keys.end()) {
			key = (*m_curr_element);
			++m_curr_element;
			return true;
		}
		//and if the end is reached, clear out the old keys and wait for the next peptide.
		m_keys_generated = false;
		m_key.clear();
		m_keys.clear();
	}
	return false;
}

void PossibleKeyGenerator::generate_keys_one_mismatch() {
	//generates all keys with one mismatch.
	for (size_t kmere_it = 0; kmere_it < GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::KMERE_LENGTH; kmere_it++) {
		char original = m_key[kmere_it];
		for (size_t aa_it = 0; aa_it < GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_AMINO_ACIDS.size(); aa_it++) {
			m_key[kmere_it] = GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_AMINO_ACIDS[aa_it];
			if (m_kmeres->contains(m_key)) {
				m_keys.push_back(m_key);
			}
		}
		m_key[kmere_it] = original;
	}
}

void PossibleKeyGenerator::generate_keys_two_mismatches() {
	//generates all keys with 2 mismatches. this function is costly and should be avoided if possible.
	static std::vector<std::pair<char, char>> combinations;
	static bool combiantions_generated(false);
	if (!combiantions_generated) {
		for (size_t i(0); i < GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_AMINO_ACIDS.size(); ++i) {
			for (size_t j = 0; j < GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_AMINO_ACIDS.size(); j++) {
				combinations.push_back(std::make_pair(GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_AMINO_ACIDS[i], GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_AMINO_ACIDS[j]));
			}
		}
		combiantions_generated = true;
	}
	
	for (size_t i = 0; i < GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::KMERE_LENGTH - 1; i++) {
		char original = m_key[i];
		for (size_t j = i; j < GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::KMERE_LENGTH; j++) {
			char original_inner = m_key[j];
			for (size_t k = 0; k < combinations.size(); k++) {
				if (m_key[i] != combinations[k].first) {
					m_key[i] = combinations[k].first;
				}
				m_key[j] = combinations[k].second;
				if (m_kmeres->contains(m_key)) {
					m_keys.push_back(m_key);
				}
			}
			m_key[j] = original_inner;
		}
		m_key[i] = original;
	}
}

void PossibleKeyGenerator::generate_keys() {
	if (GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_MISMATCHES == 1) {
		generate_keys_one_mismatch();
	} else if (GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_MISMATCHES == 2) {
		if (GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ONE_IN_FIVE_MODE) {
			generate_keys_one_mismatch();
		} else {
			generate_keys_two_mismatches();
		}
	}
	//add to here if three mismatches are required at some point.
}
