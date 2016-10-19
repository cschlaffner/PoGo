#include "KmerMap.h"

KmerMap::match_protein_fun_t			KmerMap::m_match_protein(KmerMap::match_protein);
KmerMap::match_protein_backwards_fun_t KmerMap::m_match_protein_backwards(KmerMap::match_backwards);

KmerMap::KmerMap() : m_key_gen(this) {
	if ((GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_MISMATCHES > 1) &&
		GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ONE_IN_FIVE_MODE) {
		m_match_protein = KmerMap::match_protein_one_in_five_mode;
		m_match_protein_backwards = KmerMap::match_backwards_one_in_five_mode;
	}
}
KmerMap::~KmerMap() {}

void KmerMap::add_protein(ProteinEntry const& protein) {
	//this function digests a protein sequence and builds a KmerMap with the bits.
	std::string const& protein_sequence = protein.get_sequence();
	char const* cstr_protein_sequence = protein_sequence.c_str();

	if (protein_sequence.length() >= GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::KMER_LENGTH) {
		//<= n!
		unsigned int n(protein_sequence.length() - GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::KMER_LENGTH);
		std::string key;
		kmer_map_t::iterator it;
		KmerEntry curr_kmer;

		for (unsigned int i(0); i <= n; i++) {
			key = protein_sequence.substr(i, GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::KMER_LENGTH);
			it = m_kmers.find(key);
			if (it == m_kmers.end()) {
				//if the kmer doesnt exist yet, create it
				it = (m_kmers.insert(kmer_map_t::value_type(key, std::vector<KmerEntry>()))).first;
			}
			curr_kmer = KmerEntry((cstr_protein_sequence + i), &protein, i);
			if (it != m_kmers.end()) {
				//add the current KmerEntry
				it->second.push_back(curr_kmer);
			}
		}
	}
}


gene_id_map_t const& KmerMap::find_peptide(std::string const& peptide_string) {
	//this function generates a gene_id_map.
	//this map will map a gene_id to all related transcript ids and all the peptides and their 
	//position in the proteinsequence

	if (!m_gene_id_map.empty()) {
		m_gene_id_map.clear();
	}

	int set_key_returned = m_key_gen.set_original_key(peptide_string);
	int backwards_multiplier(0);

	std::string curr_key;
	std::vector<unsigned int> mismatches;

	size_t peptide_length(peptide_string.length());
	char const* cstr_peptide = peptide_string.c_str();

	if (set_key_returned >= 0) {
		while (m_key_gen.get_next_key(curr_key)) {
			if (m_kmers.count(curr_key) > 0) {
				kmer_map_t::const_iterator curr_entry_it = (m_kmers.find(curr_key));
				if (curr_entry_it != m_kmers.end()) {
					kmer_map_t::value_type const& curr_entry = *curr_entry_it;
					for (std::vector<KmerEntry>::const_iterator kmer_it = curr_entry.second.begin(); kmer_it != curr_entry.second.end(); ++kmer_it) {
						if (set_key_returned == 0) {
							if (m_match_protein(cstr_peptide, *kmer_it, mismatches, peptide_length)) {
								insert_into_gene_id_map(*kmer_it, mismatches);
							}
						} else if (set_key_returned == 1) {
							//this mode is used when only allowed_mismatches + 1 keys are generated. 
							//see PossibleKeyGenerator::set_original_key
							int offset(backwards_multiplier * GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::KMER_LENGTH);
							if (m_match_protein_backwards(cstr_peptide, *kmer_it, mismatches, peptide_length, offset)) {
								insert_into_gene_id_map(*kmer_it, mismatches, offset);
							}
						}
						mismatches.clear();
					}
				}
			}
			backwards_multiplier++;
		}
	}
	return m_gene_id_map;
}


bool KmerMap::match_protein(char const* peptideString, KmerEntry const& kmerEntry, std::vector<unsigned int>& mismatches, size_t peptideLength) {
	size_t protein_length(kmerEntry.m_p_protein->get_sequence().length() - kmerEntry.m_pos_in_protein);
	if (peptideLength <= protein_length) {
		char const* p_protein = (kmerEntry.m_p_protein->get_sequence().c_str()) + kmerEntry.m_pos_in_protein;
		for (size_t i = 0; i < peptideLength; i++) {
			if (peptideString[i] != p_protein[i]) {
				mismatches.push_back(i);
				if (mismatches.size() > GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_MISMATCHES) {
					return false;
				}
			}
		}
		return true;
	}
	return false;
}

bool KmerMap::match_protein_one_in_five_mode(char const* peptideString, KmerEntry const& kmerEntry, std::vector<unsigned int>& mismatches, size_t peptideLength) {
	size_t protein_length(kmerEntry.m_p_protein->get_sequence().length() - kmerEntry.m_pos_in_protein);
	if (peptideLength <= protein_length) {
		char const* p_protein = (kmerEntry.m_p_protein->get_sequence().c_str()) + kmerEntry.m_pos_in_protein;
		for (size_t i = 0; i < peptideLength; i++) {
			if (peptideString[i] != p_protein[i]) {
				if (mismatches.size() != 0) {
					if ((mismatches[mismatches.size() - 1] - i) <= GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_MISMATCHES) {
						return false;
					}
				}
				mismatches.push_back(i);
				if (mismatches.size() > GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_MISMATCHES) {
					return false;
				}
			}
		}
		return true;
	}
	return false;
}

bool KmerMap::match_backwards(char const* peptideString, KmerEntry const& kmerEntry, std::vector<unsigned int>& mismatches, size_t peptideLength, int offset) {
	if (kmerEntry.m_pos_in_protein >= offset) {
		size_t protein_length(kmerEntry.m_p_protein->get_sequence().length() - kmerEntry.m_pos_in_protein + offset);
		if (peptideLength <= protein_length) {
			char const* p_protein = (kmerEntry.m_p_protein->get_sequence().c_str()) + (kmerEntry.m_pos_in_protein - offset);
			for (size_t i = 0; i < peptideLength; i++) {
				if (peptideString[i] != p_protein[i]) {
					mismatches.push_back(i);
					if (mismatches.size() > GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_MISMATCHES) {
						return false;
					}
				}
			}
			return true;
		}
	}
	return false;
}

bool KmerMap::match_backwards_one_in_five_mode(char const* peptideString, KmerEntry const& kmerEntry, std::vector<unsigned int>& mismatches, size_t peptideLength, int offset) {
	if (kmerEntry.m_pos_in_protein >= offset) {
		size_t protein_length(kmerEntry.m_p_protein->get_sequence().length() - kmerEntry.m_pos_in_protein + offset);
		if (peptideLength <= protein_length) {
			char const* p_protein = (kmerEntry.m_p_protein->get_sequence().c_str()) + (kmerEntry.m_pos_in_protein - offset);
			for (size_t i = 0; i < peptideLength; i++) {
				if (peptideString[i] != p_protein[i]) {
					if (mismatches.size() != 0) {
						if ((mismatches[mismatches.size() - 1] - i) <= GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_MISMATCHES) {
							return false;
						}
					}
					mismatches.push_back(i);
					if (mismatches.size() > GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_MISMATCHES) {
						return false;
					}
				}
			}
			return true;
		}
	}
	return false;
}


void KmerMap::insert_into_gene_id_map(KmerEntry const& entry, std::vector<unsigned int> const& mismatches, int offset) {
	//inserts a found position into the gene id map	
	int pos_in_protein(entry.m_pos_in_protein - offset);
	std::string const& gene_id((entry.m_p_protein)->get_gene_id());
	std::string const& transcript_id((entry.m_p_protein)->get_transcript_id());
	if (m_gene_id_map.count(gene_id) == 0) {
		m_gene_id_map.insert(gene_id_map_t::value_type(gene_id, transcripts_t()));
	}
	if ((m_gene_id_map[gene_id].m_entries).count(transcript_id) == 0) {
		m_gene_id_map[gene_id].m_entries.insert(transcripts_t::transcript_id_map_t::value_type(transcript_id, std::vector<position_mismatch_t>()));
	}

	//care: if the position_mismatch_t struct is changed to accomodate more than 2 mismatches this has to be updated as well.
	m_gene_id_map[gene_id].m_entries[transcript_id].push_back(position_mismatch_t(pos_in_protein,
		(mismatches.size() > 0) ? pos_in_protein + mismatches[0] : -1,
		(mismatches.size() > 1) ? pos_in_protein + mismatches[1] : -1));
}

size_t KmerMap::size() const {
	return m_kmers.size();
}

bool KmerMap::contains(const std::string& key) const {
	return m_kmers.count(key) > 0;
}
