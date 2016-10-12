//#include "KmereMap.h"

KmereMap::match_protein_fun_t			KmereMap::m_match_protein(KmereMap::match_protein);
KmereMap::match_protein_backwards_fun_t KmereMap::m_match_protein_backwards(KmereMap::match_backwards);

KmereMap::KmereMap() : m_key_gen(this) {
	if ((GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_MISMATCHES > 1) &&
		GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ONE_IN_FIVE_MODE) {
		m_match_protein = KmereMap::match_protein_one_in_five_mode;
		m_match_protein_backwards = KmereMap::match_backwards_one_in_five_mode;
	}
}
KmereMap::~KmereMap() {}

void KmereMap::add_protein(ProteinEntry const& protein) {
	//this function digests a protein sequence and builds a KmereMap with the bits.
	std::string const& protein_sequence = protein.get_sequence();
	char const* cstr_protein_sequence = protein_sequence.c_str();

	if (protein_sequence.length() >= GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::KMERE_LENGTH) {
		//<= n!
		unsigned int n(protein_sequence.length() - GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::KMERE_LENGTH);
		std::string key;
		kmere_map_t::iterator it;
		KmereEntry curr_kmere;

		for (unsigned int i(0); i <= n; i++) {
			key = protein_sequence.substr(i, GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::KMERE_LENGTH);
			it = m_kmeres.find(key);
			if (it == m_kmeres.end()) {
				//if the kmere doesnt exist yet, create it
				it = (m_kmeres.insert(kmere_map_t::value_type(key, std::vector<KmereEntry>()))).first;
			}
			curr_kmere = KmereEntry((cstr_protein_sequence + i), &protein, i);
			if (it != m_kmeres.end()) {
				//add the current KmereEntry
				it->second.push_back(curr_kmere);
			}
		}
	}
}


gene_id_map_t const& KmereMap::find_peptide(std::string const& peptide_string) {
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
			if (m_kmeres.count(curr_key) > 0) {
				kmere_map_t::const_iterator curr_entry_it = (m_kmeres.find(curr_key));
				if (curr_entry_it != m_kmeres.end()) {
					kmere_map_t::value_type const& curr_entry = *curr_entry_it;
					for (std::vector<KmereEntry>::const_iterator kmere_it = curr_entry.second.begin(); kmere_it != curr_entry.second.end(); ++kmere_it) {
						if (set_key_returned == 0) {
							if (m_match_protein(cstr_peptide, *kmere_it, mismatches, peptide_length)) {
								insert_into_gene_id_map(*kmere_it, mismatches);
							}
						} else if (set_key_returned == 1) {
							//this mode is used when only allowed_mismatches + 1 keys are generated. 
							//see PossibleKeyGenerator::set_original_key
							int offset(backwards_multiplier * GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::KMERE_LENGTH);
							if (m_match_protein_backwards(cstr_peptide, *kmere_it, mismatches, peptide_length, offset)) {
								insert_into_gene_id_map(*kmere_it, mismatches, offset);
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


bool KmereMap::match_protein(char const* peptideString, KmereEntry const& kmereEntry, std::vector<unsigned int>& mismatches, size_t peptideLength) {
	size_t protein_length(kmereEntry.m_p_protein->get_sequence().length() - kmereEntry.m_pos_in_protein);
	if (peptideLength <= protein_length) {
		char const* p_protein = (kmereEntry.m_p_protein->get_sequence().c_str()) + kmereEntry.m_pos_in_protein;
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

bool KmereMap::match_protein_one_in_five_mode(char const* peptideString, KmereEntry const& kmereEntry, std::vector<unsigned int>& mismatches, size_t peptideLength) {
	size_t protein_length(kmereEntry.m_p_protein->get_sequence().length() - kmereEntry.m_pos_in_protein);
	if (peptideLength <= protein_length) {
		char const* p_protein = (kmereEntry.m_p_protein->get_sequence().c_str()) + kmereEntry.m_pos_in_protein;
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

bool KmereMap::match_backwards(char const* peptideString, KmereEntry const& kmereEntry, std::vector<unsigned int>& mismatches, size_t peptideLength, int offset) {
	if (kmereEntry.m_pos_in_protein >= offset) {
		size_t protein_length(kmereEntry.m_p_protein->get_sequence().length() - kmereEntry.m_pos_in_protein + offset);
		if (peptideLength <= protein_length) {
			char const* p_protein = (kmereEntry.m_p_protein->get_sequence().c_str()) + (kmereEntry.m_pos_in_protein - offset);
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

bool KmereMap::match_backwards_one_in_five_mode(char const* peptideString, KmereEntry const& kmereEntry, std::vector<unsigned int>& mismatches, size_t peptideLength, int offset) {
	if (kmereEntry.m_pos_in_protein >= offset) {
		size_t protein_length(kmereEntry.m_p_protein->get_sequence().length() - kmereEntry.m_pos_in_protein + offset);
		if (peptideLength <= protein_length) {
			char const* p_protein = (kmereEntry.m_p_protein->get_sequence().c_str()) + (kmereEntry.m_pos_in_protein - offset);
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


void KmereMap::insert_into_gene_id_map(KmereEntry const& entry, std::vector<unsigned int> const& mismatches, int offset) {
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

size_t KmereMap::size() const {
	return m_kmeres.size();
}

bool KmereMap::contains(const std::string& key) const {
	return m_kmeres.count(key) > 0;
}
