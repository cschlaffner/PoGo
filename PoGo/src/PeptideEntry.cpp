#include "PeptideEntry.h"
#include "PeptideMapperUtils.h"

#include <numeric>


PeptideEntry::PeptideEntry(GeneEntry* associatedgene) :
m_sequence(),
m_numtranscripts(0),
m_geneunique(true),
m_transcriptunique(true),
m_pepcoordinates(std::set<PeptideCoordinates*, peptidecoords_p_compare>()),
m_tissuetags(std::map<std::string, std::pair<std::vector<unsigned int>, std::vector<double>>>()),
m_startcoord(0),
m_endcoord(0),
m_associatedgene(associatedgene),
m_pepforms(std::map<std::string, std::map<std::string, PTMEntry>>()) {}
PeptideEntry::~PeptideEntry(void) {}

std::map<std::string, PTMEntry> PeptideEntry::ptm_set(std::string sequence) {
	std::map<std::string, PTMEntry> map = std::map<std::string, PTMEntry>();

	size_t start_b = sequence.find_first_of("(");
	size_t end_b = sequence.find_first_of(")");
	std::string name = "";
	while (start_b != std::string::npos && end_b != std::string::npos) {
		name = sequence.substr(start_b + 1, (end_b - start_b - 1));
		sequence.erase(sequence.begin() + start_b, sequence.begin() + end_b + 1);
		if (EnumStringMapper::ptm_to_colour(name).compare("255,51,153") != 0) {
			if (start_b != 0) {
				start_b = start_b - 1;
			}
			if (name.compare("") != 0) {
				if (map.count(name) == 0) {
					map.insert(std::pair<std::string, PTMEntry>(name, PTMEntry(name, start_b, start_b)));
				} else {
					map.at(name).add_coord(start_b);
				}
			}
		}
		start_b = sequence.find_first_of("(");
		end_b = sequence.find_first_of(")");
	}
	return map;
}

bool PeptideEntry::operator<(const PeptideEntry& rhs) const {
	return (m_startcoord < rhs.m_startcoord) || (m_startcoord == rhs.m_startcoord && m_endcoord < rhs.m_endcoord)
		|| (m_startcoord == rhs.m_startcoord && m_endcoord == rhs.m_endcoord && m_sequence.compare(rhs.m_sequence) < 0);
}

std::ostream& PeptideEntry::to_gtf(const std::string& source, std::ostream& os) {
	std::string sequence_add = "";
	unsigned int count = 0;
	for (std::set<PeptideCoordinates*, peptidecoords_p_compare>::iterator it = m_pepcoordinates.begin(); it != m_pepcoordinates.end(); ++it) {
		count += 1;
		if (m_pepcoordinates.size() > 1) {
			std::stringstream ss_seq_add;
			ss_seq_add << "." << count;
			sequence_add = ss_seq_add.str();
		}
		if (count > 1) {
			os << "\n";
		}

		os << coordinates_to_gtf_string((*it)->get_transcript_coordinates(), "transcript", false, source);
		os << "gene_id \"" << m_associatedgene->get_id() << "\"; transcript_id \"" << m_associatedgene->get_id() << "." << m_sequence << sequence_add << "\"; gene_type \"" << m_associatedgene->get_type() << "\"; gene_status \"" << m_associatedgene->get_status() << "\"; gene_name \"" << m_associatedgene->get_name();
		os << "\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"" << m_associatedgene->get_name() << "." << m_sequence << sequence_add << "\";";
		os << " tag \"Transcripts:" << m_numtranscripts << "\";";

		for (std::map<std::string, std::pair<std::vector<unsigned int>, std::vector<double>>>::iterator it_tissue = m_tissuetags.begin(); it_tissue != m_tissuetags.end(); ++it_tissue) {
			os << " tag \"" << it_tissue->first << ":";
			std::stringstream ss_tissue;
			for (size_t i_tissue(0); i_tissue < it_tissue->second.first.size(); ++i_tissue) {
				if (i_tissue > 0) {
					os << "/";
					ss_tissue << "/";
				}
				os << it_tissue->second.first.at(i_tissue);
				ss_tissue << it_tissue->second.second.at(i_tissue);
			}
			os << " sig PSMs" << " " << ss_tissue.str() << " Quant" << "\";";
		}

		std::vector<GenomeCoordinates> exon_coordinates = (*it)->get_exon_coordinates();
		unsigned int exon_count = 0;
		std::string exon_add = "";
		for (size_t exit(0); exit < exon_coordinates.size(); ++exit) {
			exon_count += 1;
			if (exon_coordinates.size() > 1) {
				std::stringstream ss_exadd;
				ss_exadd << "." << exon_count;
				exon_add = ss_exadd.str();
			}

			os << "\n";
			os << coordinates_to_gtf_string(exon_coordinates.at(exit), "exon", false, source);
			os << "gene_id \"" << m_associatedgene->get_id() << "\"; transcript_id \"" << m_associatedgene->get_id() << "." << m_sequence << sequence_add << "\"; gene_type \"" << m_associatedgene->get_type() << "\"; gene_status \"" << m_associatedgene->get_status() << "\"; gene_name \"" << m_associatedgene->get_name();
			os << "\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"" << m_associatedgene->get_name() << "." << m_sequence << sequence_add << "\";";
			os << " exon_number " << exon_count << "; exon_id \"" << m_associatedgene->get_id() << "." << m_sequence << sequence_add << exon_add << "\";";
		}
	}
	return os;
}

std::ostream& PeptideEntry::to_bed(std::ostream& os) {
	for (std::set<PeptideCoordinates*, peptidecoords_p_compare>::iterator it = m_pepcoordinates.begin(); it != m_pepcoordinates.end(); ++it) {
		os << coordinates_to_bed_string((*it)->get_transcript_coordinates(), m_sequence);

		std::stringstream exon_starts;
		std::stringstream exon_lengths;
		std::vector<GenomeCoordinates> exon_coordinates = (*it)->get_exon_coordinates();
		unsigned int exon_count = 0;
		for (size_t exit(0); exit < exon_coordinates.size(); ++exit) {
			exon_count += 1;
			if (exon_count > 1) {
				exon_starts << ",";
				exon_lengths << ",";
			}
			unsigned int exon_start = exon_coordinates.at(exit).start - (*it)->get_transcript_coordinates().start;
			unsigned int exon_length = exon_coordinates.at(exit).end - exon_coordinates.at(exit).start + 1;
			exon_starts << exon_start;
			exon_lengths << exon_length;
		}
		std::string colour = "128,128,128";
		if (m_geneunique == true && m_transcriptunique == false) {
			colour = "0,0,0";
		} else if (m_geneunique == true && m_transcriptunique == true) {
			colour = "204,0,0";
		}

		os << colour << "\t" << exon_count << "\t" << exon_lengths.str() << "\t" << exon_starts.str() << "\n";
	}
	return os;
}

std::ostream& PeptideEntry::to_gct(const std::string& geneID, std::vector<std::string> const& tissuevector, std::ostream& os) {
	std::string sequence_add = "";
	unsigned int count = 0;
	for (std::set<PeptideCoordinates*, peptidecoords_p_compare>::iterator it = m_pepcoordinates.begin(); it != m_pepcoordinates.end(); ++it) {
		count += 1;
		if (m_pepcoordinates.size() > 1) {
			std::stringstream ss_seqadd;
			ss_seqadd << "." << count;
			sequence_add = ss_seqadd.str();
		}

		os << geneID << "." << m_sequence << sequence_add << "\t\"" << geneID << "|@";

		std::vector<GenomeCoordinates> exoncoords = (*it)->get_exon_coordinates();
		os << coordinates_to_gct_string(exoncoords) << "|\"";

		os << tissue_quant_to_string(tissuevector) << "\n";
	}
	return os;
}

std::string PeptideEntry::tissue_quant_to_string(std::vector<std::string> const& tissuevector) {
	std::stringstream ss;
	for (size_t i(1); i < tissuevector.size(); ++i) {
		ss << "\t";
		if (m_tissuetags.count(tissuevector.at(i)) == 1) {
			double sum = std::accumulate(m_tissuetags.at(tissuevector.at(i)).second.begin(), m_tissuetags.at(tissuevector.at(i)).second.end(), 0.0);
			double mean = sum / m_tissuetags.at(tissuevector.at(i)).second.size();
			ss << mean;
		}
	}
	return ss.str();
}

std::ostream& PeptideEntry::to_ptmbed(std::ostream& os) {
	for (std::map<std::string, std::map<std::string, PTMEntry>>::iterator ptm_it = m_pepforms.begin(); ptm_it != m_pepforms.end(); ++ptm_it) {
		for (std::map<std::string, PTMEntry>::iterator ptm_single_it = ptm_it->second.begin(); ptm_single_it != ptm_it->second.end(); ++ptm_single_it) {
			std::vector<std::pair<PeptideCoordinates*, GenomeCoordinates>> coord = ptm_single_it->second.get_genome_coordinates();
			for (std::vector<std::pair<PeptideCoordinates*, GenomeCoordinates>>::iterator coord_it = coord.begin(); coord_it != coord.end(); ++coord_it) {
				std::string bed_string = coordinates_to_short_bed_string(coord_it->first->get_transcript_coordinates(), ptm_it->first);
				std::vector<GenomeCoordinates> exon_coords = coord_it->first->get_exon_coordinates();

				std::stringstream exon_starts;
				std::stringstream exon_lengths;
				unsigned int exon_count = 0;
				for (size_t exit(0); exit < exon_coords.size(); ++exit) {
					exon_count += 1;
					if (exon_count > 1) {
						exon_starts << ",";
						exon_lengths << ",";
					}
					unsigned int exon_start = exon_coords.at(exit).start - coord_it->first->get_transcript_coordinates().start;
					unsigned int exon_length = exon_coords.at(exit).end - exon_coords.at(exit).start + 1;

					exon_starts << exon_start;
					exon_lengths << exon_length;
				}
				std::string colour = EnumStringMapper::ptm_to_colour(ptm_single_it->first);
				os << bed_string << coord_it->second.start - 1 << "\t" << coord_it->second.end << "\t" << "\t" << colour << "\t" << exon_count << "\t" << exon_lengths.str() << "\t" << exon_starts.str() << "\n";
			}
		}
	}
	return os;
}

void PeptideEntry::add_tags(std::string const& tag, unsigned int sigPSMs, double quant) {
	if (m_tissuetags.count(tag) == 0) {
		m_tissuetags.insert(std::pair<std::string, std::pair<std::vector<unsigned int>, std::vector<double>>>(tag, std::pair<std::vector<unsigned int>, std::vector<double>>(std::vector<unsigned int>(), std::vector<double>())));
	}
	m_tissuetags.at(tag).first.push_back(sigPSMs);
	m_tissuetags.at(tag).second.push_back(quant);
}

void PeptideEntry::add_ptm(std::string const& ptmsequence) {
	//adds the modfication to the peptide entry
	m_pepforms.insert(std::pair<std::string, std::map<std::string, PTMEntry>>(ptmsequence, ptm_set(ptmsequence)));
	for (std::map<std::string, PTMEntry>::iterator ptm_it = m_pepforms.at(ptmsequence).begin(); ptm_it != m_pepforms.at(ptmsequence).end(); ++ptm_it) {
		for (std::set<PeptideCoordinates*, peptidecoords_p_compare>::iterator peptide_coordinates_it = m_pepcoordinates.begin(); peptide_coordinates_it != m_pepcoordinates.end(); ++peptide_coordinates_it) {
			std::vector<GenomeCoordinates> genomic_coordinates = (*peptide_coordinates_it)->find_coordinates(ptm_it->second.get_range());
			GenomeCoordinates ptm_coordinates;
			for (size_t j(0); j < genomic_coordinates.size(); ++j) {
				if (j == 0) {
					ptm_coordinates = genomic_coordinates.at(j);
				} else {
					if (ptm_coordinates.start > genomic_coordinates.at(j).start) {
						ptm_coordinates.start = genomic_coordinates.at(j).start;
					}
					if (ptm_coordinates.end < genomic_coordinates.at(j).end) {
						ptm_coordinates.end = genomic_coordinates.at(j).end;
					}
				}
			}
			ptm_it->second.add_genome_coordinates(*peptide_coordinates_it, ptm_coordinates);
		}
	}
}

CoordinateMapType PeptideEntry::create_coordinate_map_type(std::vector<GenomeCoordinates> const& genomecoords) {
	CoordinateMapType coord_map = CoordinateMapType();
	Coordinates prot_coord = Coordinates();
	Coordinates prev_prot_coord = Coordinates();
	prev_prot_coord.Cterm = off3;
	prev_prot_coord.Nterm = off3;
	prev_prot_coord.start = 0;
	prev_prot_coord.end = 0;

	for (unsigned int i(0); i < genomecoords.size(); ++i) {
		GenomeCoordinates genome_coordinates = genomecoords.at(i);
		if (genome_coordinates.strand == rev) {
			unsigned int start = genome_coordinates.end;
			unsigned int end = genome_coordinates.start;
			genome_coordinates.start = start;
			genome_coordinates.end = end;
		}

		if (prev_prot_coord.Cterm != off3) {
			prot_coord.Nterm = Offset(3 - prev_prot_coord.Cterm);
		} else {
			prot_coord.Nterm = off3;
		}

		int length;

		if (genome_coordinates.strand == fwd) {
			length = genome_coordinates.end - genome_coordinates.start + 1;
		} else {
			length = genome_coordinates.start - genome_coordinates.end + 1;
		}
		// calc cterm
		if (length % 3 == 0) {
			if (prot_coord.Nterm != off3) {
				prot_coord.Cterm = Offset(3 - prot_coord.Nterm);
			} else {
				prot_coord.Cterm = off3;
			}
		} else if (length % 3 == 2) {
			if (prot_coord.Nterm == off3) {
				prot_coord.Cterm = off2;
			} else if (prot_coord.Nterm == off2) {
				prot_coord.Cterm = off3;
			} else if (prot_coord.Nterm == off1) {
				prot_coord.Cterm = off1;
			}
		} else if (length % 3 == 1) {
			if (prot_coord.Nterm == off3) {
				prot_coord.Cterm = off1;
			} else if (prot_coord.Nterm == off1) {
				prot_coord.Cterm = off3;
			} else if (prot_coord.Nterm == off2) {
				prot_coord.Cterm = off2;
			}
		}

		// calc protein coordinates
		if (prot_coord.Nterm != off3) {
			prot_coord.start = prev_prot_coord.end;
		} else {
			if (prev_prot_coord.end == 0 && coord_map.empty()) {
				prot_coord.start = 0;
			} else {
				prot_coord.start = prev_prot_coord.end + 1;
			}
		}

		int offsets = 0;
		if (prot_coord.Nterm != off3) {
			offsets = offsets + prot_coord.Nterm;
		}

		if (genome_coordinates.strand == fwd) {
			length = genome_coordinates.end - genome_coordinates.start + 1 - offsets;
		} else {
			length = genome_coordinates.start - genome_coordinates.end + 1 - offsets;
		}

		int pep_length = length / 3;

		int pep_end = prot_coord.start + pep_length - 1;
		if (prot_coord.Cterm != off3) {
			pep_end = pep_end + 1;
		}
		if (prot_coord.Nterm != off3) {
			pep_end = pep_end + 1;
		}

		prot_coord.end = pep_end;

		prev_prot_coord = prot_coord;

		coord_map.insert(CoordinateMapType::value_type(prot_coord, genome_coordinates));
	}
	return coord_map;
}

void PeptideEntry::add_peptide(CoordinateWrapper& coordwrapper, const std::string& sequence, const std::string& ptmSequence, const std::string& tag, unsigned int sigPSMs, transcripts_t const& transcripts, unsigned int genes, std::ofstream& ofstream, double quant) {
	if (m_sequence.compare("") == 0) {
		m_sequence = sequence;
	}
	if (m_numtranscripts == 0) {
		m_numtranscripts = transcripts.m_entries.size();
	}

	if (genes > 1) {
		m_geneunique = false;
	}
	if (transcripts.m_entries.size() > 1) {
		m_transcriptunique = false;
	}

	add_tags(tag, sigPSMs, quant);
	//iterate all found transcripts.
	for (transcripts_t::transcript_id_map_t::const_iterator it = transcripts.m_entries.begin(); it != transcripts.m_entries.end(); ++it) {
		// find all genomic coordinates
		std::vector<std::vector<GenomeCoordinates>> genomic_coordinates = coordwrapper.lookup_entry((*it).first).find_coordinates(sequence.size(), (*it).second);
		unsigned int CDS_annotation_correct = coordwrapper.lookup_entry((*it).first).is_cds_annotation_correct();
		//iterate all genomic coordinates.
		for (size_t pep_it = 0; pep_it < genomic_coordinates.size(); ++pep_it) {
			//creates PetideCoordinates
			PeptideCoordinates* pep_coord = new PeptideCoordinates(create_coordinate_map_type(genomic_coordinates.at(pep_it)), CDS_annotation_correct);
			if (genomic_coordinates.at(pep_it).size() != 0) {
				//and saves them. 
				m_pepcoordinates.insert(pep_coord);
				//sets start and end coord of the peptide entry to min/max values. these are used for comparing PeptideEntry objects.
				if (m_startcoord == 0) {
					m_startcoord = pep_coord->get_transcript_coordinates().start;
				}
				if (m_endcoord == 0) {
					m_endcoord = pep_coord->get_transcript_coordinates().end;
				}
				if (pep_coord->get_transcript_coordinates().start < m_startcoord) {
					m_startcoord = pep_coord->get_transcript_coordinates().start;
				}
				if (pep_coord->get_transcript_coordinates().end > m_endcoord) {
					m_endcoord = pep_coord->get_transcript_coordinates().end;
				}
			}
		}
	}
	add_ptm(ptmSequence);
}

void PeptideEntry::add_peptide(std::string const& ptmsequence, std::string const& tag, unsigned int sigPSMs, double quant) {
	if (m_sequence.compare(ptmsequence) != 0) {
		add_ptm(ptmsequence);
	}
	add_tags(tag, sigPSMs, quant);
}

bool peptideentry_p_compare::operator()(const PeptideEntry* lhs, const PeptideEntry* rhs) const {
	return (*lhs) < (*rhs);
}