#include "GTFParser.h"


GTFParser* GTFParser::m_instance = nullptr;

GTFParser::GTFParser() {}
GTFParser::GTFParser(const GTFParser& other) {}
GTFParser::~GTFParser() {}

GTFParser& GTFParser::operator =(const GTFParser& other) { return *this; }

GTFParser* GTFParser::get_instance() {
	if (!m_instance) {
		m_instance = new GTFParser();
	}
	return m_instance;
}

bool GTFParser::open(std::string file) {
	if (!m_ifstream.is_open()) {
		m_line = "";
		m_ifstream.open(file.c_str());
	}
	return m_ifstream.good();
}

void GTFParser::close() {
	m_line.clear();
	m_ifstream.close();
	m_ifstream.clear();
}

bool GTFParser::is_first_strand(std::vector<std::string> const& tokens) {
	return tokens.at(6).compare("+") == 0;
}

bool GTFParser::is_cds(std::vector<std::string> const& tokens) {
	return tokens.at(2).compare("CDS") == 0;
}

bool GTFParser::is_next_transcript(std::vector<std::string> const& tokens) {
	return tokens.at(2).compare("transcript") == 0;
}

bool GTFParser::is_next_gene(std::vector<std::string> const& tokens) {
	return tokens.at(2).compare("gene") == 0;
}

void GTFParser::read(const std::string& file, CoordinateWrapper& coordwrapper, MappedPeptides& mapping) {
	if (!open(file)) {
		throw GTFParser__file_not_found_exception();
	}

	ProteinEntry* p_protein_entry = nullptr;
	CoordinateMapType coordinates_map = CoordinateMapType();
	Coordinates protein_coordinates = Coordinates();
	Coordinates prev_proteint_coordinates = Coordinates();

	std::vector<std::string> tokens;
	while (std::getline(m_ifstream, m_line)) {
		if ((m_line[0] != '#')) {
			tokenize(m_line, tokens, "\t");

			if (is_next_gene(tokens)) {
				mapping.add_gene_from_gtf(m_line);
			}
			if (is_next_transcript(tokens)) {
				mapping.add_transcript_id_to_gene(m_line);
				if (p_protein_entry != nullptr) {
					p_protein_entry->set_coordinate_map(coordinates_map);
				}
				p_protein_entry = &coordwrapper.lookup_entry(GeneEntry::extract_transcript_id(m_line));
				if (p_protein_entry == nullptr) {
					std::cout << "ERROR: No entry for with transcript ID: " << GeneEntry::extract_transcript_id(m_line) << "\n";
					continue;
				}
				protein_coordinates = Coordinates();
				prev_proteint_coordinates = Coordinates();
				prev_proteint_coordinates.Cterm = off3;
				prev_proteint_coordinates.Nterm = off3;
				prev_proteint_coordinates.start = 0;
				prev_proteint_coordinates.end = 0;
				coordinates_map = CoordinateMapType();
			} else if (is_cds(tokens)) {
				GenomeCoordinates genCoord = extract_coordinates_from_gtf_line(tokens);
				protein_coordinates = Coordinates();
				// get nterm from prev exon
				if (genCoord.frame != unknown) {
					protein_coordinates.Nterm = Offset(int(genCoord.frame));
				} else {
					if (prev_proteint_coordinates.Cterm != off3) {
						protein_coordinates.Nterm = Offset(3 - prev_proteint_coordinates.Cterm);
					} else {
						protein_coordinates.Nterm = off3;
					}
				}

				int length = 0;

				if (is_first_strand(tokens)) {
					length = genCoord.end - genCoord.start + 1;
				} else if (!is_first_strand(tokens)) {
					length = genCoord.start - genCoord.end + 1;
				}

				// calc cterm
				if (length % 3 == 0) {
					if (protein_coordinates.Nterm != off3) {
						protein_coordinates.Cterm = Offset(3 - protein_coordinates.Nterm);
					} else {
						protein_coordinates.Cterm = off3;
					}
				} else if (length % 3 == 2) {
					if (protein_coordinates.Nterm == off3) {
						protein_coordinates.Cterm = off2;
					} else if (protein_coordinates.Nterm == off2) {
						protein_coordinates.Cterm = off3;
					} else if (protein_coordinates.Nterm == off1) {
						protein_coordinates.Cterm = off1;
					}
				} else if (length % 3 == 1) {
					if (protein_coordinates.Nterm == off3) {
						protein_coordinates.Cterm = off1;
					} else if (protein_coordinates.Nterm == off1) {
						protein_coordinates.Cterm = off3;
					} else if (protein_coordinates.Nterm == off2) {
						protein_coordinates.Cterm = off2;
					}
				}

				// calc protein coordinates
				if (protein_coordinates.Nterm != off3) {
					protein_coordinates.start = prev_proteint_coordinates.end;
				} else {
					if (prev_proteint_coordinates.end == 0 && coordinates_map.empty()) {
						protein_coordinates.start = 0;
					} else {
						protein_coordinates.start = prev_proteint_coordinates.end + 1;
					}
				}

				int offsets = 0;
				if (protein_coordinates.Nterm != off3) {
					offsets = offsets + protein_coordinates.Nterm;
				}

				if (is_first_strand(tokens)) {
					length = genCoord.end - genCoord.start + 1 - offsets;
				} else if (!is_first_strand(tokens)) {
					length = genCoord.start - genCoord.end + 1 - offsets;
				}

				int peplength = length / 3;

				int pepend = protein_coordinates.start + peplength - 1;
				if (protein_coordinates.Cterm != off3) {
					pepend = pepend + 1;
				}
				if (protein_coordinates.Nterm != off3) {
					pepend = pepend + 1;
				}

				protein_coordinates.end = pepend;

				prev_proteint_coordinates = protein_coordinates;

				coordinates_map.insert(CoordinateMapType::value_type(protein_coordinates, genCoord));
			}
		}
		tokens.clear();
	}
	if (p_protein_entry != nullptr) {
		p_protein_entry->set_coordinate_map(coordinates_map);
	}
	close();
}
