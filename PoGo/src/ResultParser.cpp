#include "ResultParser.h"

ResultParser::ResultParser(void) {}
ResultParser::~ResultParser(void) {}

//this function will set the wheels in motion to find the peptides in the proteins.
void ResultParser::read(std::string file, CoordinateWrapper& coordwrapper, MappedPeptides& mapping, std::string unmappedoutput, KmerMap& k) {
	std::vector <std::string>	tokens;
	std::ifstream	ifs(file.c_str());
	if (ifs.good()) {
		std::ofstream	ofs(unmappedoutput.c_str());
		if (!ofs.good()) {
			ifs.close();
			throw ResultParser__output_file_exception();
		}
		std::string		line;
		std::string		peptide_string;
		std::string		tissue;
		std::string		iso_seq_without_ptms;

		unsigned int	sig_PSMs;
		double			quant;
		gene_id_map_t gene_id_map;

		while (std::getline(ifs, line)) {
			if (line.substr(0, 10).compare("Experiment") != 0 && line.substr(0, 6).compare("Sample") != 0) {
				tokenize(line, tokens, "\t", false);
				//using only the tokens needed.
				tissue = tokens.at(0);
				peptide_string = tokens.at(1);
				//std::cout << "Mapping following peptide: " << peptide_string << std::endl;
				sig_PSMs = atoi(tokens.at(2).c_str());
				quant = atof(tokens.at(3).c_str());

				//clearing the tokens vector for the next iteration.
				tokens.clear();

				if (sig_PSMs > 0) {
					//the matching will only use the amino acids.
					iso_seq_without_ptms = make_iso_sequence(remove_ptms(peptide_string));

					if (!coordwrapper.peptide_already_exists(iso_seq_without_ptms)) {
						//the gene_id_map.find_peptide function will match the peptide.
						gene_id_map = k.find_peptide(iso_seq_without_ptms);
						for (gene_id_map_t::iterator it = gene_id_map.begin(); it != gene_id_map.end(); ++it) {
							mapping.add_peptide(coordwrapper, peptide_string, tissue, sig_PSMs, gene_id_map.size(), ofs, quant, *it);
						}
					} else {
						//if the peptide already exists its genomic coordinates dont have to be recalculated.
						//only the tags and PTMs have to be added
						std::vector<PeptideEntry*>& refVec = coordwrapper.get_existing_peptides_at(iso_seq_without_ptms);
						for (std::vector<PeptideEntry*>::iterator it = refVec.begin(); it != refVec.end(); ++it) {
							(*it)->add_peptide(peptide_string, tissue, sig_PSMs, quant);
						}
					}
				}
			}
		}
		ofs.close();
	}
	ifs.close();

}
