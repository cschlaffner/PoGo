#include "main_functions.h"
#include "ExistingPeptides.h"

//exit codes --------------------------------
const int __GENOME_MAPPER_EXIT_HELP = 1;
const int __GENOME_MAPPER_EXIT_TOO_FEW_ARGS = 2;
const int __GENOME_MAPPER_EXIT_INVALID_ARG = 3;
//-------------------------------------------

// -------------------------
// Main
// -------------------------
int main(int argc, char* argv[]) {
	//setting all the possible command line options
	std::vector<std::pair<std::string, std::string>> param_list;
	param_list.push_back(std::pair<std::string, std::string>("-fasta", "Filepath for file containing protein sequences in FASTA format"));
	param_list.push_back(std::pair<std::string, std::string>("-gtf", "Filepath for file containing genome annotation in GTF format"));
	param_list.push_back(std::pair<std::string, std::string>("-in", "Filepaths for files containing peptide identifications in tab seperated format. (File format: four columns: SampleName\t\tPeptideSequence\t\tPSMs\tQuant)"));
	param_list.push_back(std::pair<std::string, std::string>("-merge", "Set 'true' to merge mappings from all files from input (default 'false')"));
	param_list.push_back(std::pair<std::string, std::string>("-format", "Select the output formats from gtf, gct, bed, ptmbed, all or combinations thereof separated by ',' (default all)"));
	param_list.push_back(std::pair<std::string, std::string>("-source", "Please give a source name which will be used in the second column in the output gtf file (default: PoGo)"));
	param_list.push_back(std::pair<std::string, std::string>("-mm", "Allowed mismatches (0, 1 or 2; default: 0)"));
	param_list.push_back(std::pair<std::string, std::string>("-mmmode", "Mismatch mode (true or false): if true mismatching with two mismaches will only allow 1 mismatch every kmersize (default: 5) positions. (default: false)"));
	param_list.push_back(std::pair<std::string, std::string>("-species", "Please give species using common or scientific name (default human). For a full list of supported species please go to https://github.com/cschlaffner/PoGo"));
	
	std::vector<std::pair<std::string, std::string>> param_list_back_sorted = param_list;

	//sorting the parameter list alphabetically
	std::sort(param_list.begin(), param_list.end());
	//keeping the -fasta -gtf and -in commands at the top of the list.
	std::sort(param_list_back_sorted.begin() + 3, param_list_back_sorted.end());

	//help option
	if (cmd_option_exists(argv, argv + argc, "-h")) {
		print_parameter_list(param_list);
		return __GENOME_MAPPER_EXIT_HELP;
	}

	//argc is 1 at least (exePoGo), and needs 3 parameters to work so if -fasta, -gtf and -in are specified argc is 4
	if (argc < 4) {
		std::cout << "Too few arguments! The mapper needs at least the -fasta, -gtf and -in parameters" << std::endl;
		print_parameter_list(param_list_back_sorted);
		return __GENOME_MAPPER_EXIT_TOO_FEW_ARGS;
	}

	std::string fasta_file_path("");
	std::string gtf_file_path("");
	std::vector<std::string> peptide_input_file_paths;
	bool merge = false;
	bool gtfout = true;
	bool gctout = true;
	bool bedout = true;
	bool ptmbedout = true;
	std::string source = "PoGo";

	//Receives the command line arguments in key-value-pairs
	std::vector<std::pair<std::string, std::string>> args = get_argument_list(argc, argv);

	//finds the corresponding pairs in the args.
	std::vector<std::pair<std::string, std::string>>::iterator fasta_it = std::find_if(args.begin(), args.end(), [](std::pair<std::string, std::string>const &p) -> bool { return (p.first == "-fasta"); });
	std::vector<std::pair<std::string, std::string>>::iterator gtf_it = std::find_if(args.begin(), args.end(), [](std::pair<std::string, std::string>const &p) -> bool { return (p.first == "-gtf"); });
	std::vector<std::pair<std::string, std::string>>::iterator in_it = std::find_if(args.begin(), args.end(), [](std::pair<std::string, std::string>const &p) -> bool { return (p.first == "-in"); });

	//check if the crucial parameters (-fasta -gtf and -in are set)
	bool exit_wrong_args(false);
	if (fasta_it == args.end() || ((fasta_it->second.find(".fa", fasta_it->second.size() - 3, 3) == std::string::npos) && (fasta_it->second.find(".fasta", fasta_it->second.size() - 6, 6) == std::string::npos))) {
		std::cout << "Please provide valid input for -fasta. The input filename has to end with .fa or .fasta" << std::endl;
		exit_wrong_args = true;
	}
	if (gtf_it == args.end() || (gtf_it->second.find(".gtf", gtf_it->second.size() - 4, 4) == std::string::npos)) {
		std::cout << "Please provide valid input for -gtf. The input filename has to end with .gtf" << std::endl;
		exit_wrong_args = true;
	}
	if (in_it == args.end() || (!(isInLastPosition(in_it->second, ".txt")) && !(isInLastPosition(in_it->second,".tsv")) && !(isInLastPosition(in_it->second,".pogo")))) {
		std::cout << "Please provide valid input for -in. Allowed file extentions are .txt, .tsv or .pogo (e.g. filename.txt or filename1.txt,filename2.txt)" << std::endl;
		exit_wrong_args = true;
	}
	//exiting if the user enters the wrong number or invalid parameters
	if (exit_wrong_args) {
		print_parameter_list(param_list_back_sorted);
		return __GENOME_MAPPER_EXIT_INVALID_ARG;
	}

	//parsing the input parameters.
	for (size_t i(0); i < param_list.size(); ++i) {
		std::string key = param_list[i].first;
		to_lowercase(key); //compares only lower cases, this way case sensitivity is disabled
		std::string param = "";
		if (!parameter_is_set(args, key, param)) {
			continue; //if the parameter is not set, its default value is used for the run
		}
		int pos = param.find(",");
		//The following part assigns the parsed values to the MS SMiV settings
		if (key == "-fasta") {
			fasta_file_path = param;
		} else if (key == "-gtf") {
			gtf_file_path = param;
		} else if (key == "-in") {
			if (pos == -1) {
				peptide_input_file_paths.push_back(param);
			} else {
				tokenize(param, peptide_input_file_paths, ",", true);
			}
		} else if (key == "-merge") {
			std::string par(param);
			to_lowercase(par);
			if (par[0] == 't') {
				merge = true;
			} else {
				std::cout << "merge could not be set to true. default (false) assumed" << std::endl;
			}
		} else if (key == "-format") {
			gtfout = false;
			gctout = false;
			ptmbedout = false;
			bedout = false;

			std::vector<std::string> outformats = std::vector<std::string>();
			tokenize(param, outformats, ",", true);
			if (outformats.size() != 0) {
				for (size_t format = 0; format < outformats.size(); ++format) {
					std::string value = outformats.at(format);
					to_lowercase(value);
					if (value == "gtf") { gtfout = true; }
					if (value == "gct") { gctout = true; }
					if (value == "bed") { bedout = true; }
					if (value == "ptmbed") { ptmbedout = true; }
					if (value == "all") {
						gtfout = true;
						gctout = true;
						bedout = true;
						ptmbedout = true;
					}
				}
			} else {
				gtfout = true;
				gctout = true;
				bedout = true;
				ptmbedout = true;
			}
		} else if (key == "-source") {
			source = param;
		} else if (key == "-mm") {
			unsigned int par = atoi(param.c_str());
			if (par >= 0 && par <= 2) {
				GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_MISMATCHES = par;
			} else {
				std::cout << "-mm: allowed mismatches need to be between 0 and 2" << ": default (0) assumed" << std::endl;
			}
		} else if (key == "-mmmode") {
			std::string par(param);
			to_lowercase(par);
			if (par[0] == 't') {
				if (GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_MISMATCHES > 1) {
					GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ONE_IN_FIVE_MODE = true;
				} else {
					std::cout << "-mmmode: cannot use mode with less than 2 mismatches. default (false) assumed." << std::endl;
				}
			} else {
				std::cout << "-mmmode: invalid input. default (F) assumed" << std::endl;
			}
		} else if (key == "-species") {
			std::string tmpparam = param;
			std::transform(tmpparam.begin(), tmpparam.end(), tmpparam.begin(), ::tolower);
			if (GENOME_MAPPER_GLOBALS::TAX.count(tmpparam) > 0) {
				GENOME_MAPPER_GLOBALS::ID::GENE_ID = GENOME_MAPPER_GLOBALS::TAX[tmpparam]->GENE_ID;
				GENOME_MAPPER_GLOBALS::ID::TRANSCRIPT_ID = GENOME_MAPPER_GLOBALS::TAX[tmpparam]->TRANSCRIPT_ID;
				GENOME_MAPPER_GLOBALS::ID::EXON_ID = GENOME_MAPPER_GLOBALS::TAX[tmpparam]->EXON_ID;
				GENOME_MAPPER_GLOBALS::ID::LENGTH = GENOME_MAPPER_GLOBALS::TAX[tmpparam]->LENGTH;
			}
			else {
				std::cout << "Error: Species not in list. For a full list of suppoerted species please go to https://github.com/cschlaffner/PoGo \n";
				return 1;
			}
		} else {
			std::cout << "Error: Could not assign parameter: " + key + "!\n"; //just in case of modifications of the param list
			return 1;
		}
	}

	//files cannot be merged if there is only one input file
	if (merge && peptide_input_file_paths.size() == 1) {
		std::cout << "cannot merge output files for one input file, default (-merge false) assumed" << std::endl;
		merge = false;
	}


	std::string plural_string((GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_MISMATCHES == 1) ? " mismatch" : " mismatches");

	std::cout << "Start: allowing " << GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_MISMATCHES << plural_string << std::endl
		<< "reading FASTA: " << fasta_file_path << std::endl;

	try {
		CoordinateWrapper coordinate_wrapper;
		coordinate_wrapper.read_fasta_file(fasta_file_path);

		std::cout << "Fasta done: " << coordinate_wrapper.size() << " proteins read." << std::endl
			<< "building KmerMap..." << std::endl;
		KmerMap kmer_map;
		coordinate_wrapper.add_all_proteins_to_kmer_map(kmer_map);

		std::cout << "KmerMap done: " << kmer_map.size() << " unique " << GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::KMER_LENGTH << "-mers created." << std::endl
			<< "reading GTF: " << gtf_file_path << std::endl;

		MappedPeptides mapped_peptides = MappedPeptides();
		assembly assem = GTFParser::get_instance()->read(gtf_file_path, coordinate_wrapper, mapped_peptides);
		std::cout << "GTF done!\nComputing genomic coordinates for: " << std::endl;

	
	for (size_t i(0); i < peptide_input_file_paths.size(); ++i) {
		std::cout << peptide_input_file_paths.at(i) << std::endl;
		std::string curr_input_file_path = peptide_input_file_paths.at(i);

        std::string final_peptide_path_results = curr_input_file_path;
		if(isInLastPosition(curr_input_file_path, ".txt")){
			final_peptide_path_results = removeExtensionOutput(final_peptide_path_results, ".txt");
		}else if(isInLastPosition(curr_input_file_path, ".tsv")){
            final_peptide_path_results = removeExtensionOutput(final_peptide_path_results, ".tsv");
        }

		std::vector<std::string> tokens;
		tokenize(curr_input_file_path, tokens, ".");

		std::string path6 = final_peptide_path_results+ "_unmapped.txt";

		ResultParser::read(curr_input_file_path, coordinate_wrapper, mapped_peptides, path6, kmer_map);
		std::cout << "Results done! (" << peptide_input_file_paths.at(i) << ")" << std::endl
			<< "writing output files" << std::endl;

		if (!merge) {
			//the gtf overwrites the input gtf if they are in the same folder

			std::string path4  = final_peptide_path_results + "_out.gtf";
			std::string path5  = final_peptide_path_results + ".bed";
			std::string path7  = final_peptide_path_results + ".gct";
			std::string path8  = final_peptide_path_results + "_ptm.bed";
			std::string path81 = final_peptide_path_results +"_no-ptm.bed";
			std::string path9 = "";
			std::string path10 = "";
			std::string path11 = "";
			std::string path12 = "";
			std::string path121 = "";

			if (assem == patchhaploscaff) {
				path9   =  final_peptide_path_results  + "_patch_hapl_scaff_out.gtf";
				path10  = final_peptide_path_results + "_patch_hapl_scaff.bed";
				path11  = final_peptide_path_results + "_patch_hapl_scaff.gct";
				path12  = final_peptide_path_results + "_patch_hapl_scaff_ptm.bed";
				path121 = final_peptide_path_results + "_patch_hapl_scaff_no-ptm.bed";
			}
			if (gtfout) {
				mapped_peptides.to_gtf(path4, source);
				mapped_peptides.to_gtf(path9, source, assem);
			}
			if (bedout) {
				mapped_peptides.to_bed(path5);
				mapped_peptides.to_bed(path10, assem);
			}
			if (gctout) {
				mapped_peptides.to_gct(path7);
				mapped_peptides.to_gct(path11, assem);
			}
			if (ptmbedout) {
				mapped_peptides.to_ptmbed(path8,path81);
				mapped_peptides.to_ptmbed(path12, path121, assem);
			}
			mapped_peptides.remove_all_peptides();
		}
	}
	if (merge) {
		std::vector<std::string> tokens;
		tokenize(peptide_input_file_paths.at(0), tokens, ".");

        std::string final_peptide_path_results = peptide_input_file_paths.at(0);
        if(isInLastPosition(final_peptide_path_results, ".txt")){
            final_peptide_path_results = removeExtensionOutput(final_peptide_path_results, ".txt");
        }else if(isInLastPosition(final_peptide_path_results, ".tsv")){
            final_peptide_path_results = removeExtensionOutput(final_peptide_path_results, ".tsv");
        }

		std::string path4 = final_peptide_path_results  + "_merged.gtf";
		std::string path5 = final_peptide_path_results  + "_merged.bed";
		std::string path7 = final_peptide_path_results  + "_merged.gct";
		std::string path8 = final_peptide_path_results  + "_merged_ptm.bed";
		std::string path81 = final_peptide_path_results + "_merged_no-ptm.bed";
		std::string path9 = "";
		std::string path10 = "";
		std::string path11 = "";
		std::string path12 = "";
		std::string path121 = "";

		if (assem == patchhaploscaff) {

			path9   = final_peptide_path_results   + "_patch_hapl_scaff_merged.gtf";
			path10  = final_peptide_path_results  + "_patch_hapl_scaff_merged.bed";
			path11  = final_peptide_path_results  + "_patch_hapl_scaff_merged.gct";
			path12  = final_peptide_path_results  + "_patch_hapl_scaff_merged_ptm.bed";
			path121 = final_peptide_path_results + "_patch_hapl_scaff_merged_no-ptm.bed";
		}

		if (gtfout) {
			mapped_peptides.to_gtf(path4, source);
			mapped_peptides.to_gtf(path9, source, assem);
		}
		if (bedout) {
			mapped_peptides.to_bed(path5);
			mapped_peptides.to_bed(path10, assem);
		}
		if (gctout) {
			mapped_peptides.to_gct(path7);
			mapped_peptides.to_gct(path11, assem);
		}
		if (ptmbedout) {
			mapped_peptides.to_ptmbed(path8, path81);
			mapped_peptides.to_ptmbed(path12, path121, assem);
		}
	}
	//if there is a problem with the reading of crucial files the program will end prematurely.
	} 
	catch (FastaParser__file_not_found_exception) {	std::cout << "FASTA file could not be opened." << std::endl; }
	catch (GTFParser__file_not_found_exception) { std::cout << "GTF file could not be opened." << std::endl; }
	catch (ResultParser__output_file_exception) { std::cout << "An error occured when trying to write a file" << std::endl; }
	

	std::cout << "done, cleaning up..." << std::endl;
	return 0;
}

