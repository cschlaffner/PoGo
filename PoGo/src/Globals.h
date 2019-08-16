#if !defined	GLOBALS_H
#define			GLOBALS_H
#include <string>
#include <vector>
#include <map>

//this struct and the contained structs hold global variables. some of these can be modified with parameters (e.g. -mm) 
//others can only be modified in the Globals.cpp file, and the executable has to be recompiled after that.
struct GENOME_MAPPER_GLOBALS {
	//globals for the peptide mapper.
	struct PEPTIDE_MAPPER {
		//KMER_LENGTH holds the size of the kmers in the KmerMap. 
		//the default value is 5. tests showed that any other kmersize slows down 
		//significantly if the number of mappings is high.
		//the executable has to be recompiled after changing this.
		static unsigned int			KMER_LENGTH;

		//allowed mismatches holds the number of allowed mismatches and has to be between 0 and 2. 
		//this can be modified with the -mm input parameter.
		static unsigned int			ALLOWED_MISMATCHES;

		//allowed amino acids holds the aminoa acids that the PossibleKeyGenerator will use to 
		//generate Keys. 
		//the executable has to be recompiled after changing this.
		static std::vector<char>	ALLOWED_AMINO_ACIDS;

		//toggles wheter 1 in 5 mode is on. this mode only works with 
		//two mismatches. if 1 in 5 mode is on, only one mismatch is 
		//allowed in every 5 amino acids. (this only works if KmerLength == 5)
		//can be toggled with the -mmmode switch.
		static bool					ONE_IN_FIVE_MODE;

		//toggles whether chromosomes and scaffolds are extracted from genome FASTA file and
		//used in output for order of chromosomes and separation of assembly and scaffolds.
		//if genome fasta file is provided through parameter -genome CHR_FROM_GENOME_FASTA==true
		//otherwise chromosome order is extracted from GTF and no separation of assembly and scaffold enabled.
		static bool					CHR_FROM_GENOME_FASTA;
	};

	//this struct holds information on the used gene and transcript ids.
	struct ID {
		//the gene id is the string pattern that preceedes the gene ID in the GTF file.
		//the executable has to be recompiled after changing this.
		static std::string		GTF_GENE_ID;

		//the transcript id is the string pattern that preceedes the transcript ID in the GTF file.
		//the executable has to be recompiled after changing this.
		static std::string		GTF_TRANSCRIPT_ID;

		//the exon id is the string pattern that preceedes the exon ID in the GTF file.
		//the executable has to be recompiled after changing this.
		static std::string		GTF_EXON_ID;

		//the gene id is the string pattern that preceedes the gene ID in the FASTA file.
		//the executable has to be recompiled after changing this.
		static std::string		FASTA_GENE_ID;

		//the transcript id is the string pattern that preceedes the transcript ID in the FASTA file.
		//the executable has to be recompiled after changing this.
		static std::string		FASTA_TRANSCRIPT_ID;


		static bool		ID_VERSION_INCLUDE;

	};
};

#endif
