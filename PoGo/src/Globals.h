#if !defined	GLOBALS_H
#define			GLOBALS_H
#include <string>
#include <vector>
#include <map>


struct TAXONOMY_IDENTIFIERS {
	std::string	GENE_ID;
	std::string	TRANSCRIPT_ID;
	std::string	EXON_ID;
	int			LENGTH;
	std::string GTF_GENE_ID;
	std::string GTF_TRANSCRIPT_ID;
	std::string GTF_EXON_ID;
	std::string FASTA_GENE_ID;
	std::string FASTA_TRANSCRIPT_ID;
};

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
	};

	//this struct holds information on the used gene and transcript ids.
	struct ID {
		//the gene id is the prefix (e. g. ENSG or ENSMUSG) used for the gene identifier.
		//the executable has to be recompiled after changing this.
		static std::string		GENE_ID;
		
		//the transcript id is the prefix (e. g. ENST or ENSMUST) used for the transcript identifier.
		//the executable has to be recompiled after changing this.
		static std::string		TRANSCRIPT_ID;

		//the exon id is the prefix (e. g. ENSE or ENSMUSE) used for the exon identifier.
		//the executable has to be recompiled after changing this.
		static std::string		EXON_ID;
		
		//LENGTH holds the combined length of the prefix (see above) and the number. 
		//(a default length of 11 is assumed for the number making for example ensembl numbers 15 characters long. (ENSG+11))
		//the executable has to be recompiled after changing this.
		static int				LENGTH;
		
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
		
	};

	static std::map<std::string, TAXONOMY_IDENTIFIERS*> TAX;
};

#endif
