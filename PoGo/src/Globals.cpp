#include "Globals.h"

unsigned int GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::KMER_LENGTH(5);
unsigned int GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_MISMATCHES(0);
std::vector<char> GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_AMINO_ACIDS = { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'J', 'K', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' };
bool GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ONE_IN_FIVE_MODE(false);

std::string GENOME_MAPPER_GLOBALS::ID::GENE_ID("ENSG");
std::string GENOME_MAPPER_GLOBALS::ID::TRANSCRIPT_ID("ENST");
std::string GENOME_MAPPER_GLOBALS::ID::EXON_ID("ENSE");
int GENOME_MAPPER_GLOBALS::ID::LENGTH(GENE_ID.size() + 11);
std::string GENOME_MAPPER_GLOBALS::ID::GTF_GENE_ID("gene_id \"");
std::string GENOME_MAPPER_GLOBALS::ID::GTF_TRANSCRIPT_ID("transcript_id \"");
std::string GENOME_MAPPER_GLOBALS::ID::GTF_EXON_ID("exon_id \"");
std::string GENOME_MAPPER_GLOBALS::ID::FASTA_GENE_ID("gene:");
std::string GENOME_MAPPER_GLOBALS::ID::FASTA_TRANSCRIPT_ID("transcript:");

TAXONOMY_IDENTIFIERS cow = { "ENSBTAG", "ENSBTAT", "ENSBTAE", 18};
TAXONOMY_IDENTIFIERS marmoset = { "ENSCJAG", "ENSCJAT", "ENSCJAE", 18};
TAXONOMY_IDENTIFIERS dog = { "ENSCAFG", "ENSCAFT", "ENSCFAE", 18};
TAXONOMY_IDENTIFIERS vervet_AGM = { "ENSCSAG", "ENSCSAT", "ENSCSAE", 18};
TAXONOMY_IDENTIFIERS cintestinalis = { "ENSCING", "ENSCINT", "ENSCINE", 18};
TAXONOMY_IDENTIFIERS horse = { "ENSECAG", "ENSECAT", "ENSECAE", 18};
TAXONOMY_IDENTIFIERS cat = { "ENSFCAG", "ENSFCAT", "ENSFCAE", 18};
TAXONOMY_IDENTIFIERS chicken = { "ENSGALG", "ENSGALT", "ENSGALE", 18};
TAXONOMY_IDENTIFIERS gorilla = { "ENSGGOG", "ENSGGOT", "ENSGGOE", 18};
TAXONOMY_IDENTIFIERS human = { "ENSG", "ENST", "ENSE", 15};
TAXONOMY_IDENTIFIERS macaque = { "ENSMMUG", "ENSMMUT", "ENSMMUE", 18};
TAXONOMY_IDENTIFIERS turkey = { "ENSMGAG", "ENSMGAT", "ENSMGAE", 18};
TAXONOMY_IDENTIFIERS opossum = { "ENSMODG", "ENSMODT", "ENSMODE", 18};
TAXONOMY_IDENTIFIERS mouse = { "ENSMUSG", "ENSMUST", "ENSMUSE", 18};
TAXONOMY_IDENTIFIERS platypus = { "ENSOANG", "ENSOANT", "ENSOANE", 18};
TAXONOMY_IDENTIFIERS rabbit = { "ENSOCUG", "ENSOCUT", "ENSOCUE", 18};
TAXONOMY_IDENTIFIERS medaka = { "ENSORLG", "ENSORLT", "ENSORLE", 18};
TAXONOMY_IDENTIFIERS sheep = { "ENSOARG", "ENSOART", "ENSOARE", 18};
TAXONOMY_IDENTIFIERS chimp = { "ENSPTRG", "ENSPTRT", "ENSPTRE", 18};
TAXONOMY_IDENTIFIERS olivebaboon = { "ENSPANG", "ENSPANT", "ENSPANE", 18};
TAXONOMY_IDENTIFIERS orangutan = { "ENSPPYG", "ENSPPYT", "ENSPPYE", 18};
TAXONOMY_IDENTIFIERS rat = { "ENSRNOG", "ENSRNOT", "ENSRNOE", 18};
TAXONOMY_IDENTIFIERS pig = { "ENSSSCG", "ENSSSCT", "ENSSSCE", 18};
TAXONOMY_IDENTIFIERS zebrafinch = { "ENSTGUG", "ENSTGUT", "ENSTGUE", 18};
TAXONOMY_IDENTIFIERS zebrafish = { "ENSDARG", "ENSDART", "ENSDARE", 18};
TAXONOMY_IDENTIFIERS tetraodon = { "ENSTNIG", "ENSTNIT", "ENSTNIE", 18};
TAXONOMY_IDENTIFIERS yeast = { "Y", "Y", "Y", 8 };


std::map<std::string, TAXONOMY_IDENTIFIERS*> GENOME_MAPPER_GLOBALS::TAX = {
	{ "bos taurus", &cow },
	{ "cow", &cow },
	{ "9913", &cow },
	{ "callithrix jacchus", &marmoset },
	{ "marmoset", &marmoset },
	{ "9483", &marmoset },
	{ "canis lupus familiaris", &dog },
	{ "dog", &dog },
	{ "9615", &dog },
	{ "chlorocebus sabaeus", &vervet_AGM },
	{ "vervet_agm", &vervet_AGM },
	{ "60711", &vervet_AGM },
	{ "ciona intestinalis", &cintestinalis },
	{ "cintestinalis", &cintestinalis },
	{ "7719", &cintestinalis },
	{ "equus caballus", &horse },
	{ "horse", &horse },
	{ "9796", &horse },
	{ "felis catus", &cat },
	{ "cat", &cat },
	{ "9685", &cat },
	{ "gallus gallus", &chicken },
	{ "chicken", &chicken },
	{ "9031", &chicken },
	{ "gorilla gorilla gorilla", &gorilla },
	{ "gorilla", &gorilla },
	{ "9595", &gorilla },
	{ "homo sapiens", &human },
	{ "human", &human },
	{ "9606", &human },
	{ "macaca mulatta", &macaque },
	{ "macaque", &macaque },
	{ "9544", &macaque },
	{ "meleagris gallopavo", &turkey },
	{ "turkey", &turkey },
	{ "9103", &turkey },
	{ "monodelphis domestica", &opossum },
	{ "opossum", &opossum },
	{ "13616", &opossum },
	{ "mus musculus", &mouse },
	{ "mouse", &mouse },
	{ "10090", &mouse },
	{ "ornithorhynchus anatinus", &platypus },
	{ "platypus", &platypus },
	{ "9258", &platypus },
	{ "oryctolagus cuniculus", &rabbit },
	{ "rabbit", &rabbit },
	{ "9986", &rabbit },
	{ "oryzias latipes", &medaka },
	{ "medaka", &medaka },
	{ "8090", &medaka },
	{ "ovis aries", &sheep },
	{ "sheep", &sheep },
	{ "9940", &sheep },
	{ "pan troglodytes", &chimp },
	{ "chimpanzee", &chimp },
	{ "9598", &chimp },
	{ "papio anubis", &olivebaboon },
	{ "olivebaboon", &olivebaboon },
	{ "9555", &olivebaboon },
	{ "pongo abelii", &orangutan },
	{ "orangutan", &orangutan },
	{ "9601", &orangutan },
	{ "rattus norvegicus", &rat },
	{ "rat", &rat },
	{ "10116", &rat },
	{ "sus scrofa", &pig },
	{ "pig", &pig },
	{ "9823", &pig },
	{ "taeniopygia guttata", &zebrafinch },
	{ "zebra finch", &zebrafinch },
	{ "59729", &zebrafinch },
	{ "danio rerio", &zebrafish },
	{ "zebrafish", &zebrafish },
	{ "7955", &zebrafish },
	{ "tetraodon nigroviridid", &tetraodon },
	{ "tetraodon", &tetraodon },
	{ "99883", &tetraodon },
	{ "saccharomyces cerevisiae", &yeast},
	{ "yeast", &yeast },
	{ "4932", &yeast }
};
