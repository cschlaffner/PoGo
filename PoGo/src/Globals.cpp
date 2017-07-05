#include "Globals.h"

unsigned int GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::KMER_LENGTH(5);
unsigned int GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_MISMATCHES(0);
std::vector<char> GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ALLOWED_AMINO_ACIDS = { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'J', 'K', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' };
bool GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::ONE_IN_FIVE_MODE(false);

std::string GENOME_MAPPER_GLOBALS::ID::GENE_ID("ENSG");
std::string GENOME_MAPPER_GLOBALS::ID::TRANSCRIPT_ID("ENST");
std::string GENOME_MAPPER_GLOBALS::ID::EXON_ID("ENSE");
int GENOME_MAPPER_GLOBALS::ID::LENGTH(GENE_ID.size() + 11);

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
TAXONOMY_IDENTIFIERS tetradon = { "ENSTNIG", "ENSTNIT", "ENSTNIE", 18};


std::map<std::string, TAXONOMY_IDENTIFIERS*> GENOME_MAPPER_GLOBALS::TAX = {
	{ "bos taurus", &cow },
	{ "cow", &cow },
	{ "callithrix jacchus", &marmoset },
	{ "marmoset", &marmoset },
	{ "canis lupus familiaris", &dog },
	{ "dog", &dog },
	{ "chlorocebus sabaeus", &vervet_AGM },
	{ "vervet_agm", &vervet_AGM },
	{ "ciona intestinalis", &cintestinalis },
	{ "cintestinalis", &cintestinalis },
	{ "equus caballus", &horse },
	{ "horse", &horse },
	{ "felis catus", &cat },
	{ "cat", &cat },
	{ "gallus gallus", &chicken },
	{ "chicken", &chicken },
	{ "gorilla gorilla gorilla", &gorilla },
	{ "gorilla", &gorilla },
	{ "homo sapiens", &human },
	{ "human", &human },
	{ "macaca mulatta", &macaque },
	{ "macaque", &macaque },
	{ "meleagris gallopavo", &turkey },
	{ "turkey", &turkey },
	{ "monodelphis domestica", &opossum },
	{ "opossum", &opossum },
	{ "mus musculus", &mouse },
	{ "mouse", &mouse },
	{ "ornithorhynchus anatinus", &platypus },
	{ "platypus", &platypus },
	{ "oryctolagus cuniculus", &rabbit },
	{ "rabbit", &rabbit },
	{ "oryzias latipes", &medaka },
	{ "medaka", &medaka },
	{ "ovis aries", &sheep },
	{ "sheep", &sheep },
	{ "pan troglodytes", &chimp },
	{ "chimp", &chimp },
	{ "papio anubis", &olivebaboon },
	{ "olivebaboon", &olivebaboon },
	{ "pongo abelii", &orangutan },
	{ "orangutan", &orangutan },
	{ "rattus norvegicus", &rat },
	{ "rat", &rat },
	{ "sus scrofa", &pig },
	{ "pig", &pig },
	{ "taeniopygia guttata", &zebrafinch },
	{ "zebra finch", &zebrafinch },
	{ "tetradon nigroviridid", &tetradon },
	{ "tetradon", &tetradon }	
};
