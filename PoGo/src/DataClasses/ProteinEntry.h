#ifndef PROTEINENTRY_H
#define PROTEINENTRY_H

#include "FastaParser.h"
struct position_mismatch_t;

class ProteinEntry {
public:
	//ctr /dtr
	ProteinEntry(void);
	ProteinEntry(std::string const&fastaHeader, std::string const&AAsequence);
	explicit ProteinEntry(FastaEntry const&fastaEntry);
	~ProteinEntry(void);
	//ctr/dtr

	//returns the transcript_id number of the current protein
	std::string const& get_transcript_id() const;
	//returns the gene_id number of the current protein
	std::string const& get_gene_id() const;
	//returns the sequence. the sequences are iso-sequences (I and L are converted to J) 
	std::string const& get_sequence() const;
	//returns the genomic coordinates.
	
	//mapping function. takes the positions calculated in the KmereMap and generates the genomic coordinates for all peptides of this protein.
	std::vector<std::vector<GenomeCoordinates>> find_coordinates(int peptideseqSize, std::vector<position_mismatch_t> const& positions);

	//setter for the coordinatesMap
	void set_coordinate_map(CoordinateMapType coordinatesMap);
	//getter for CDS_annotation_correct.
	unsigned int is_cds_annotation_correct() const;

private:
	//whole fasta header
	std::string m_fasta_header;
	//transcript number
	std::string m_transcript_id;
	//gene number
	std::string m_gene_id;
	//the AA sequence.
	std::string m_aa_sequence;

	//the coordinateMapType is a multimap: 
	//std::multimap <Coordinates (protein coordinates), GenomeCoordinates(corresponding genomic coordinates), Coordinates (passing this as third argument will use the Coordinates::operator() as comparator)>  
	//the first Coordinate are the coordinates of exons within the protein and the GenomeCoordinate is its corresponding location in the genome.
	CoordinateMapType m_coordinates_map;
	//check if the coding sequence is dividable by 3bp and not offset due to incomplete transcript annotation.
	unsigned int m_cds_annotation_correct;

	//sets all crucial values. 
	void init(std::string fastaHeader, std::string AAsequence);
	//QOL: delegates to void init(std::string fastaHeader, std::string AAsequence);
	void init(FastaEntry const&fastaEntry);
	//gets the transcriptId from a fasta header
	static std::string extract_transcript_id_fasta(std::string str);
	//gets the gene id from a fasta header
	static std::string extract_gene_id_fasta(std::string str);
};

#endif
