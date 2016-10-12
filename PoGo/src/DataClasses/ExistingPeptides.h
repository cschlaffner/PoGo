#if !defined	EXISTING_PEPTIDES_H
#define			EXISTING_PEPTIDES_H

//#include "PeptideEntry.h"

class ExistingPeptides {
public:
	//ctr/dtr
	ExistingPeptides() {}
	~ExistingPeptides() {}

	//evaluates if a certain peptidestring is already present
	bool contains(std::string const& peptideString) const {
		return m_existing_peptides.count(peptideString) > 0;
	}

	// add an entry. creates the entry in the map if it does not already exist
	void add(std::string const& peptideString, PeptideEntry* peptideEntry) {
		if (!contains(peptideString)) {
			m_existing_peptides.insert(existing_peptides_map_t::value_type(peptideString, std::vector<PeptideEntry*>()));
		}
		m_existing_peptides[peptideString].push_back(peptideEntry);
	}

	//access to the elements that are saved in the map.
	std::vector<PeptideEntry*>& operator[](std::string const& peptideString) {
		return m_existing_peptides[peptideString];
	}

private:
	typedef std::map<std::string, std::vector<PeptideEntry*>> existing_peptides_map_t;
	//this map holds the existing pepides.
	existing_peptides_map_t m_existing_peptides;

};


#endif
