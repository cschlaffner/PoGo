#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include "Globals.h"

class Chromosome {
public:
	//ctr / dtr
	Chromosome(void);
	Chromosome(std::string const& name);
	~Chromosome(void);
	//end ctr / dtr

	static std::string forValue(int const& value);
	static int forName(std::string const& name);
	static void addChr(std::string const& name);
	static void addScaffold(std::string const& name);

	int getValue(void) const;
	std::string getName(void) const;
	void setName(std::string const& name);
	bool isNA(void) const;
	bool isScaffold(void) const;

private:
	//holds mapping of chromosome names to order values
	static std::map<std::string, int> chr_to_int;
	//holds mapping of order values to chromosome names
	static std::map<int, std::string> int_to_chr;
	//holds scaffold names
	static std::map<std::string, int> scaffold_names;

	std::string name;

	static std::map<std::string, int>& getChr_to_int(void);
	static std::map<int, std::string>& getInt_to_chr(void);
	static std::map<std::string, int>& getScaffold_names(void);
};

#endif