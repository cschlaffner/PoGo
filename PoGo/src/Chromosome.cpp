#include "Chromosome.h"

std::map<std::string, int> Chromosome::chr_to_int;
std::map<int, std::string> Chromosome::int_to_chr;
std::map<std::string, int> Chromosome::scaffold_names;

Chromosome::Chromosome(void) {
	this->name = "NA";
}

Chromosome::Chromosome(std::string const & name) {
	this->name = name;
}

Chromosome::~Chromosome(void) {

}

std::map<std::string, int>& Chromosome::getChr_to_int(void) {
	return chr_to_int;
}

std::map<int, std::string>& Chromosome::getInt_to_chr(void) {
	return int_to_chr;
}

std::map<std::string, int>& Chromosome::getScaffold_names(void) {
	return scaffold_names;
}

std::string Chromosome::forValue(int const & value) {
	std::map<int, std::string>::iterator it = getInt_to_chr().find(value);
	if (it != getInt_to_chr().end()) {
		return it->second;
	}
	return "NA";
}

int Chromosome::forName(std::string const & name) {
	std::map<std::string, int>::iterator it = getChr_to_int().find(name);
	if (it != getChr_to_int().end()) {
		return it->second;
	}
	return -1;
}

void Chromosome::addChr(std::string const & name) {
	std::string tmpname = name;
	std::string substring = tmpname.substr(0, 3);
	if (substring == "chr" || substring == "Chr") {
		tmpname = tmpname.substr(3);
	}
	std::map<std::string, int>::iterator it = getChr_to_int().find(tmpname);
	if (it == getChr_to_int().end())
	{
		it = getScaffold_names().find(tmpname);
		if (it == getScaffold_names().end()) {
			int newrank = getChr_to_int().size() + 1;
			getChr_to_int().insert(std::pair<std::string, int>(tmpname, newrank));
			getInt_to_chr().insert(std::pair<int, std::string>(newrank, tmpname));
			if (tmpname == "M") {
				newrank = newrank + 1;
				getChr_to_int().insert(std::pair<std::string, int>("MT", newrank));
				getInt_to_chr().insert(std::pair<int, std::string>(newrank, "MT"));
			}
			else if (tmpname == "MT") {
				newrank = newrank + 1;
				getChr_to_int().insert(std::pair<std::string, int>("M", newrank));
				getInt_to_chr().insert(std::pair<int, std::string>(newrank, "M"));
			}
		}
	}
}

void Chromosome::addScaffold(std::string const & name) {
	std::map<std::string, int>::iterator it = getChr_to_int().find("scaffold");
	if (it != getChr_to_int().end()) {
		int newrank = getChr_to_int().size() + 1;
		getChr_to_int().insert(std::pair<std::string, int>("scaffold", newrank));
	}
}

int Chromosome::getValue(void) const {
	std::map<std::string, int>::iterator it = getChr_to_int().find(name);
	if (it != getChr_to_int().end()) {
		return it->second;
	}
	return -1;
}

std::string Chromosome::getName(void) const {
	return name;
}

void Chromosome::setName(std::string const& name) {
	this->name = name;
}

bool Chromosome::isNA(void) const {
	return name == "NA";
}

bool Chromosome::isScaffold(void) const {
	std::map<std::string, int>::iterator it = getScaffold_names().find(name);
	return it != getScaffold_names().end();
}
