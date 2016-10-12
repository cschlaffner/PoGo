#include "PeptideMapperUtils.h"

position_mismatch_t::position_mismatch_t(int posInProtein, int firstMismatch, int secondMismatch)
	: m_position_in_protein(posInProtein), m_first_mismatch_positon(firstMismatch), m_second_mismatch_positon(secondMismatch) {}

int position_mismatch_t::position_in_protein() const {
	return m_position_in_protein;
}

int position_mismatch_t::first() const {
	return m_first_mismatch_positon;
}

int position_mismatch_t::second() const {
	return m_second_mismatch_positon;
}
