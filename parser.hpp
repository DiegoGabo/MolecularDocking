#ifndef parser_h
#define parser_h

#include "structures_molecule.hpp"
#include <string>

Molecule* parseFile(std::string n, int l);

int getDimension(std::string file_name);

#endif
