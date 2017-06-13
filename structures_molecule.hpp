#ifndef structures_molecule_hpp
#define structures_molecule_hpp

#include <string>

#define MAX_SIZE 50
#define MAX_LINKS 15

typedef struct 
{
	float x, y, z;
}Atom;

typedef struct
{
	int atoms[MAX_LINKS];
	int numberOfLinks;
}Link;

typedef struct 
{
	Atom atoms[MAX_SIZE];
	Link links[MAX_SIZE];
    int numberOfAtoms;
	std::string name;
}Molecule;

typedef struct
{
	int first;
	int second;
}Rotamer;

void addAtom(Molecule* molecule, float x, float y, float z);
void setEdge(Molecule* molecule, int src, int dest);
#endif

