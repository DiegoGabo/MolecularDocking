#ifndef structures_molecule_hpp
#define structures_molecule_hpp

#define MAX_SIZE 100
#define MAX_LINKS 100
#define NAME_DIMENSION 13 //each name has 12 characters, one more for the termination character


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
    char name[NAME_DIMENSION];
	Atom atoms[MAX_SIZE];
	Link links[MAX_SIZE];
    int numberOfAtoms;
	
}Molecule;

typedef struct
{
	int first;
	int second;
    
}Rotamer;

void addAtom(Molecule* molecule, float x, float y, float z);
void setEdge(Molecule* molecule, int src, int dest);
#endif

