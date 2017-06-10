#include "structures_molecule.hpp"
void addAtom(Molecule* molecule, Atom atom)
{
	int index = molecule->numberOfAtoms;
	molecule->atoms[index]=atom;
	molecule->links[index].numberOfLinks=0;
	molecule->numberOfAtoms += 1;
}


void setEdge(Molecule * molecule, int src, int dest)
{
	int index = molecule->links[src].numberOfLinks;
	molecule->links[src].atoms[index]=dest;
	molecule->links[src].numberOfLinks+=1;
	
	index = molecule->links[dest].numberOfLinks;
	molecule->links[dest].atoms[index]=src;
	molecule->links[dest].numberOfLinks+=1;
}


