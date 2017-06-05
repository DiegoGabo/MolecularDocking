#include "structures_molecule.hpp"
#include <string>
#include <vector>
#include <list>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace boost::numeric::ublas;


/*
Implementations of class Atom functions
*/
Atom::Atom(float x, float y, float z)
{
    this->x = x;
    this->y = y;
    this->z = z;
}

/*
 @return the x coordinate
 */
float Atom::getX() const
{
    return x;
}

/*
 @return the y coordinate
 */
float Atom::getY() const
{
    return y;
}

/*
 @return the z coordinate
 */
float Atom::getZ() const
{
    return z;
}

/*
return a string that describe the atom
*/
string Atom::to_string()
{
    return "x = " + std::to_string(x) + " y = " + std::to_string(y) + " z = " + std::to_string(z);
}

/*
The coordinates of the atom are modified in function of the matrix given as parameter
@param transformationMatrix the matrix that describes the trasformation
*/
void Atom::transform(const matrix<float> & transformationMatrix)
{
    //transform the original vector in omogeneous coordinates in order to do the transformation
    matrix<float> homogeneusCoordinatesPoint(4, 1);
    homogeneusCoordinatesPoint(0, 0) = x;
    homogeneusCoordinatesPoint(1, 0) = y;
    homogeneusCoordinatesPoint(2, 0) = z;
    homogeneusCoordinatesPoint(3, 0) = 1;
    
    homogeneusCoordinatesPoint = prod(transformationMatrix, homogeneusCoordinatesPoint);
    
    this->x = homogeneusCoordinatesPoint(0, 0) / homogeneusCoordinatesPoint(3, 0);
    this->y = homogeneusCoordinatesPoint(1, 0) / homogeneusCoordinatesPoint(3, 0);
    this->z = homogeneusCoordinatesPoint(2, 0) / homogeneusCoordinatesPoint(3, 0);
}

/*
 Centres the molecule in its mass center
 */
void Molecule::centre(){
	
	float x_cm=0;
	float y_cm=0;
	float z_cm=0;
	
	for(Atom atom: molecule.getAtoms()){
		
		x_cm += atom.getX();
		y_cm += atom.getY();
		z_cm += atom.getZ();
		
	}
	
	x_cm= x_cm/molecule.getAtoms().size();
	y_cm= y_cm/molecule.getAtoms().size();
	z_cm= z_cm/molecule.getAtoms().size();
	
	matrix<float> translationMatrix= createTranslationMatrix(x_cm, y_cm, z_cm);
	
	for(Atom atom: molecule.getAtoms()){
		
		atom.transform(translationMatrix);
	}

}

/*
 Creates the translation matrix of the molecule in its mass center
 @param xcm, ycm, zcm the coordinates of the mass center
 @return the translation matrix
 */
matrix<float> createTranslationMatrix(float xcm, float ycm, float zcm)
{
	matrix<float>  translationMatrix(4, 4);

	translationMatrix(0,0) = 1.0;
	translationMatrix(0,1) = 0.0;
	translationMatrix(0,2) = 0.0;
	translationMatrix(0,3) = -xcm;
	translationMatrix(1,0) = 0.0;
	translationMatrix(1,1) = 1.0;
	translationMatrix(1,2) = 0.0;
	translationMatrix(1,3) = -ycm;
	translationMatrix(2,0) = 0.0;
	translationMatrix(2,1) = 0.0;
	translationMatrix(2,2) = 1.0;
	translationMatrix(2,3) = -zcm;
	translationMatrix(3,0) = 0.0;
	translationMatrix(3,1) = 0.0;
	translationMatrix(3,2) = 0.0;
	translationMatrix(3,3) = 1.0;

	return translationMatrix;
}

/*
Implementations of class Molecule functions
*/
int Molecule::getAtomIndex(const Atom atom) const
{
    for (int i = 0; i<atoms.size(); i++)
    {
        Atom a = atoms.at(i);
        
        if (a.getX() == atom.getX() && a.getY() == atom.getY() && a.getZ() == atom.getZ())
            return i;
    }
    return 0;
}

/*
@param rotamer the rotamer that we want to verify if it is in a cycle
@return a boolean that indicates if the rotamer is in a cycle
*/
bool Molecule::Molecule::isRotamerInCycle(const std::pair<Atom, Atom> rotamer)
{
    int first = getAtomIndex(rotamer.first);
    int second = getAtomIndex(rotamer.second);
    
    std::vector<int> atomsToBeConsidered;
    std::vector<int> successorsIndex;
    
    for (int a : getSuccessor(second))
        if (a != first && a != second)
            atomsToBeConsidered.push_back(a);
    
    while (!atomsToBeConsidered.empty())
    {
        int nextAtom = atomsToBeConsidered.back();
        atomsToBeConsidered.pop_back();
        if (std::find(successorsIndex.begin(), successorsIndex.end(), nextAtom) == successorsIndex.end() && nextAtom != second)
        {
            successorsIndex.push_back(nextAtom);
            for (int a : getSuccessor(nextAtom))
                atomsToBeConsidered.push_back(a);
        }
    }
    for (int successor : successorsIndex)
        if (successor == first)
            return true;
    return false;
}

/*
creates a new molecule
@param name the name of the molecule
*/
Molecule::Molecule(string name)
{
    this->name = name;
    atoms.clear();
    links.clear();
}

/*
sets the name of the molecule
@param name the name of the molecule
*/
void Molecule::setName(string name)
{
    this->name = name;
}

/*
gets the name of the molecule
@return the name of the molecule
*/
string Molecule::getName() const
{
    return name;
}

/*
 it adds an atom to the molecule
 @param atom the atom that has to be added
 */
void Molecule::addAtom(const Atom atom)
{
    atoms.push_back(atom);
    list<int> link;
    links.push_back(link);
}

/*
 it adds an edge to the molecule
 @param src the first element of the edge
 @param dest the second element of the edge
 */
void Molecule::setEdge(int src, int dest)
{
    links.at(src).push_back(dest);
    links.at(dest).push_back(src);
}

/*
 It returns the list of rotamers in the molecule and so all the non-terminal and not-in-cycle links
 @return the list of rotamer in the molecule
 */
std::vector<pair<Atom, Atom>> Molecule::getRotamers() 
{
    std::vector<pair<Atom, Atom>> rotamersWithCycles;
    std::vector<pair<Atom, Atom>> rotamers;
    std::vector<Atom>::iterator firstAtomIt = atoms.begin();
    std::vector<list<int> >::iterator firstListElement = links.begin();
    
    //cycle that identifies all non-terminal links
    for (list<int> link : links)
    {
        for (list<int>::iterator it = link.begin(); it != link.end(); ++it)
        {
            std::vector<Atom>::iterator secondAtomIt = atoms.begin();
            
            bool isRotamer = (link.size() > 1) && ((*(firstListElement + *it)).size() > 1);
            if (isRotamer)/*link non finali e non in un ciclo*/
            {
                secondAtomIt += *it;
                pair<Atom, Atom> rotamer = make_pair(*firstAtomIt, *secondAtomIt);
                rotamersWithCycles.push_back(rotamer);
            }
        }
        firstAtomIt += 1;
    }
    
    //cycles that identifies all rotamer deleting those that are in a cycle
    for (pair<Atom, Atom> rotamer : rotamersWithCycles)
    {
        if (!isRotamerInCycle(rotamer) && getAtomIndex(rotamer.first) < getAtomIndex(rotamer.second))
            rotamers.push_back(rotamer);
    }
    
    return rotamers;
}

/*
 @param atom the atom of which we want to know the successor
 @return the list of atoms connected to the atom passed as a parameter
 */
std::vector<int> Molecule::getSuccessor(int atom) const
{
    std::vector<int> successors;
    list<int> link = links.at(atom);
    
    for (list<int>::iterator it = link.begin(); it != link.end(); ++it)
        successors.push_back(*it);
    
    return successors;
}

/*
@rotamer the rotamer of which we want to know the successors
@return the successors of the rotamer
*/
std::vector<int> Molecule::getRotamerSuccessors(const std::pair<Atom, Atom> rotamer) const
{
	int first = getAtomIndex(rotamer.first);
	int second = getAtomIndex(rotamer.second);

	std::vector<int> atomsToBeConsidered;
	std::vector<int> successorsIndex;

	for (int a : getSuccessor(second))
		atomsToBeConsidered.push_back(a);

	while (!atomsToBeConsidered.empty())
	{
		int nextAtom = atomsToBeConsidered.back();
		atomsToBeConsidered.pop_back();
		if (std::find(successorsIndex.begin(), successorsIndex.end(), nextAtom) == successorsIndex.end() &&
			nextAtom != first && nextAtom != second)
		{
			successorsIndex.push_back(nextAtom);
			for (int a : getSuccessor(nextAtom))
				atomsToBeConsidered.push_back(a);
		}
	}
	return successorsIndex;
}

/*
gets the atoms of the molecule
@return the atoms of the molecule
*/
std::vector<Atom> Molecule::getAtoms() const
{
	return atoms;
}

/*
gets the links of the molecule
@return the links of the molecule
*/
std::vector<std::list<int>> Molecule::getLinks() const
{
	return links;
}

/*
@param transformationMatrix the matrix that describe the transformation
@param index the index of the atom
Trasforms the point passed as a parameter
*/
void Molecule::transform(const matrix<float> & transformationMatrix, int index)
{
	atoms.at(index).transform(transformationMatrix);
}

/*
@return a string that describes the molecule
*/
string Molecule::to_string()
{
    //creates a string with the list of atoms
    string molecule = "The molecule has the following atoms";
    for (Atom atom : atoms)
    {
        molecule += "\n";
        molecule += atom.to_string();
    }
    
    //creates a string with the structure of the graph
    molecule += "\n\nGraph structure:";
    for (int atom = 0; atom < links.size(); atom++)
    {
        std::vector<int> successors = getSuccessor(atom);
        molecule += "\nAtom ";
        molecule += std::to_string(atom);
        molecule += " linked to: ";
        for (int succ : successors)
        {
            molecule += std::to_string(succ);
            molecule += ", ";
        }
    }
    return molecule;
}

