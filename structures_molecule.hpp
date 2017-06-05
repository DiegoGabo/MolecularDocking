#ifndef structures_molecule_hpp
#define structures_molecule_hpp

#include <string>
#include <vector>
#include <list>
#include <boost/numeric/ublas/matrix.hpp>



class Atom
{
private:
    
    float x, y, z;
    
public:
    
    /*
     creates an atom and initializes the 3 coordinates x, y and z
     */
    Atom(float x, float y, float z);
    
    /*
     @return the x coordinate
     */
    float getX() const;
    
    /*
     @return the y coordinate
     */
    float getY() const;
    
    /*
     @return the z coordinate
     */
    float getZ() const;

     /*
     return a string that describe the atom
     */
    std::string to_string();
    
    /*
     The coordinates of the atom are modified in function of the matrix given as parameter
     @param transformationMatrix the matrix that describes the trasformation
     */
    void transform(const boost::numeric::ublas::matrix<float> &  transformationMatrix);
};


/*
 Class in which there is a molecule. It is a graph that is composed by
 atoms which is a vector of Atom,
 links which is a vector of List of inters used to describe how atoms are linked
 */
class Molecule
{
private:
    
    std::string name;
    std::vector<Atom> atoms;
    std::vector<std::list<int> > links;
    
    /*
     @rotamer the rotamer that we want to verify if it is in a cycle
     @return a boolean that indicates if the rotamer is in a cycle
     */
    bool isRotamerInCycle(const std::pair<Atom, Atom> rotamer);
    
public:
    
    Molecule(std::string name);
    
    void setName(std::string name);
    
    std::string getName() const;
    
    /*
     Centres the molecule in its mass center
     */
    void centre();
    
    /*
     it adds an atom to the molecule
     @atom the atom that has to be added
     */
    void addAtom(const Atom atom);
    
    /*
     returns the index of vector atoms that corresponds to the atom passed as a parameter
     @param atom the atom of which we want to know the index
     @return the index of atoms in which there is the atom passed as a paramater
     */
    int getAtomIndex(const Atom atom) const;
    
    /*
     it adds an edge to the molecule
     @src the first element of the edge
     @dest the second element of the edge
     */
    void setEdge(int src, int dest);
    
    /*
     It returns the list of rotamers in the molecule and so all the non-terminal and not-in-cycle links
     @return the list of rotamer in the molecule
     */
    std::vector<std::pair<Atom, Atom>> getRotamers();
    
    /*
     @param atom the atom of which we want to know the successor
     @return the list of atoms connected to the atom passed as a parameter
     */
    std::vector<int> getSuccessor(int atom) const;
    
    /*
     @rotamer the rotamer of which we want to know the successors
     @return the successors of the rotamer
     */
    std::vector<int> getRotamerSuccessors(const std::pair<Atom, Atom> rotamer) const;
    
    /*
     	@return the list of atoms of the molecule
     */
	std::vector<Atom> getAtoms() const;
    
    /*
     	@return the list of links of the molecule
     */
	std::vector<std::list<int>> getLinks() const;

    /*
     * @param transformationMatrix the matrix that describe the transformation
     * @param index the index of the atom
     	Trasforms the point passed as a parameter
     */
    void transform(const boost::numeric::ublas::matrix<float> & transformationMatrix, int index);

    /*
     @return a string that describes the molecule
     */
    std::string to_string();
};

/*
     Creates the translation matrix of the molecule in its mass center
     @param xcm, ycm, zcm the coordinates of the mass center
     @return the translation matrix
     */
boost::numeric::ublas::matrix<float> createTranslationMatrix(float xcm, float ycm, float zcm);
#endif

