#include <iostream>
#include <vector>
#include <list>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <math.h>

using namespace std;
using namespace boost::numeric::ublas;

/*
Class in which there is an atom. It is composed by the 3 coordinates x, y and z
*/
class Atom
{
	private:
		float x, y, z;

	public:

		/*
		creates an atom and initializes the 3 coordinates x, y and z
		*/
		Atom(float x, float y, float z)
		{
			this->x = x;
			this->y = y;
			this->z = z;
		}

		/*
		@return the x coordinate
		*/
		float getX()
		{
			return x;
		}

		/*
		@return the y coordinate
		*/
		float getY()
		{
			return y;
		}

		/*
		@return the z coordinate
		*/
		float getZ()
		{
			return z;
		}

		/*
		return a string that describe the atom
		*/
		string toString()
		{
			return "x = " + std::to_string(x) + " y = " + std::to_string(y) + " z = " + std::to_string(z);
		}

		/*
		The coordinates of the atom are modified in function of the matrix given as parameter
		@param transformationMatrix the matrix that describes the trasformation
		*/
		void transform(matrix<float>  transformationMatrix)
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
};

/*
Class in which there is a molecule. It is a graph that is composed by
atoms which is a vector of Atom,
links which is a vector of List of inters used to describe how atoms are linked
*/
class Molecule
{
	private:
		std::vector<Atom> atoms;
		std::vector<std::list<int>> links;

		/*
		returns the index of vector atoms that corresponds to the atom passed as a parameter
		@param atom the atom of which we want to know the index
		@return the index of atoms in which there is the atom passed as a paramater
		*/
		int getAtomIndex(Atom atom)
		{
			for (int i = 0; i<atoms.size(); i++)
			{
				Atom a = atoms.at(i);
				if (a.getX() == atom.getX() && a.getY() == atom.getY() && a.getZ() == atom.getZ())
					return i;
			}
		}

		/*
		@rotator the rotator that we want to verify if it is in a cycle
		@return a boolean that indicates if the rotator is in a cycle
		*/
		bool isRotatorInCycle(std::pair<Atom, Atom> rotator)
		{
			int first = getAtomIndex(rotator.first);
			int second = getAtomIndex(rotator.second);

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

	public:

		/*
		it adds an atom to the molecule
		@atom the atom that has to be added
		*/
		void addAtom(Atom atom)
		{
			atoms.push_back(atom);
			std::list<int> link;
			links.push_back(link);
		}

		/*
		it adds an edge to the molecule
		@src the first element of the edge
		@dest the second element of the edge
		*/
		void setEdge(int src, int dest)
		{
			links.at(src).push_back(dest);
			links.at(dest).push_back(src);
		}

		/*
		It returns the list of rotators in the molecule and so all the non-terminal and not-in-cycle links
		@return the list of rotator in the molecule
		*/
		std::vector<std::pair<Atom, Atom>> getRotators()
		{
			std::vector<std::pair<Atom, Atom>> rotatorsWithCycles;
			std::vector<std::pair<Atom, Atom>> rotators;
			std::vector<Atom>::iterator firstAtomIt = atoms.begin();
			std::vector<std::list<int>>::iterator firstListElement = links.begin();

			//cycle that identifies all non-terminal links
			for (std::list<int> link : links)
			{
				for (std::list<int>::iterator it = link.begin(); it != link.end(); ++it)
				{
					std::vector<Atom>::iterator secondAtomIt = atoms.begin();

					bool isRotator = (link.size() > 1) && ((*(firstListElement + *it)).size() > 1);
					if (isRotator)/*link non finali e non in un ciclo*/
					{
						secondAtomIt += *it;
						std::pair<Atom, Atom> rotator = std::make_pair(*firstAtomIt, *secondAtomIt);
						rotatorsWithCycles.push_back(rotator);
					}
				}
				firstAtomIt += 1;
			}

			//cycles that identifies all rotator deleting those that are in a cycle
			for (std::pair<Atom, Atom> rotator : rotatorsWithCycles)
			{
				if (!isRotatorInCycle(rotator))
					rotators.push_back(rotator);
			}

			return rotators;
		}

		/*
		@param atom the atom of which we want to know the successor 
		@return the list of atoms connected to the atom passed as a parameter
		*/
		std::vector<int> getSuccessor(int atom)
		{
			std::vector<int> successors;
			std::list<int> link = links.at(atom);

			for (std::list<int>::iterator it = link.begin(); it != link.end(); ++it)
				successors.push_back(*it);

			return successors;
		}

		/*
		@rotator the rotator of which we want to know the successors
		@return the successors of the rotator
		*/
		std::vector<Atom> getRotatorSuccessors(std::pair<Atom, Atom> rotator)
		{
			int first = getAtomIndex(rotator.first);
			int second = getAtomIndex(rotator.second);

			std::vector<int> atomsToBeConsidered;
			std::vector<int> successorsIndex;
			std::vector<Atom> successors;

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

			for (int a : successorsIndex)
				successors.push_back(atoms.at(a));

			return successors;
		}

		/*
		@return a string that describes the molecule
		*/
		string toString()
		{
			//creates a string with the list of atoms
			string molecule = "The molecule has the following atoms";
			for (Atom atom : atoms)
			{
				molecule += "\n";
				molecule += atom.toString();
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
};

/*
@param angle the angle that specify the rotation
@param first the first point that identify the ax of the rotation
@param second the second point that identify the ax of the rotation
@return the matrix that describes the rotation
*/
matrix<float> createRotationMatrix(int angle, Atom first, Atom second)
{
	matrix<float>  rotationMatrix(4, 4);

	//creates the rotation matrix
	float u, v, w, L, u2, v2, w2, theta;
	u = second.getX() - first.getX();
	v = second.getY() - first.getY();
	w = second.getZ() - first.getZ();
	u2 = u * u;
	v2 = v * v;
	w2 = w * w;
	L = u2 + v2 + w2;
	theta = angle * M_PI / 180.0;

	//assigns to each element of the matrix the proper value
	rotationMatrix(0, 0) = (u2 + (v2 + w2) * cos(theta)) / L;
	rotationMatrix(0, 1) = (u * v * (1 - cos(theta)) - w * sqrt(L) * sin(theta)) / L;
	rotationMatrix(0, 2) = (u * w * (1 - cos(theta)) + v * sqrt(L) * sin(theta)) / L;
	rotationMatrix(0, 3) = 0.0;
	rotationMatrix(1, 0) = (u * v * (1 - cos(theta)) + w * sqrt(L) * sin(theta)) / L;
	rotationMatrix(1, 1) = (v2 + (u2 + w2) * cos(theta)) / L;
	rotationMatrix(1, 2) = (v * w * (1 - cos(theta)) - u * sqrt(L) * sin(theta)) / L;
	rotationMatrix(1, 3) = 0.0;
	rotationMatrix(2, 0) = (u * w * (1 - cos(theta)) - v * sqrt(L) * sin(theta)) / L;
	rotationMatrix(2, 1) = (v * w * (1 - cos(theta)) + u * sqrt(L) * sin(theta)) / L;
	rotationMatrix(2, 2) = (w2 + (u2 + v2) * cos(theta)) / L;
	rotationMatrix(2, 3) = 0.0;
	rotationMatrix(3, 0) = 0.0;
	rotationMatrix(3, 1) = 0.0;
	rotationMatrix(3, 2) = 0.0;
	rotationMatrix(3, 3) = 1.0;

	return rotationMatrix;
}


int main()
{

	//creation of a molecule in order to text the rotations
	Molecule molecule;

	//creation of all the atoms
	Atom a1(0, 1, 0);
	Atom a2(0, -1, 0);
	Atom a3(1, 0, 0);
	Atom a4(2, 0, 0);
	Atom a5(3, 1, 0);
	Atom a6(3, -1, 0);
	Atom a7(4, 0, 0);
	Atom a8(-1, -1, 0);

	//atoms added to molecule
	molecule.addAtom(a1);
	molecule.addAtom(a2);
	molecule.addAtom(a3);
	molecule.addAtom(a4);
	molecule.addAtom(a5);
	molecule.addAtom(a6);
	molecule.addAtom(a7);
	molecule.addAtom(a8);

	//edges added to molecule
	molecule.setEdge(0, 2);
	molecule.setEdge(1, 2);
	molecule.setEdge(2, 3);
	molecule.setEdge(3, 4);
	molecule.setEdge(3, 5);
	molecule.setEdge(4, 6);
	molecule.setEdge(5, 6);
	molecule.setEdge(1, 7);

	std::cout << molecule.toString();

	//get all the rotators of the molecule
	std::vector<std::pair<Atom, Atom>> rotators = molecule.getRotators();
	int rotatorNumber = 1;
	std::cout << "\n\nList of rotators";
	for (std::pair<Atom, Atom> rotator : rotators)
	{
		cout << "\nRotator number " << to_string(rotatorNumber) << " " << rotator.first.toString() << rotator.second.toString();
		rotatorNumber++;
	}

	for (std::pair<Atom, Atom> rotator : rotators)
	{
		std::cout << "\n\nI Consider the rotator " << rotator.first.toString() << rotator.second.toString();
		//cicle in which all rotations are performed
		for (int angle = 0; angle<360; angle += 15)
		{
			//Molecule moleculeToRotate = new Molecule(molecule);
			std::vector<Atom> pointsToRotate;
			std::vector<Atom> pointsRotated;

			pointsToRotate = molecule.getRotatorSuccessors(rotator);

			cout << "\nI want to rotate the following points of " << std::to_string(angle) << " degree:";

			for (Atom point : pointsToRotate)
			{
				std::cout << std::endl;
				std::cout << point.toString();
			}

			matrix<float> rotationMatrix = createRotationMatrix(angle, a3, a4);
			cout << std::endl << rotationMatrix << std::endl;

			//the for loop applies the rotation to all points in pointsToRotate
			for (Atom point : pointsToRotate)
			{
				point.transform(rotationMatrix);
				pointsRotated.push_back(point);
			}

			cout << "\nAfter rotation:";
			for (Atom point : pointsRotated)
			{
				cout << point.toString();
			}
		}
	}
}