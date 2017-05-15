#include <iostream>
#include <vector>
#include <list>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <math.h>
using namespace std;
using namespace boost::numeric::ublas;

class Atom
{
	float x, y, z;
public:
	Atom(float x, float y, float z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}
	float getX()
	{
		return x;
	}
	float getY()
	{
		return y;
	}
	float getZ()
	{
		return z;
	}
	string toString()
	{
		return "x = " + std::to_string(x) + " y = " + std::to_string(y) + " z = " + std::to_string(z);
	}

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

class Molecule
{
private:
	std::vector<Atom> atoms;
	std::vector<std::list<int>> links;
	int getAtomIndex(Atom atom)
	{
		for (int i = 0; i<atoms.size(); i++)
		{
			Atom a = atoms.at(i);
			if (a.getX() == atom.getX() && a.getX() == atom.getX() && a.getX() == atom.getX())
				return i;
		}
	}

public:

	void addAtom(Atom atom)
	{
		atoms.push_back(atom);
		std::list<int> link;
		links.push_back(link);
	}

	void setEdge(int src, int dest)
	{
		links.at(src).push_back(dest);
		links.at(dest).push_back(src);
	}

	std::vector<std::pair<Atom, Atom>> getRotators()
	{
		std::vector<std::pair<Atom, Atom>> rotators;
		std::vector<Atom>::iterator firstAtomIt = atoms.begin();
		std::vector<std::list<int>>::iterator firstListElement = links.begin();

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
					rotators.push_back(rotator);
				}
			}
			firstAtomIt += 1;
		}

		return rotators;
	}

	std::vector<int> getSuccessor(int atom)
	{
		std::vector<int> successors;
		std::list<int> link = links.at(atom);
		for (std::list<int>::iterator it = link.begin(); it != link.end(); ++it)
			successors.push_back(*it);
		return successors;
	}

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

	std::vector<Atom> getAllSuccessors(Atom firstAtom, Atom secondAtom)
	{
		int first = getAtomIndex(firstAtom);
		int second = getAtomIndex(secondAtom);

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
};

matrix<float> createRotationMatrix(int angle, Atom first, Atom second)
{
	matrix<float>  rotationMatrix(4, 4);

	//create the rotation matrix
	float u, v, w, L, u2, v2, w2, theta;
	u = second.getX() - first.getX();
	v = second.getY() - first.getY();
	w = second.getZ() - first.getZ();
	u2 = u * u;
	v2 = v * v;
	w2 = w * w;
	L = u2 + v2 + w2;
	theta = angle * M_PI / 180.0;

	//assign to each element of the matrix the proper value
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

	//atoms added to molecule
	molecule.addAtom(a1);
	molecule.addAtom(a2);
	molecule.addAtom(a3);
	molecule.addAtom(a4);
	molecule.addAtom(a5);
	molecule.addAtom(a6);
	molecule.addAtom(a7);

	//edges added to molecule
	molecule.setEdge(0, 2);
	molecule.setEdge(1, 2);
	molecule.setEdge(2, 3);
	molecule.setEdge(3, 4);
	molecule.setEdge(3, 5);
	molecule.setEdge(4, 6);
	molecule.setEdge(5, 6);

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

			pointsToRotate = molecule.getAllSuccessors(rotator.first, rotator.second);

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