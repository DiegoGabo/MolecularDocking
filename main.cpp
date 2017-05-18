#include <iostream>
#include <vector>
#include <list>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <math.h>
#include "structures_molecule.hpp"
#include "structures_pocket.hpp"
using namespace std;
using namespace boost::numeric::ublas;
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
	////creation of a molecule in order to text the rotations
	//Molecule molecule("TESTING MOLECULE");

	////creation of all the atoms
	//Atom a1(0, 1, 0);
	//Atom a2(0, -1, 0);
	//Atom a3(1, 0, 0);
	//Atom a4(2, 0, 0);
	//Atom a5(3, 1, 0);
	//Atom a6(3, -1, 0);
	//Atom a7(4, 0, 0);
	//Atom a8(-1, -1, 0);

	////atoms added to molecule
	//molecule.addAtom(a1);
	//molecule.addAtom(a2);
	//molecule.addAtom(a3);
	//molecule.addAtom(a4);
	//molecule.addAtom(a5);
	//molecule.addAtom(a6);
	//molecule.addAtom(a7);
	//molecule.addAtom(a8);

	////edges added to molecule
	//molecule.setEdge(0, 2);
	//molecule.setEdge(1, 2);
	//molecule.setEdge(2, 3);
	//molecule.setEdge(3, 4);
	//molecule.setEdge(3, 5);
	//molecule.setEdge(4, 6);
	//molecule.setEdge(5, 6);
	//molecule.setEdge(1, 7);

	//std::cout << molecule.toString();

	std::vector<Molecule> molecules = parseFile();

	for (Molecule molecule : molecules)
	{
		//get all the rotators of the molecule
		std::vector<std::pair<Atom, Atom>> rotators = molecule.getRotators();
		int rotatorNumber = 1;
		std::cout << "\n\nList of rotators";
		for (std::pair<Atom, Atom> rotator : rotators)
		{
			cout << "\nRotator number " << to_string(rotatorNumber) << " " << molecule.getAtomIndex(rotator.first) << " " << molecule.getAtomIndex(rotator.second);
			rotatorNumber++;
		}

		//cycle for each rotator of the molecule
		for (std::pair<Atom, Atom> rotator : rotators)
		{
			std::cout << "\n\nI Consider the rotator " << rotator.first.toString() << rotator.second.toString();

			//cicle in which all rotations are performed
			for (int angle = 0; angle<360; angle += 15)
			{

				//create a copy of the molecule that has to be rotated
				Molecule moleculeRotated;
				for (Atom atom : molecule.getAtoms())
					moleculeRotated.addAtom(atom);

				std::vector<std::list<int>> links = molecule.getLinks();
				for (int src = 0; src < links.size(); src++)
				{
					for (int dest : links.at(src))
						moleculeRotated.setEdge(src, dest);
				}

				std::vector<int> pointsToRotate;
				pointsToRotate = molecule.getRotatorSuccessors(rotator);

				cout << "\n\nI want to rotate the following points of " << std::to_string(angle) << " degree:\t";
				for (int point : pointsToRotate)
					cout << std::to_string(point) << " ";

				matrix<float> rotationMatrix = createRotationMatrix(angle, rotator.first, rotator.second);
				cout << "\n\nThis is the rotation matrix:\n" << rotationMatrix << std::endl;

				//the for loop applies the rotation to all points in pointsToRotate
				for (int point : pointsToRotate)
					molecule.transform(rotationMatrix, point);

				cout << "\nAfter rotation:\n" << molecule.toString();
			}
		}
	}
}