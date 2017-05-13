#include <iostream>
#include <vector>
#include <matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <math.h>
#include "basicgraph.h"

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
		return "\nx = " + std::to_string(x) + " y = " + std::to_string(y) + " z = " + std::to_string(z);
	}

	void transform(matrix<double>  transformationMatrix)
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
	std::vector<Atom> atoms;
	//std::list<Link> links;

public:
	Molecule(std::vector<Atom> atoms)
	{
		this->atoms = atoms;
	}

	string toString()
	{
		string molecule = "The molecule has the following atoms";
		for (Atom atom : atoms)
		{
			molecule += atom.toString();
		}
		return molecule;
	}
};

matrix<double> createRotationMatrix(int angle, Atom first, Atom second)
{
	matrix<double>  rotationMatrix(4, 4);
	float u, v, w, L, u2, v2, w2, theta;
	u = second.getX() - first.getX();
	v = second.getY() - first.getY();
	w = second.getZ() - first.getZ();
	u2 = u * u;
	v2 = v * v;
	w2 = w * w;
	L = u2 + v2 + w2;
	theta = angle * M_PI / 180.0;

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
	Atom a1(0, 1, 0);
	Atom a2(0, -1, 0);
	Atom a3(1, 0, 0);
	Atom a4(2, 0, 0);
	Atom a5(3, 1, 0);
	Atom a6(3, -1, 0);
	Atom a7(4, 0, 0);

	std::vector<Atom> listOfAtoms(7, a1);
	listOfAtoms.at(0) = a1;
	listOfAtoms.at(1) = a2;
	listOfAtoms.at(2) = a3;
	listOfAtoms.at(3) = a4;
	listOfAtoms.at(4) = a5;
	listOfAtoms.at(5) = a6;
	listOfAtoms.at(6) = a7;

	Molecule molecule(listOfAtoms);
	std::cout << molecule.toString();

	//cicle in which all rotations are performed
	for (int angle = 0; angle<360; angle++)
	{
		//Molecule moleculeToRotate = new Molecule(molecule);
		std::vector<Atom> pointsToRotate(3, a5);
		std::vector<Atom> pointsRotated(3, a5);
		pointsToRotate.at(0) = a5;
		pointsToRotate.at(1) = a6;
		pointsToRotate.at(2) = a7;

	std:cout << "\n\nI want to rotate the following points of " << std::to_string(angle) << " degree:";
		for (Atom point : pointsToRotate)
		{
			std::cout << point.toString();
		}

		matrix<float> rotationMatrix = createRotationMatrix(angle, a3, a4);
		std::cout << std::endl << rotationMatrix << std::endl;

		//the for loop applies the rotation to all points in pointsToRotate
		int i = 0;
		for (Atom point : pointsToRotate)
		{
			point.transform(rotationMatrix);
			pointsRotated.at(i) = point;
			i++;
		}

		std::cout << "\nAfter rotation:";
		for (Atom point : pointsRotated)
		{
			std::cout << point.toString();
		}

	}
}
