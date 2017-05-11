#include <iostream>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>

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
		
		string toString()
		{
			return "\nx = " + std::to_string(x) + " y = " + std::to_string(y) + " z = " + std::to_string(z);
		}

};

class Molecule
{
	std::vector<Atom> atoms;
	//std::list<Link> links;

public:
	Molecule(vector<Atom> atoms)
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

int main()
{

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

	for (int angle = 0; angle<360; angle++)
	{
		//Molecule moleculeToRotate = new Molecule(molecule);
		std::vector<Atom> pointsToRotate(3, a5);
		listOfAtoms.at(0) = a5;
		listOfAtoms.at(1) = a6;
		listOfAtoms.at(2) = a7;

		std:cout << "\n\nI want to rotate the following points of " << std::to_string(angle) << " degree";
		for (Atom point : pointsToRotate)
		{
			std::cout << point.toString();
		}
		
		
		matrix<double> transformationMatrix(4, 4);
		//create the matrix that rotate the points
		//int transformationMatrix[4][4] = createTransformationMatrix(angle, rotator);

		//apply the transformation to all points
		//for (point : pointsToRotate)
			//point.transform(transformationMatrix);
	}
}