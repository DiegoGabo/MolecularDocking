#include <iostream>
#include <vector>
#include <list>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <math.h>
#include "structures_molecule.hpp"
#include "structures_pocket.hpp"
#include "parser.hpp"

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

/*
@param a1 the first atom
@param a2 the second atom
return the euclidean distance of the atoms passed as a parameter
*/
float euclideanDistance(Atom a1, Atom a2)
{
	return sqrt(pow((a1.getX() - a2.getX()), 2) + pow((a1.getY() - a2.getY()), 2) + pow((a1.getZ() - a2.getZ()), 2));
}

/*
@param ligand the molecole 
@param pocket the pocket in which the ligand has to merge
@return a score that indicates how the ligand fix the pocket
*/
float calcolateScore(Molecule ligand, Pocket pocket)
{
	const float my_epsilon = 1.0e-16f;
	float score = 0.0f;

	for (const Atom atom_l : ligand.getAtoms())
	{
		float distance_min = 1.0e37f;

		for (const Atom atom_p : pocket.getAtoms())
		{
			float d = euclideanDistance(atom_l, atom_p);

			if (d < distance_min)
				distance_min = d;
		}

		score += distance_min;
	}

	if (score < my_epsilon)
		score = my_epsilon;

	return static_cast<float>(ligand.getAtoms().size()) / score;
}

Molecule copyMolecule(Molecule molecule)
{
	Molecule newMolecule(molecule.getName());
	for (Atom atom : molecule.getAtoms())
		newMolecule.addAtom(atom);

	std::vector<std::list<int>> links = molecule.getLinks();
	for (int src = 0; src < links.size(); src++)
	{
		for (int dest : links.at(src))
			if(src < dest)
				newMolecule.setEdge(src, dest);
	}
	return newMolecule;
}

float rotateMolecule(Molecule & molecule, Molecule & moleculeRotated, int angle, std::pair<Atom, Atom> rotamer, Pocket pocket)
{
	moleculeRotated = copyMolecule(molecule);
	std::vector<int> pointsToRotate;
	pointsToRotate = molecule.getRotamerSuccessors(rotamer);

	cout << "\n\nI want to rotate the following points of " << std::to_string(angle) << " degree:\t";
	for (int point : pointsToRotate)
		cout << std::to_string(point) << " ";

	matrix<float> rotationMatrix = createRotationMatrix(angle, rotamer.first, rotamer.second);
	cout << "\n\nThis is the rotation matrix:\n" << rotationMatrix << std::endl;

	//the for loop applies the rotation to all points in pointsToRotate
	for (int point : pointsToRotate)
		moleculeRotated.transform(rotationMatrix, point);

	cout << "\nAfter rotation:\n" << moleculeRotated.to_string();
	float score = calcolateScore(moleculeRotated, pocket);
	cout << "\nScore" << std::to_string(score);
	
	return score;
}

int main()
{
	////creation of a molecule in order to text the rotations
	Molecule molecule("TESTING MOLECULE");

	//creation of all the atoms
	Atom a1(0, 1, 0);
	Atom a2(0, -1, 0);
	Atom a3(1, 0, 0);
	Atom a4(2, 0, 0);
	Atom a5(3, 1, 0);
	Atom a6(3, -1, 0);
	Atom a7(4, 0, 0);
	Atom a8(-1, -1, 0);

	////atoms added to molecule
	molecule.addAtom(a1);
	molecule.addAtom(a2);
	molecule.addAtom(a3);
	molecule.addAtom(a4);
	molecule.addAtom(a5);
	molecule.addAtom(a6);
	molecule.addAtom(a7);
	molecule.addAtom(a8);

	////edges added to molecule
	molecule.setEdge(0, 2);
	molecule.setEdge(1, 2);
	molecule.setEdge(2, 3);
	molecule.setEdge(3, 4);
	molecule.setEdge(3, 5);
	molecule.setEdge(4, 6);
	molecule.setEdge(5, 6);
	molecule.setEdge(1, 7);

	std::cout << molecule.to_string();

	std::vector<Molecule> molecules;
	molecules.push_back(molecule);
	float bestScore = 0;
	Molecule bestMolecule("");

	//create the pocket
	Pocket pocket(5, 5, 0.2);
	pocket.transformation();
	cout << "\n\n" << pocket.to_string();

	for (Molecule molecule : molecules)
	{
		//get all the rotamers of the molecule
		std::vector<std::pair<Atom, Atom>> rotamers = molecule.getRotamers();
		int rotamerNumber = 1;
		std::cout << "\n\nList of rotamers";
		for (std::pair<Atom, Atom> rotamer : rotamers)
		{
			cout << "\nRotamer number " << to_string(rotamerNumber) << " " << molecule.getAtomIndex(rotamer.first) << " " << molecule.getAtomIndex(rotamer.second);
			rotamerNumber++;
		}

		//cycle for each rotamer of the molecule
		for (std::pair<Atom, Atom> rotamer : rotamers)
		{
			std::cout << "\n\nI Consider the rotamer " << rotamer.first.to_string() << rotamer.second.to_string();
			
			float bestLocalScore = 0;
			Molecule bestLocalMolecule(molecule.getName());

			//cicles in which all rotations are performed 
			for (int angle = 0; angle<360; angle += 120)
			{
				Molecule moleculeRotated = Molecule(molecule.getName());
				float score = rotateMolecule(molecule, moleculeRotated, angle, rotamer, pocket);
				if (score > bestLocalScore)
				{
					bestLocalMolecule = copyMolecule(moleculeRotated);
					bestLocalScore = score;
				}
			}
			
			for (int angle = 0; angle<360; angle += 120)
			{
				Molecule moleculeRotated = Molecule(bestLocalMolecule.getName());
				float score = rotateMolecule(bestLocalMolecule, moleculeRotated, angle, make_pair(rotamer.second, rotamer.first), pocket);
				if (score > bestLocalScore)
				{
					bestLocalMolecule = copyMolecule(moleculeRotated);
					bestLocalScore = score;
				}
			}

			if(bestLocalScore > bestScore)
			{
				bestMolecule = copyMolecule(bestLocalMolecule);
				bestScore = bestLocalScore;			
			}
		}
	}
	cout << "\n\nThe best score is: " << std::to_string(bestScore);
	cout << "\n\nBest molecule:\n " << bestMolecule.to_string() << std::endl;
}
