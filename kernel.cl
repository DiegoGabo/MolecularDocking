#define MAX_SIZE 20
#include "structures_molecule.hpp"

inline void addAtom(Molecule* molecule, Atom atom)
{
	int index = molecule->numberOfAtoms;
	molecule->atoms[index]=atom;
	molecule->links[index].numberOfLinks=0;
	molecule->numberOfAtoms += 1;
}


inline void setEdge(Molecule * molecule, int src, int dest)
{
	int index = molecule->links[src].numberOfLinks;
	molecule->links[src].atoms[index]=dest;
	molecule->links[src].numberOfLinks+=1;
	
	index = molecule->links[dest].numberOfLinks;
	molecule->links[dest].atoms[index]=src;
	molecule->links[dest].numberOfLinks+=1;
}

inline void printMolecule(Molecule* molecule)
{
	int i, j;
	for(i = 0; i<molecule->numberOfAtoms; i++)
		printf("\nAtom number %d: x = %f, y = %f, z = %f", i, molecule->atoms[i].x, molecule->atoms[i].y, molecule->atoms[i].z);

	for(i=0; i<molecule->numberOfAtoms; i++)
	{
		printf("\nAtom number %d linked to:  ", i);
		for(j=0; j<molecule->links[i].numberOfLinks; j++)
			printf("%d ", molecule->links[i].atoms[j]);
	}
}

inline void createRotationMatrix(Molecule molecule, int angle, Rotamer rotamer, float (*rotationMatrix)[4])
{
	
	//creates the rotation matrix
	float u, v, w, L, u2, v2, w2, theta;
	u = molecule.atoms[rotamer.second].x - molecule.atoms[rotamer.first].x;
	v = molecule.atoms[rotamer.second].y - molecule.atoms[rotamer.first].y;
	w = molecule.atoms[rotamer.second].z - molecule.atoms[rotamer.first].z;
	u2 = u * u;
	v2 = v * v;
	w2 = w * w;
	L = u2 + v2 + w2;
	theta = angle * M_PI / 180.0;

	//assigns to each element of the matrix the proper value
	rotationMatrix[0][0] = (u2 + (v2 + w2) * cos(theta)) / L;
	rotationMatrix[0][1] = (u * v * (1 - cos(theta)) - w * sqrt(L) * sin(theta)) / L;
	rotationMatrix[0][2] = (u * w * (1 - cos(theta)) + v * sqrt(L) * sin(theta)) / L;
	rotationMatrix[0][3] = 0.0;
	rotationMatrix[1][0] = (u * v * (1 - cos(theta)) + w * sqrt(L) * sin(theta)) / L;
	rotationMatrix[1][1] = (v2 + (u2 + w2) * cos(theta)) / L;
	rotationMatrix[1][2] = (v * w * (1 - cos(theta)) - u * sqrt(L) * sin(theta)) / L;
	rotationMatrix[1][3] = 0.0;
	rotationMatrix[2][0] = (u * w * (1 - cos(theta)) - v * sqrt(L) * sin(theta)) / L;
	rotationMatrix[2][1] = (v * w * (1 - cos(theta)) + u * sqrt(L) * sin(theta)) / L;
	rotationMatrix[2][2] = (w2 + (u2 + v2) * cos(theta)) / L;
	rotationMatrix[2][3] = 0.0;
	rotationMatrix[3][0] = 0.0;
	rotationMatrix[3][1] = 0.0;
	rotationMatrix[3][2] = 0.0;
	rotationMatrix[3][3] = 1.0;
}

inline void copyMolecule(Molecule* molecule, Molecule* moleculeCopied)
{
	int i,j;
	moleculeCopied->numberOfAtoms=0;
	for (i=0; i<molecule->numberOfAtoms; i++)
	{
		addAtom(moleculeCopied, molecule->atoms[i]);
		for (j=0; j < molecule->links[i].numberOfLinks; j++)
		{
			if(i > molecule->links[i].atoms[j])
				setEdge(moleculeCopied, i, molecule->links[i].atoms[j]);
		}
	}
}

inline void trasform(Atom* atom, float (*matrix)[4])
{
	float x, y, z;
	x = atom->x;
	y = atom->y;
	z = atom->z;
	atom->x = matrix[0][0]*x + matrix[0][1]*y + matrix[0][2]*z;
	atom->y = matrix[1][0]*x + matrix[1][1]*y + matrix[1][2]*z;
	atom->z = matrix[2][0]*x + matrix[2][1]*y + matrix[2][2]*z;
}

inline float rotateMolecule(Molecule* molecule, Molecule* moleculeRotated, int angle, Rotamer rotamer, float (*rotationMatrix)[4], int *pointsToRotate, int numberOfPointsToRotate)
{
	copyMolecule(molecule, moleculeRotated);
	int i; 

	for(i=0; i<numberOfPointsToRotate; i++)
		trasform(&(moleculeRotated->atoms[pointsToRotate[i]]), rotationMatrix);

	printMolecule(moleculeRotated);
	return 0.5;
}

inline void addSuccessor(Molecule* molecule, int* successors, int atom, int* numberOfSuccessors)
{
	int i;
	for(i=0; i<molecule->links[atom].numberOfLinks; i++)
	{
		successors[*numberOfSuccessors] = molecule->links[atom].atoms[i];
		(*numberOfSuccessors) += 1;
	}
}

inline int getRotamerSuccessor(Molecule* molecule, Rotamer rotamer, int* successors)
{
	int numberOfSuccessors=0;
	int numberOfNotConsidered=0;
	int atomToBeConsidered[MAX_SIZE];
	addSuccessor(molecule, &atomToBeConsidered, rotamer.second, &numberOfNotConsidered);

	while(numberOfNotConsidered>0)
	{
		numberOfNotConsidered -= 1;
		int nextAtom = atomToBeConsidered[numberOfNotConsidered];
		int i, alreadyTaken = 0;

		for(i=0; i<numberOfSuccessors; i++)
			if(nextAtom == successors[i])
				alreadyTaken = 1;

		if(nextAtom == rotamer.first || nextAtom == rotamer.second)
			alreadyTaken = 1;

		if(alreadyTaken == 0)
		{
			successors[numberOfSuccessors] = nextAtom;
			numberOfSuccessors += 1;
			addSuccessor(molecule, &atomToBeConsidered, nextAtom, &numberOfNotConsidered);
		}
	}
	return numberOfSuccessors;
}

kernel void doAllRotation(global Molecule * molecules) {
    const int idx = get_global_id(0);
	Molecule molecule = molecules[idx];
	int i, j;
	printMolecule(&molecule);
	
	Rotamer listOfRotamer[MAX_SIZE];
	int numberOfRotamer=1;
	Rotamer r1;
	r1.first=1; r1.second=2;
	listOfRotamer[0]=r1;

	for (i=0; i<numberOfRotamer; i++)
	{
		Rotamer currentRotamer = listOfRotamer[0];
		int rotamerSuccessor[MAX_SIZE];
		int pointsToTrasform = getRotamerSuccessor(&molecule, currentRotamer, &rotamerSuccessor);
		
		printf("\n%d points to trasform:  ", pointsToTrasform);
		for(j=0; j<pointsToTrasform; j++)
			printf("%d  ", rotamerSuccessor[j]);

		int angle;
		for (angle=0; angle<360; angle+=60)
		{
			float rotationMatrix[4][4];
			createRotationMatrix(molecule, angle, r1, &rotationMatrix);

			printf("\n\nRotation matrix for a %d degree rotation\n", angle);
			for (i=0; i<4; i++)
			{
				for (j=0; j<4; j++)
				{
					printf("%f  ", rotationMatrix[i][j]);
				}
				printf("\n");
			}

			Molecule moleculeRotated;
			float score = rotateMolecule(&molecule, &moleculeRotated, angle, currentRotamer, &rotationMatrix, &rotamerSuccessor, pointsToTrasform);
		}
	}
}
