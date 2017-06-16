#define MAX_SIZE 100
#define SIZE_POCKET 5 
#include "structures_molecule.hpp"
#define DEBUG
#undef DEBUG

void addAtom(Molecule* molecule, float x, float y, float z)
{
	int index = molecule->numberOfAtoms;
	molecule->atoms[index].x = x;
	molecule->atoms[index].y = y;
	molecule->atoms[index].z = z;
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
	#ifdef DEBUG
	for(i = 0; i<molecule->numberOfAtoms; i++)
		printf("\nAtom number %d: x = %f, y = %f, z = %f", i, molecule->atoms[i].x, molecule->atoms[i].y, molecule->atoms[i].z);

	for(i=0; i<molecule->numberOfAtoms; i++)
	{
		printf("\nAtom number %d linked to:  ", i);
		for(j=0; j<molecule->links[i].numberOfLinks; j++)
			printf("%d ", molecule->links[i].atoms[j]);
	}
	#endif
}

inline void printRotationMatrix(float (*rotationMatrix)[4] , int angle)
{
	#ifdef DEBUG
	printf("\n\nRotation matrix for a %d degree rotation\n", angle);
	int i, j;
	for (i=0; i<4; i++)
	{
		for (j=0; j<4; j++)
		{
			printf("%f  ", rotationMatrix[i][j]);
		}
		printf("\n");
	}
	#endif
}

inline void createRotationMatrix(int angle, Atom first, Atom second, float (*rotationMatrix)[4])
{
	
	//creates the rotation matrix
	float u, v, w, L, u2, v2, w2, theta;
	u = second.x - first.x;
	v = second.y - first.y;
	w = second.z - first.z;
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
		addAtom(moleculeCopied, molecule->atoms[i].x, molecule->atoms[i].y, molecule->atoms[i].z);
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

inline float euclideanDistance(Atom a1, Atom a2)
{
	return sqrt(pow((a1.x - a2.x), 2) + pow((a1.y - a2.y), 2) + pow((a1.z - a2.z), 2));
}

inline float calculateScore(Molecule* molecule, Atom pocket[])
{
	const float my_epsilon = 1.0e-16f;
	float score = 0.0f;
	int i, j;
	for (i=0; i< molecule->numberOfAtoms; i++)
	{
		Atom atom_l = molecule->atoms[i]; 
		float distance_min = 1.0e37f;

		for (j=0; j<SIZE_POCKET*SIZE_POCKET ; j++)
		{
			Atom atom_p = pocket[j];
			float d = euclideanDistance(atom_l, atom_p);

			if (d < distance_min)
				distance_min = d;
		}
		score += distance_min;
	}

	if (score < my_epsilon)
		score = my_epsilon;

	return (float)molecule->numberOfAtoms / score;
}

inline float rotateMolecule(Molecule* molecule, Molecule* moleculeRotated, float (*rotationMatrix)[4], int *pointsToRotate, int numberOfPointsToRotate, Atom pocket[])
{ 
	copyMolecule(molecule, moleculeRotated);
	int i; 

	for(i=0; i<numberOfPointsToRotate; i++)
	{
		Atom* atomToTransform = &(moleculeRotated->atoms[pointsToRotate[i]]);
		trasform(atomToTransform, rotationMatrix);
	}
	return calculateScore(moleculeRotated, pocket);
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

inline int isInCycle(Molecule* molecule, int first, int second)
{

	int numberOfSuccessors=0, numberOfNotConsidered=0, i;
	int successors[MAX_SIZE];
	int atomToBeConsidered[MAX_SIZE];
	addSuccessor(molecule, &atomToBeConsidered, second, &numberOfNotConsidered);

	for(i=0; i<numberOfNotConsidered; i++)
		if(atomToBeConsidered[i]==first)
			atomToBeConsidered[i]=-1;

	while(numberOfNotConsidered>0)
	{
		numberOfNotConsidered -= 1;
		int nextAtom = atomToBeConsidered[numberOfNotConsidered];
		int alreadyTaken = 0;
	
		if(nextAtom == first)
			return 1;

		for(i=0; i<numberOfSuccessors; i++)
			if(nextAtom == successors[i])
				alreadyTaken = 1;

		if(nextAtom == -1 || nextAtom == second)
			alreadyTaken = 1;
	
		if(alreadyTaken == 0)
		{
			successors[numberOfSuccessors] = nextAtom;
			numberOfSuccessors += 1;
			addSuccessor(molecule, &atomToBeConsidered, nextAtom, &numberOfNotConsidered);
		}
	}
	return 0;
}

inline int getRotaimer(Molecule* molecule, Rotamer* listOfRotamer)
{
	int i, j, numberOfRotamer=0;
	int first, second;
	for(i=0; i<molecule->numberOfAtoms; i++)
	{
		for(j=0; j<molecule->links[i].numberOfLinks; j++)
		{

			first = i;
			second = molecule->links[i].atoms[j];
			if(first < second && molecule->links[first].numberOfLinks>1 && molecule->links[second].numberOfLinks>1 && isInCycle(molecule, first, second)==0)
			{	
				Rotamer rotamer;
				rotamer.first = first;
				rotamer.second = second;
				listOfRotamer[numberOfRotamer] = rotamer;
				numberOfRotamer += 1;
			}
		}
	}
	return numberOfRotamer;
}

inline void centre(Molecule* molecule){

	float x_cm=0;
	float y_cm=0;
	float z_cm=0;
	int i;
	
	for(i=0;i<molecule->numberOfAtoms;i++){
		
		x_cm = x_cm + molecule->atoms[i].x;
		y_cm = y_cm + molecule->atoms[i].y;
		z_cm = z_cm + molecule->atoms[i].z;
		
	}
	
	x_cm= x_cm/molecule->numberOfAtoms;
	y_cm= y_cm/molecule->numberOfAtoms;
	z_cm= z_cm/molecule->numberOfAtoms;
	
	for(i=0;i<molecule->numberOfAtoms;i++){
		
		molecule->atoms[i].x= molecule->atoms[i].x - x_cm;
		molecule->atoms[i].y= molecule->atoms[i].y - y_cm;
		molecule->atoms[i].z= molecule->atoms[i].z - z_cm;
		
	}

}

kernel void doAllRotation(global Molecule * molecules, global Atom p[], global int* bestMoleculeIndex, global float* bestScore) 
{

    const int idx = get_global_id(0);
	Molecule molecule = molecules[idx];

	int i, j, k;
	printMolecule(&molecule);

	Atom pocket[SIZE_POCKET * SIZE_POCKET];
	for (i=0; i<25; i++)
	{
		pocket[i] = p[i];
	}

	centre(&molecule);
	
	Rotamer listOfRotamer[MAX_SIZE];
	int numberOfRotamer = getRotaimer(&molecule, &listOfRotamer);
	
	#ifdef DEBUG
	printf("\nList of rotamers");
	for(i=0; i<numberOfRotamer; i++)
	{
		printf("\nRotamer %d: %d %d", i, listOfRotamer[i].first, listOfRotamer[i].second);
	}
	#endif
	
	Molecule bestMolecule;

	for (k=0; k<numberOfRotamer; k++)
	{
		Rotamer currentRotamer = listOfRotamer[k];
		int rotamerSuccessor[MAX_SIZE];
		int pointsToTrasform = getRotamerSuccessor(&molecule, currentRotamer, &rotamerSuccessor);
		
		#ifdef DEBUG
		printf("\n%d points to trasform:  ", pointsToTrasform);
		for(j=0; j<pointsToTrasform; j++)
			printf("%d  ", rotamerSuccessor[j]);
		#endif
		
		int angle;
		float bestLocalScore = 0;
        Molecule bestLocalMolecule;
		float rotationMatrix[4][4];

		for (angle=0; angle<360; angle+=60)
		{
			createRotationMatrix(angle, molecule.atoms[currentRotamer.first], molecule.atoms[currentRotamer.second], &rotationMatrix);
			printRotationMatrix(&rotationMatrix, angle);

			Molecule moleculeRotated;
			float score = rotateMolecule(&molecule, &moleculeRotated, &rotationMatrix, &rotamerSuccessor, pointsToTrasform, &pocket);
			printMolecule(&moleculeRotated);
			
			#ifdef DEBUG
			printf("\nWith score %f", score);
			#endif
			
			if(score > bestLocalScore)
			{
				copyMolecule(&moleculeRotated, &bestLocalMolecule);
				bestLocalScore = score;
			}
		}
		
		#ifdef DEBUG
		printf("\n\nNow the best molecule is ");
		printMolecule(&bestLocalMolecule);
		printf("\nWith best local score %f", bestLocalScore);
		#endif
		
		Rotamer oppositeRotamer;
		oppositeRotamer.first = currentRotamer.second;
		oppositeRotamer.second = currentRotamer.first;
		Molecule secondStepMolecule;
		copyMolecule(&bestLocalMolecule, &secondStepMolecule);

		pointsToTrasform = getRotamerSuccessor(&bestLocalMolecule, oppositeRotamer, &rotamerSuccessor);
		
		#ifdef DEBUG
		printf("\n%d points to trasform:  ", pointsToTrasform);
		for(j=0; j<pointsToTrasform; j++)
			printf("%d  ", rotamerSuccessor[j]);
		#endif
		
		for (angle=0; angle<360; angle+=60)
		{
			createRotationMatrix(angle, molecule.atoms[currentRotamer.second], molecule.atoms[currentRotamer.first], &rotationMatrix);
			printRotationMatrix(&rotationMatrix, angle);

			Molecule moleculeRotated;

			float score = rotateMolecule(&secondStepMolecule, &moleculeRotated, &rotationMatrix, &rotamerSuccessor, pointsToTrasform, &pocket);

			#ifdef DEBUG
			printMolecule(&moleculeRotated);
			printf("\nWith score %f", score);
			#endif
			
			if(score > bestLocalScore)
			{
				copyMolecule(&moleculeRotated, &bestLocalMolecule);
				bestLocalScore = score;
			}
		}
		#ifdef DEBUG
		printf("\n\nAnd now the best molecule is ");
		printMolecule(&bestLocalMolecule);
		printf("\nWith best local score %f", bestLocalScore);
		#endif
		
		if(bestLocalScore > bestScore[0])
		{
			copyMolecule(&bestLocalMolecule, &bestMolecule);
			bestScore[0] = bestLocalScore;
			bestMoleculeIndex = idx;
		}
	}
}
