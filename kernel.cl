#define MAX_SIZE 20
#include "structures_molecule.hpp"

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

kernel void vecadd(global Molecule * molecules) {
    const int idx = get_global_id(0);
	int i, j;
	for(i = 0; i<molecules[idx].numberOfAtoms; i++)
	{
		printf("\nAtom number %d: x = %f, y = %f, z = %f", i, molecules[idx].atoms[i].x, molecules[idx].atoms[i].y, molecules[idx].atoms[i].z);
	}

	for(i=0; i<molecules[idx].numberOfAtoms; i++)
	{
		printf("\nAtom number %d linked to:  ", i);
		for(j=0; j<molecules[idx].links[i].numberOfLinks; j++)
		{
			printf("%d ", molecules[idx].links[i].atoms[j]);
		}
	}
	Rotamer r1;
	r1.first=1; r1.second=2;
	float rotationMatrix[4][4];
	createRotationMatrix(molecules[idx], 30, r1, &rotationMatrix);
	
	printf("\n\nRotation matrix\n");
	for (i=0; i<4; i++)
	{
		for (j=0; j<4; j++)
		{
			printf("%f  ", rotationMatrix[i][j]);
		}
		printf("\n");
	}
}


