#define MAX_SIZE 20
typedef struct 
{
	float x, y, z;
}Atom;

typedef struct 
{
	Atom atoms[MAX_SIZE];
    int numberOfAtoms;
}Molecule;

kernel void vecadd(global Molecule * molecules) {
    const int idx = get_global_id(0);
    //molecules[idx].atoms[0] = 0;
}
