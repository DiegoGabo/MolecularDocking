#include <iostream>
#include <vector>
#include <list>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/program_options.hpp>
#include <math.h>
#include "structures_molecule.hpp"
#include "structures_pocket.hpp"
#include "parser.hpp"
#include <time.h>

#define DB_DIMENSION 3961

#define DEBUG
#undef DEBUG

using namespace std;
using namespace boost::numeric::ublas;
namespace po = boost::program_options;


/*
@param angle the angle that specify the rotation
@param first the first point that identify the ax of the rotation
@param second the second point that identify the ax of the rotation
@return the matrix that describes the rotation
*/
matrix<float> createRotationMatrix(int angle, const Atom first, const Atom second)
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
float euclideanDistance(const Atom a1, const Atom a2)
{
	return sqrt(pow((a1.getX() - a2.getX()), 2) + pow((a1.getY() - a2.getY()), 2) + pow((a1.getZ() - a2.getZ()), 2));
}

/*
@param ligand the molecole 
@param pocket the pocket in which the ligand has to merge
@return a score that indicates how the ligand fix the pocket
*/
float calcolateScore(const Molecule & ligand, const Pocket & pocket)
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

/*
@param the molecule that has to be copied
@the new molecule created
*/
Molecule copyMolecule(const Molecule & molecule)
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

/*
@molecule the molecule that has to be rotated
@moleculeRotated the molecule after the rotation
@angle the angle of the rotation
@rotamer the ax or the rotation
@pocket the pocket necessary for calculating the score
@score the score of the rotation
*/
float rotateMolecule(const Molecule & molecule, Molecule & moleculeRotated, int angle, const std::pair<Atom, Atom> rotamer, const Pocket & pocket)
{
	moleculeRotated = copyMolecule(molecule);
	std::vector<int> pointsToRotate;
	pointsToRotate = molecule.getRotamerSuccessors(rotamer);
	
	#ifdef DEBUG
	cout << "\n\nI want to rotate the following points of " << std::to_string(angle) << " degree:\t";
	for (int point : pointsToRotate)
		cout << std::to_string(point) << " ";
    #endif 
	
	matrix<float> rotationMatrix = createRotationMatrix(angle, rotamer.first, rotamer.second);
	
    #ifdef DEBUG
	cout << "\n\nThis is the rotation matrix:\n" << rotationMatrix << std::endl;
    #endif 
	
	//the for loop applies the rotation to all points in pointsToRotate
	for (int point : pointsToRotate)
		moleculeRotated.transform(rotationMatrix, point);
	
	#ifdef DEBUG
	cout << "\nAfter rotation:\n" << moleculeRotated.to_string();
	#endif 
	
	float score = calcolateScore(moleculeRotated, pocket);
	
	#ifdef DEBUG
	cout << "\nScore" << std::to_string(score);
	#endif 
	
	return score;
}


int main(int argc, char *argv[])
{
    string file_name = "NULL NAME";
    string n_string = "NULL NUMBER";
    int n = 0;
    
    po::options_description desc;
    
    desc.add_options()
        ("help, h", "Shows description of the options")
        ("file_name, f", po::value<string>(&file_name)->default_value("db.mol2"), "Set file name; if not setted <db.mol2> will be read.")
        ("number, n", po::value<string>(&n_string)->default_value("all"), "Set the number of the elements to be read; default value is <all>");
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if(vm.count("help"))
    {
        cout << desc;
        return 0;
    }
    
    if (file_name.compare("db.mol2") == 0 && n_string.compare("all") == 0)
        n = DB_DIMENSION;
    else if(n_string.compare("all") == 0)
        n = getDimension(file_name);
    else
        n = stoi(n_string);
    
    std::vector<Molecule> molecules = parseFile(file_name, n);
    
    float bestScore = 0;
    Molecule bestMolecule("");

    //create the pocket
    Pocket pocket(5, 5, 0.2);
    pocket.transformation();
    
    float numberOfProcessedAtoms=0;
    float throughput;
    
    //for calculating execution time 
    clock_t start,end;
    double executionTime;
    start=clock();
	
    for (Molecule molecule : molecules)
    {
    	numberOfProcessedAtoms+= molecule.getAtoms().size();
    	molecule.centre();

	//get all rotamers of the molecule
        std::vector<std::pair<Atom, Atom>> rotamers = molecule.getRotamers();
        int rotamerNumber = 1;
		#ifdef DEBUG
        std::cout << "\n\nList of rotamers";
		#endif 
        for (std::pair<Atom, Atom> rotamer : rotamers)
        {	
			#ifdef DEBUG
            cout << "\nRotamer number " << to_string(rotamerNumber) << " " << molecule.getAtomIndex(rotamer.first) << " " << molecule.getAtomIndex(rotamer.second);
			#endif 
            rotamerNumber++;
        }

        //cycle for each rotamer of the molecule
        for (std::pair<Atom, Atom> rotamer : rotamers)
        {
        	#ifdef DEBUG 
            std::cout << "\n\nI Consider the rotamer " << rotamer.first.to_string() << rotamer.second.to_string();
			#endif
            
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
			
			Molecule secondStepMolecule(bestLocalMolecule.getName());
			secondStepMolecule = copyMolecule(bestLocalMolecule);
            
            for (int angle = 0; angle<360; angle += 120)
            {
                Molecule moleculeRotated = Molecule(bestLocalMolecule.getName());
                float score = rotateMolecule(secondStepMolecule, moleculeRotated, angle, make_pair(rotamer.second, rotamer.first), pocket);
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
    //calculate the execution time from start to end and then the number of atoms processed per second
    end=clock();
    executionTime=((double)(end-start))/CLOCKS_PER_SEC;
    cout << "\n\nThe best score is: " << std::to_string(bestScore);
    cout << "\n\nBest molecule:\n " << bestMolecule.getName() << std::endl;
    cout << "\n\nExecution time : "<< executionTime;
    cout << "\n\nAtom processed : " << numberOfProcessedAtoms;
    throughput= numberOfProcessedAtoms/executionTime;
    cout << "\n\nThroughput : "<< throughput << "\n";
}
