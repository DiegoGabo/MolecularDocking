// how to understand mol2 contents -> http://thegrantlab.org/bio3d/html/read.mol2.html

#include "parser.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include "structures_molecule.hpp"

using namespace std;

const string MOLECULE_BEGIN ("@<TRIPOS>MOLECULE");
const string MOLECULE_ATOMS ("@<TRIPOS>ATOM");
const string MOLECULE_BONDS ("@<TRIPOS>BOND");
const string MOLECULE_NULL ("EMPTY MOLECULE");
const

/*
 This method takes a string and returns a vector containing each word that composes the string
 @param the string to analyse
 @return a vector containing string
 */
vector<string> splitLine(string& line)
{
    vector<string> elements;
    string temp = "";
    
    for (int i = 0; i < line.size(); i++)
    {
        if (line.at(i) != ' ')
            temp.append(line.substr(i,1));
        else
        {
            if (temp.compare("") != 0)
            {
                elements.push_back(temp);
                temp.clear();
            }
        }
    }
    
    if (temp.compare("") != 0)
        elements.push_back(temp);
    
    return elements;
    
}


/*
This method takes a file as parameter and returns a vectors containing the molecules described in the file
 @param a file to be parsed
 @return a vector containing molecules
 */
Molecule* parseFile(string name, int l)
{
    int limit = 0;
    int dimension = 0;
    
    //if limit is zero the parser will read all the ligands
    if (l < 0)
    {
        limit = l;
        dimension = -l;
    }
    else if( l == 0)
    {
        limit = 1;
        dimension = 50;
    }
    else
    {
        limit = -l;
        dimension = l;
    }
    
    Molecule* molecules = new Molecule[dimension];
    
    fstream ligands(name);
    
    if( ligands.is_open())
    {
        string line;
        vector<string> text_elements;
        
        bool name_flag = false;
        bool atoms_flag = false;
        bool bonds_flag = false;
        
        int mol_index = -1;
        
        while (getline(ligands, line) and limit)
        {
            text_elements = splitLine(line);
            
            if (text_elements.empty())
                continue;
            
            if (text_elements[0].compare(MOLECULE_BEGIN) == 0)
            {
                //if the first element we meet is "@<TRIPOS> MOLECULE" this means that we are at the begin of the description of a molecule
                //in this case we set the relative flag to true and, if it's possible, we insert the new molecule in the vector and instanziates a new empty molecule.
                
                name_flag = true;
                atoms_flag = false;
                bonds_flag = false;

                mol_index++;

            }
            else if (text_elements[0].compare(MOLECULE_ATOMS) == 0)
            {
                //if the first element we meet is "@<TRIPOS> ATOM" means that the next lines describes the atoms of the molecule
                name_flag = false;
                atoms_flag = true;
                bonds_flag = false;
            }
            else if (text_elements[0].compare(MOLECULE_BONDS) == 0)
            {
                //if the first element we meet is "@<TRIPOS> BONDS" means that the next lines describes the bonds between the atoms of the molecule
                name_flag = false;
                atoms_flag = false;
                bonds_flag = true;
            }
            else    //we follow this branch if we are exctracting information about the name, the atoms or the bonds of the molecule
            {
                if (!name_flag && !atoms_flag && !bonds_flag)
                    continue;   //we are skipping the part between "@<TRIPOS> MOLECULE" and "@<TRIPOS> ATOM"
                else if (name_flag && !atoms_flag && !bonds_flag)
                {
                    molecules[mol_index].name = text_elements[0];
                    name_flag = false;
                }
                else if (!name_flag && atoms_flag && !bonds_flag)
                {
                    float x = stof(text_elements[2]);
                    float y = stof(text_elements[3]);
                    float z = stof(text_elements[4]);
                    
                    addAtom(&molecules[mol_index], x, y, z);
                }
                else if (!name_flag && !atoms_flag && bonds_flag)
                {
                    int source = stoi(text_elements[1]) -1 ;
                    int destination = stoi(text_elements[2]) - 1;
                    
                    setEdge(&molecules[mol_index], source, destination);
                }
            }
        }
        
        ligands.close();
    }
    else
        cout << "Unable to open the file\n";
    
    return &molecules;
}
