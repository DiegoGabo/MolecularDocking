//
//  main.cpp
//  parser
//
//  Created by Danilo Labanca on 14/05/17.
//  Copyright © 2017 Danilo Labanca. All rights reserved.
//
// how to understand mol2 contents -> http://thegrantlab.org/bio3d/html/read.mol2.html

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>

using namespace std;

const string MOLECULE_BEGIN ("@<TRIPOS>MOLECULE");
const string MOLECULE_ATOMS ("@<TRIPOS>ATOM");
const string MOLECULE_BONDS ("@<TRIPOS>BOND");
const string MOLECULE_NULL ("EMPTY MOLECULE");


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
        return "x = " + to_string(x) + " y = " + to_string(y) + " z = " + to_string(z);
    }
    
};

class Molecule
{
    vector<Atom> atoms;
    vector<list<int>> links;
    string name;
    
    
public:
    
    Molecule(string name)
    {
        this->name = name;
        atoms.clear();
        links.clear();
    }
    
    void setName(string name)
    {
        this->name = name;
    }
    
    string getName()
    {
        return name;
    }
    
    void addAtom(Atom atom)
    {
        atoms.push_back(atom);
        list<int> link;
        links.push_back(link);
    }
    
    void setEdge(int src, int dest)
    {
        links.at(src).push_back(dest);
    }
    
    
    string toString()
    {
        
        cout << "Molecule name: " << name << '\n';
        
        //creates a string with the list of atoms
        string molecule = "The molecule has the following atoms";
        for (Atom atom : atoms)
        {
            molecule += "\n";
            molecule += atom.toString();
        }
        
        //creates a string with the structure of the graph
        molecule += "\n\nGraph structure:";
        int atom = 0;
        for (list<int> link : links)
        {
            molecule += "\nAtom:  ";
            molecule += std::to_string(atom);
            molecule += " linked to: ";
            for (list<int>::iterator it = link.begin(); it != link.end(); ++it)
            {
                molecule += to_string(*it);
                molecule += ", ";
            }
            atom++;
        }
        return molecule;
    }
};

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


int main()
{
    
    fstream ligands("ace_ligands.mol2");
    
    if( ligands.is_open())
    {
        cout << "File ok\n\n\n";
        
        string line;
        vector<string> text_elements;
        
        Molecule* temp_molecule = new Molecule(MOLECULE_NULL);
        vector<Molecule> molecules;
        
        bool name_flag = false;
        bool atoms_flag = false;
        bool bonds_flag = false;
        
        while (getline (ligands, line))
        {
            text_elements = splitLine(line);
            
            if (text_elements.empty())
                continue;
            
            if (text_elements[0].compare(MOLECULE_BEGIN) == 0)
            {
                //if the first element we meet is @<TRIPOS> MOLECULE this means that we are at the begin of the description of a molecule
                //in this case we set the relative flag to true and, if it's possible, we insert the new molecule in the vector and instanziates a new empty molecule.
                
                name_flag = true;
                atoms_flag = false;
                bonds_flag = false;
                
                if (temp_molecule->getName().compare(MOLECULE_NULL) != 0)
                {
                    molecules.push_back(*temp_molecule);
                    temp_molecule = new Molecule(MOLECULE_NULL);
                }
            }
            else if (text_elements[0].compare(MOLECULE_ATOMS) == 0)
            {
                cout << "These are its atoms\n";    //just for testing
                
                //if the first element we meet is @<TRIPOS> ATOM means that the next lines describes the atoms of the molecule
                name_flag = false;
                atoms_flag = true;
                bonds_flag = false;
            }
            else if (text_elements[0].compare(MOLECULE_BONDS) == 0)
            {
                cout << "These are its bonds\n";    //just for testing
                
                //if the first element we meet is @<TRIPOS> BONDS means that the next lines describes the bonds between the atoms of the molecule
                name_flag = false;
                atoms_flag = false;
                bonds_flag = true;
            }
            else    //we follow this branch if we are exctracting information about the name, the atoms or the bonds of the molecule
            {
                if (!name_flag && !atoms_flag && !bonds_flag)
                {
                    continue;
                }
                else if (name_flag && !atoms_flag && !bonds_flag)
                {
                    temp_molecule->setName(text_elements[0]);
                    name_flag = false;
                }
                else if (!name_flag && atoms_flag && !bonds_flag)
                {
                    float x = stof(text_elements[2]);
                    float y = stof(text_elements[3]);
                    float z = stof(text_elements[4]);
                    
                    temp_molecule->addAtom(Atom (x, y, z));
                }
                else if (!name_flag && !atoms_flag && bonds_flag)
                {
                    int source = stoi(text_elements[1]);
                    int destination = stoi(text_elements[2]);
                    
                    temp_molecule->setEdge(source, destination);
                }
            }
            
            cout << '\n';
            
//            understand how the lines int the mol2 file are made up
//
//            for (int i = 0; i<elements.size(); i++)
//                cout << i << ": " << elements.at(i) << '\n';
//            
//            cout << '\n';
        }
        
        ligands.close();
        
        cout << "\n\nThese are the molecules:\n\n";
        
        for (int i = 0; i < molecules.size(); i++)
        {
            cout << molecules[i].toString() << "\n\n";
            cout << "-----------------------------------------------------------------\n";
        }
        
            
    }
    else
        cout << "Unable to open the file\n";
    
    return 0;
}
