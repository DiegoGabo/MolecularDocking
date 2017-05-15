//
//  main.cpp
//  parser
//
//  Created by Danilo Labanca on 14/05/17.
//  Copyright © 2017 Danilo Labanca. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>

using namespace std;

vector<string> split_line(string& line)
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
    
    
    return elements;
}

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
        return "x = " + std::to_string(x) + " y = " + std::to_string(y) + " z = " + std::to_string(z);
    }
    
};

class Molecule
{
    std::vector<Atom> atoms;
    std::vector<std::list<int>> links;
    
    
public:
    
    void addAtom(Atom atom)
    {
        atoms.push_back(atom);
        std::list<int> link;
        links.push_back(link);
    }
    
    void setEdge(int src, int dest)
    {
        links.at(src).push_back(dest);
    }
    
    
    string toString()
    {
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
        for (std::list<int> link : links)
        {
            molecule += "\nAtom:  ";
            molecule += std::to_string(atom);
            molecule += " linked to: ";
            for (std::list<int>::iterator it = link.begin(); it != link.end(); ++it)
            {
                molecule += std::to_string(*it);
                molecule += ", ";
            }
            atom++;
        }
        return molecule;
    }
};


int main()
{
    fstream ligand("ligand_dataset.txt");
    
    
    if( ligand.is_open())
    {
        cout << "File ok\n";
        
        string line;
        vector<string> elements;
        Molecule molecule;
        int atoms_index = 0;
        
        int count_atoms = 0;    //just for testing
        
        while (getline (ligand, line))
        {
            elements = split_line(line);
            
            if(elements[0].compare("ATOM") == 0)
            {
                cout << "I find a new atom #" << ++count_atoms << '\n';
                
                float x = stof(elements[6]);
                float y = stof(elements[7]);
                float z = stof(elements[8]);
                
                Atom new_atom (x, y, z);
                
                cout << "It's composition " << new_atom.toString() << '\n';
                
                molecule.addAtom(new_atom);
            }
            else
            {
                cout << "I'm connecting some atoms\n";
                
                for (int i = 1; i < elements.size(); i++)
                {
                    cout << "Edge: " << atoms_index << "->" << elements[i] << '\n';
                    molecule.setEdge(atoms_index, stoi(elements[i]));
                }
                
                atoms_index++;
            }
            
            cout << '\n';
            
//            for (int i = 0; i<elements.size(); i++)
//                cout << i << ": " << elements.at(i) << '\n';
//            
//            cout << '\n';
        }
        
        ligand.close();
        
        cout << "The resultant molecule:\n";
        cout << molecule.toString();
    }
    else cout << "Unable to open the file\n";
    
    return 0;
}
