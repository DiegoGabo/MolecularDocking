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



int main()
{
    fstream ligand("ligand_dataset.txt");
    string line;
    vector<string> elements;
    vector<Atom> atoms;
    Molecule molecules;
    
    if( ligand.is_open())
    {
        cout << "File ok\n";
        
        while (getline (ligand, line))
        {
            int atoms_index = 0;
            elements = split_line(line);
            
            if(elements[0].compare(="ATOM"))
            {
                Atom new_atom (elements[6], elements[7], elements[8]);
                atoms.push_back(new_atom);
            }
            else
            {
                
                
            }
            
//            for (int i = 0; i<elements.size(); i++)
//                cout << i << ": " << elements.at(i) << '\n';
//            
//            cout << '\n';
        }
        
        ligand.close();
    }
    else cout << "Unable to open the file\n";
    
    return 0;
}
