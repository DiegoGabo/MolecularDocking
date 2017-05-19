#include "structures_pocket.hpp"
#include "structures_molecule.hpp"
#include <iostream>
#include <vector>
#include <list>
#include <math.h>

using namespace std;

const float PI = 3.14159265;


//Implementations of Vertex class functions

void Vertex::setLatitude(float latitude)
{
    this->latitude = latitude;
}

void Vertex::setLongitude(float longitude)
{
    this->longitude = longitude;
}

float Vertex::getLatitude()
{
    return this->latitude;
}

float Vertex::getLongitude()
{
    return this->longitude;
}

std::vector<Atom> getAtoms()
{
	return spherePoints;
}

string Vertex::toString()
{
    string vertex;
    vertex="Vertex with latitude ";
    vertex += to_string(latitude);
    vertex += " and longitude ";
    vertex += to_string(longitude);
    
    return vertex;
}


//Implementations of Pocket class functions

Pocket::Pocket(float latmax, float longmax)
{
    this->latmax=latmax;
    this->longmax=longmax;
    
    for(i=0;i<latmax;i++)
        for (j = 0; j < longmax; j++)
        {
            Vertex vertex;
            vertex.setLatitude(i);
            vertex.setLongitude(j);
            vertexMatrix.push_back(vertex);
        }
}

void Pocket::transformation()
{
    
    for(Vertex vertex: vertexMatrix)
    {
        Atom vertex3d(0,0,0);
        vertex3d.setVariableX(sin(PI * vertex.getLatitude() / latmax) *cos(2*PI * vertex.getLongitude() / longmax));
        vertex3d.setVariableY(sin(PI * vertex.getLatitude() / latmax) *sin(2*PI * vertex.getLongitude() / longmax));
        vertex3d.setVariableZ(cos(PI * vertex.getLatitude() / latmax));
        spherePoints.push_back(vertex3d);
    }

    vector<Atom>::iterator begin = spherePoints.begin();
    vector<Atom>::iterator end = begin + latmax - 1;
    spherePoints.erase(begin, end);
    
    Atom vertex;
    vertex.setVariableX(0);
    vertex.setVariableY(0);
    vertex.setVariableZ(-1);
    spherePoints.push_back(vertex);
}

string Pocket::toString()
{
    string pocket = "The vector2d:";
    for (Vertex vertex: vertexMatrix)
    {
        pocket += "\n";
        pocket += vertex.toString();
    }
    
    pocket += "\nThe sphere:";
    for (Atom vertex3D: spherePoints)
    {
        pocket += "\n";
        pocket += vertex3D.toString();
    }
    return pocket;
}
