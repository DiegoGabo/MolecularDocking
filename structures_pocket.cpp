#include "structures_pocket.hpp"
#include "structures_molecule.hpp"
#include <iostream>
#include <vector>
#include <list>
#include <math.h>

using namespace std;

const float PI = 3.14159265;


//Implementations of Vertex class functions

/*
set the latitude
*/
void Vertex::setLatitude(float latitude)
{
    this->latitude = latitude;
}
/*
	 set the longitude
 */
void Vertex::setLongitude(float longitude)
{
    this->longitude = longitude;
}
/*
 @return the latitude of the 2d coordinate
 */
float Vertex::getLatitude()
{
    return this->latitude;
}
/*
 @return the longitude of the 2d coordinate
 */
float Vertex::getLongitude()
{
    return this->longitude;
}

/*
 @return a string that describes the 2d vertex
 */
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

/*
Creates a mesh of equidistant (latmax*longmax) 2d points 
*/
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

/*
 @return the set of equidistant atoms of the sphere
 */
vector<Atom> Pocket::getAtoms()
{
    return spherePoints;
}

/*
 Transforms the bidimensional points of a mesh into coordinates of equidistant atoms in the sphere
 */
void Pocket::transformation()
{
    
    for(Vertex vertex: vertexMatrix)
    {
	int x, y, z;
        x = sin(PI * vertex.getLatitude() / latmax) *cos(2*PI * vertex.getLongitude() / longmax);
        y = sin(PI * vertex.getLatitude() / latmax) *sin(2*PI * vertex.getLongitude() / longmax);
        z = cos(PI * vertex.getLatitude() / latmax);
	Atom vertex3d(x, y, z);	
        spherePoints.push_back(vertex3d);
    }

    vector<Atom>::iterator begin = spherePoints.begin();
    vector<Atom>::iterator end = begin + latmax - 1;
    spherePoints.erase(begin, end);
    
    Atom vertex(0,0,-1);
    spherePoints.push_back(vertex);
}

/*
 @return a string that describes the pocket and the transformation
 */
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
