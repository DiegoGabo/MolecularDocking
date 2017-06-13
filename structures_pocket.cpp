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
string Vertex::to_string()
{
    string vertex;
    vertex="Vertex with latitude ";
    vertex += std::to_string(latitude);
    vertex += " and longitude ";
    vertex += std::to_string(longitude);
    
    return vertex;
}


//Implementations of Pocket class functions

/*
Creates a mesh of equidistant (latmax*longmax) 2d points 
*/
Pocket::Pocket(float latmax, float longmax,float distance)
{	
	this->distance=distance;
    this->latmax=latmax;
    this->longmax=longmax;
    
    for(int i=0;i<latmax;i++)
        for (int j = 0; j < longmax; j++)
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
vector<Atom> Pocket::getAtoms() const
{
    return spherePoints;
}

/*
 Transforms the bidimensional points of a mesh into coordinates of equidistant atoms in the sphere
 */
void Pocket::transformation()
{
    float latAlfa,latBeta,lonAlfa,lonBeta,phi;
    
    latAlfa= (vertexMatrix.at(2*latmax).getLatitude()*PI)/180;
    latBeta= (vertexMatrix.at(2*latmax+1).getLatitude()*PI)/180;
	lonAlfa= (vertexMatrix.at(2*latmax).getLongitude()*PI)/180;
	lonBeta= (vertexMatrix.at(2*latmax+1).getLongitude()*PI)/180;
	phi=fabs(lonAlfa-lonBeta);
    float radius= distance/(acos(sin(latBeta)*sin(latAlfa)+cos(latBeta)*cos(latAlfa)*cos(phi)));
    
    for(Vertex vertex: vertexMatrix)
    {
    	float x, y, z;
        x = sin(PI * vertex.getLatitude() / latmax) *cos(2*PI * vertex.getLongitude() / longmax);
        y = sin(PI * vertex.getLatitude() / latmax) *sin(2*PI * vertex.getLongitude() / longmax);
        z = cos(PI * vertex.getLatitude() / latmax);
        Atom vertex3d(radius * x, radius * y, radius * z);	
        spherePoints.push_back(vertex3d);
    }

}

/*
 @return a string that describes the pocket and the transformation
 */
string Pocket::to_string()
{
    string pocket = "The vector2d:";
    for (Vertex vertex: vertexMatrix)
    {
        pocket += "\n";
        pocket += vertex.to_string();
    }
    
    pocket += "\nThe sphere:";
    for (Atom vertex3D: spherePoints)
    {
        pocket += "\n";
        pocket += vertex3D.to_string();
    }
    return pocket;
}
