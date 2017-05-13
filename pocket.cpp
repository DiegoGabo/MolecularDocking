#include <iostream>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <math.h>

using namespace std;
using namespace boost::numeric::ublas;

class Vertex {

	float latitude;
	float longitude;

};
class Vertex3D {

	float x;
	float y;
	float z;
};
class Pocket
{
	float latmax,longmax;
	Vertex vertexMatrix[longmax * latmax];
	Vertex3D spherePoints[];
	float standardDistance;
	int i, j; 
public:
	Pocket()
	{
		for(i=0;i<latmax;i++)
			for (j = 0; j < longmax; j++)
			{
				vertexMatrix[i*j]->latitude = i*standardDistance;
				vertexMatrix[i*j]->longitude = j*standardDistance;
			}
		
	}
};

void Pocket::transformation()
{
	for (i = 0; i<latmax; i++)
		for (j = 0; j < longmax; j++)
		{
			spherePoints[i*j]->x = sin(Pi * vertexMatrix[i*j]->latitude / latmax) *cos(2Pi * vertexMatrix[i*j]->longitude / longmax);
			spherePoints[i*j]->y = sin(Pi * vertexMatrix[i*j]->latitude / latmax) *sin(2Pi * vertexMatrix[i*j]->longitude / longmax);
			spherePoints[i*j]->z = cos(Pi * vertexMatrix[i*j]->latitude / latmax);
		}

};