#include <iostream>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <math.h>
#define PI 3.14159265

class Vertex {

	float latitude;
	float longitude;
	public:
		void setLatitude(float l)
		{
			latitude = l;
		}
		void setLongitude(float l)
		{
			latitude = l;
		}

};
class Vertex3D {

	float x;
	float y;
	float z;
};
class Pocket
{
	float latmax=10,longmax=10;
	Vertex vertexMatrix[longmax * latmax];
	Vertex3D spherePoints[];
	int i, j;
public:
	Pocket(float standardDistance)
	{
		for(i=0;i<latmax;i++)
			for (j = 0; j < longmax; j++)
			{
				vertexMatrix[latmax*i+j].setLatitude(i*standardDistance);
				vertexMatrix[latmax*i+j].setLongitude(j*standardDistance);
			}

	}

	public:
	
		void transformation()
		{
			for (i = 0; i<latmax; i++)
				for (j = 0; j < longmax; j++)
				{
					spherePoints[i*j]->x = sin(PI * vertexMatrix[i*j]->latitude / latmax) *cos(2Pi * vertexMatrix[i*j]->longitude / longmax);
					spherePoints[i*j]->y = sin(PI * vertexMatrix[i*j]->latitude / latmax) *sin(2Pi * vertexMatrix[i*j]->longitude / longmax);
					spherePoints[i*j]->z = cos(PI * vertexMatrix[i*j]->latitude / latmax);
				}

		}
		
};

int main()
{
	Pocket pocket(1);
}