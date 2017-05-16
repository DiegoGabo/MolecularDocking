#include <iostream>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <math.h>
#include <string>
#define PI 3.14159265

class Vertex {

	float latitude;
	float longitude;
	public:
		void setLatitude(float latitude)
		{
			this->latitude = latitude;
		}
		void setLongitude(float longitude)
		{
			this->longitude = longitude;
		}
		float getLatitude()
		{
			return this->latitude;
		}
		float getLongitude()
		{
			return this->longitude;
		}
		std::string toString()
		{
			std::string vertex;
			vertex="Vertex with latitude ";
			vertex += std::to_string(latitude);
			vertex += " and longitude ";
			vertex += std::to_string(longitude);
			return vertex;
		}
		

};
class Vertex3D {

	float x;
	float y;
	float z;
	
		public:
		void setVariableX(float x)
		{
			this->x = x;
		}
		void setVariableY(float y)
		{
			this->y = y;
		}
		void setVariableZ(float z)
		{
			this->z = z;
		}
		std::string toString()
		{
			std::string vertex="Vertex with x = ";
			vertex += std::to_string(x);
			vertex += " and y = ";
			vertex += std::to_string(y);
			vertex += " and z = ";
			vertex += std::to_string(z);
			return vertex;
		}
	
};
class Pocket
{
	float longmax,latmax;
	std::vector<Vertex> vertexMatrix;
	std::vector<Vertex3D> spherePoints;
	int i, j;
	public:
		Pocket(float latmax, float longmax)
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

	public:
	
		void transformation()
		{
			
			for(Vertex vertex: vertexMatrix)
			{
				Vertex3D vertex3d;
				vertex3d.setVariableX(sin(PI * vertex.getLatitude() / latmax) *cos(2*PI * vertex.getLongitude() / longmax));
				vertex3d.setVariableY(sin(PI * vertex.getLatitude() / latmax) *sin(2*PI * vertex.getLongitude() / longmax));
				vertex3d.setVariableZ(cos(PI * vertex.getLatitude() / latmax));
				spherePoints.push_back(vertex3d);
			}

		}
		
		std::string toString()
		{
			std::string pocket = "The vector2d:";
			for (Vertex vertex: vertexMatrix)
			{
				pocket += "\n";
				pocket += vertex.toString();
			}
			
			pocket += "\nThe sphere:";
			for (Vertex3D vertex3D: spherePoints)
			{
				pocket += "\n";
				pocket += vertex3D.toString();
			}
			return pocket;
		}	
};

int main()
{
	Pocket pocket(10,10);
	pocket.transformation();
	std::cout << "ciao";
	std::cout << pocket.toString();
}