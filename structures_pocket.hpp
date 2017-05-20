#ifndef structures_h
#define structures_h

#include <string>
#include "structures_molecule.hpp"
#include <vector>
#include <list>

class Vertex
{
    private:
    
        float latitude;
        float longitude;
    
    public:
        
        /*
        set the latitude
        */
        void setLatitude(float latitude);
        
        /*
        set the longitude
        */
        void setLongitude(float longitude);

        /*
        @return the latitude of the 2d point
        */
        float getLatitude();
  
        /*
        @return the longitude of the 2d point
        */
        float getLongitude();

        /*
        @return a string that describes the 2d vertex
        */
        std::string to_string();
};



class Pocket
{
    private:
    
        float longmax,latmax;
        
        std::vector<Vertex> vertexMatrix;
        std::vector<Atom> spherePoints;
        
        int i, j;
    
    public:

        /*
        Creates a pocket
        */

        Pocket(float latmax, float longmax);
        /*
         Transforms the bidimensional points of a mesh into coordinates of equidistant atoms in the sphere
         */
        void transformation();
       
        /*
         @return a string that describes the pocket and the transformation
         */
        std::string to_string();

	/*
	@return the set of equidistant atoms of the sphere
	*/
	std::vector<Atom> getAtoms();
};


#endif
