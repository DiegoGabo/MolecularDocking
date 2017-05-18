#ifndef structures_h
#define structures_h

#include <string>
#include <vector>
#include <list>

class Vertex
{
    private:
    
        float latitude;
        float longitude;
    
    public:
    
        void setLatitude(float latitude);

        void setLongitude(float longitude);

        float getLatitude();

        float getLongitude();

        std::string toString();
};


class Vertex3D
{
    
    private:
    
        float x;
        float y;
        float z;
    
    public:
    
        void setVariableX(float x);
        
        void setVariableY(float y);
        
        void setVariableZ(float z);
        
        std::string toString();
    
};

class Pocket
{
    private:
    
        float longmax,latmax;
        
        std::vector<Vertex> vertexMatrix;
        std::vector<Vertex3D> spherePoints;
        
        int i, j;
    
    public:
    
        Pocket(float latmax, float longmax);
        
        void transformation();
        
        std::string toString();
};


#endif
