#include "structures_molecule.hpp"
#include <CL/cl.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include "parser.hpp"
#include <math.h>
#define SIZE_POCKET 6 
const float PI = 3.14159265;

using namespace std;
namespace po = boost::program_options;
using namespace boost::algorithm;

void createPocket(Atom* pocket,float distance){
	
	struct Vertex2D{
		
		float latitude;
		float longitude;
		
	}vertex2D[SIZE_POCKET*SIZE_POCKET];
	
	float latAlfa,latBeta,lonAlfa,lonBeta,phi;
	
	//Creates a mesh of equidistant 2d points 
	
	for(int i=0;i<SIZE_POCKET;i++){
		for(int j=0;j<SIZE_POCKET;j++){
		
			vertex2D[i*SIZE_POCKET+j].latitude=i;
			vertex2D[i*SIZE_POCKET+j].longitude=j;
			
		}
	}
	
	   latAlfa= ((vertex2D[2*latmax]->latitude)*PI)/180;
	   latBeta= ((vertex2D[2*latmax+1]->latitude)*PI)/180;
	   lonAlfa= ((vertex2D[2*latmax]->longitude)*PI)/180;
	   lonBeta= ((vertex2D[2*latmax+1]->longitude)*PI)/180;
	   phi=fabs(lonAlfa-lonBeta); 
	   
	   float radius= distance/(acos(sin(latBeta)*sin(latAlfa)+cos(latBeta)*cos(latAlfa)*cos(phi)));
	   
	   //Transforms the bidimensional points of a mesh into coordinates of equidistant atoms in the sphere
	   
	   for(int i=0;i<SIZE_POCKET;i++){
			for(int j=0;j<SIZE_POCKET;j++){
				
				pocket->atom[i*SIZE_POCKET+j].x= radius*(sin(PI * vertex2D->latitude / SIZE_POCKET) *cos(2*PI * vertex.getLongitude() / SIZE_POCKET));
				pocket->atom[i*SIZE_POCKET+j].y= radius*(sin(PI * vertex2D->latitude / SIZE_POCKET) *sin(2*PI * vertex.getLongitude() / SIZE_POCKET));
				pocket->atom[i*SIZE_POCKET+j].z= radius*(cos(PI * vertex2D->latitude / SIZE_POCKET));
	  
			}
		}	  
}


int main( int argc, char** argv ) {

    const int N_ELEMENTS=1;
    unsigned int platform_id=0, device_id=0;

	Atom pocket[36];

    createPocket(pocket,0.2);
    
//	Atom a1, a2, a3, a4, a5, a6, a7, a8;
//	a1.x=0.0;  a1.y=1.0;  a1.z=0.0;
//	a2.x=0.0;  a2.y=-1.0; a2.z=0.0;
//	a3.x=1.0;  a3.y=0.0;  a3.z=0.0;
//	a4.x=2.0;  a4.y=0.0;  a4.z=0.0;
//	a5.x=3.0;  a5.y=1.0;  a5.z=0.0;
//	a6.x=4.0;  a6.y=-1.0; a6.z=0.0;
//	a7.x=4.0;  a7.y=0.0;  a7.z=0.0;
//	a8.x=-1.0; a8.y=-1.0;  a8.z=0.0;
//	
//	Molecule m1;
//	m1.numberOfAtoms=0;
//	addAtom(&m1, a1);
//	addAtom(&m1, a2);
//	addAtom(&m1, a3);
//	addAtom(&m1, a4);
//	addAtom(&m1, a5);
//	addAtom(&m1, a6);
//	addAtom(&m1, a7);
//	addAtom(&m1, a8);
//
//	setEdge(&m1, 0, 2);
//	setEdge(&m1, 1, 2);
//	setEdge(&m1, 1, 7);
//	setEdge(&m1, 2, 3);
//	setEdge(&m1, 3, 4);
//	setEdge(&m1, 3, 5);
//	setEdge(&m1, 4, 6);
//	setEdge(&m1, 5, 6);

    //aggiungere device GPU o CPU, togliere membro OR
    
    string file_name = "NULL";
    int n = 0;
    string device = "NULL";
    
    po::options_description desc;
    
    desc.add_options()
    ("help, h", "Shows description of the options")
    ("file_name, f", po::value<string>(&file_name)->default_value("ace_ligands.mol2"), "Set file name")
    ("number, n", po::value<int>(&n)->default_value(1), "Set the number of the elements to be read")
    ("device, d", po::value<string>(&device)->default_value("cpu"), "Set the type of device you want use. Available option: <gpu> or <cpu>");
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if(vm.count("help"))
        return 0;

    if(!vm.count("device"))
    {
        cout << "Device not setted\n";
        return 0;
    }
    
    
	//std::unique_ptr<Molecule[]> A(new Molecule[N_ELEMENTS]);
	Molecule* A = new Molecule[N_ELEMENTS];
//  A[0] = m1;
    
    A = parseFile(file_name, n);

	// Query for platforms
	std::vector<cl::Platform> platforms;
	cl::Platform::get(&platforms);
	// Get a list of devices on this platform
	std::vector<cl::Device> devices;
    
    // Select the platform.
    to_lower(device);
    if(device.compare("gpu") == 0)
        platforms[platform_id].getDevices(CL_DEVICE_TYPE_GPU, &devices);
    else if(device.compare("cpu") == 0)
        platforms[platform_id].getDevices(CL_DEVICE_TYPE_CPU, &devices);
    else
    {
        cout << "Unable to select the platform";
        return 0;
    }

	// Create a context
	cl::Context context(devices);

	// Create a command queue
	cl::CommandQueue queue = cl::CommandQueue( context, devices[device_id] );   // Select the device.
	// Create the memory buffers
	cl::Buffer bufferA=cl::Buffer(context, CL_MEM_READ_ONLY, N_ELEMENTS * sizeof(Molecule));

	// Copy the input data to the input buffers using the command queue.
	queue.enqueueWriteBuffer( bufferA, CL_FALSE, 0, N_ELEMENTS * sizeof(Molecule), A );

	// Read the program source
	std::ifstream sourceFile("kernel.cl");
	std::string sourceCode( std::istreambuf_iterator<char>(sourceFile), (std::istreambuf_iterator<char>()));
	cl::Program::Sources source(1, std::make_pair(sourceCode.c_str(), sourceCode.length()));

	// Make program from the source code
	cl::Program program=cl::Program(context, source);

	// Build the program for the devices
	program.build(devices);

	// Make kernel
	cl::Kernel vecadd_kernel(program, "doAllRotation");

	// Set the kernel arguments
	vecadd_kernel.setArg( 0, bufferA );

	// Execute the kernel
	cl::NDRange global( N_ELEMENTS );
	cl::NDRange local( 1 );
	queue.enqueueNDRangeKernel( vecadd_kernel, cl::NullRange, global, local );

	std::cout<< "Success!\n";
	return( EXIT_SUCCESS );
}
