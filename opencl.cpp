#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <stdlib.h>

#define __CL_ENABLE_EXCEPTIONS
#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.cpp>
#else
#include <CL/cl.hpp>
#endif

#define MAX_SIZE 3

struct Atom
{
	float x, y, z;
};

struct Molecule
{
	Atom atoms[MAX_SIZE];
    int numberOfAtoms;
};


int main( int argc, char** argv ) {

    const int N_ELEMENTS=1;
    unsigned int platform_id=0, device_id=0;

	Atom a1, a2, a3;
	a1.x=1.0;
	a1.y=1.0;
	a1.z=1.0;
	a2.x=2.0;
	a2.y=2.0;
	a2.z=2.0;
	a3.x=3.0;
	a3.y=3.0;
	a3.z=3.0;

	Molecule m;
	m.atoms[0]=a1;
	m.atoms[1]=a2;
	m.atoms[2]=a3;
	m.numberOfAtoms = 3;

	try 
	{
		std::unique_ptr<Atom[]> A(new Atom[N_ELEMENTS]); // Or you can use simple dynamic arrays like: int* A = new int[N_ELEMENTS];
		A[0] = a1;

		// Query for platforms
		std::vector<cl::Platform> platforms;
		cl::Platform::get(&platforms);
		// Get a list of devices on this platform
		std::vector<cl::Device> devices;
		platforms[platform_id].getDevices(CL_DEVICE_TYPE_GPU|CL_DEVICE_TYPE_CPU, &devices); // Select the platform.

		// Create a context
		cl::Context context(devices);

		// Create a command queue
		cl::CommandQueue queue = cl::CommandQueue( context, devices[device_id] );   // Select the device.
		// Create the memory buffers
		cl::Buffer bufferA=cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(Atom));

		// Copy the input data to the input buffers using the command queue.
		queue.enqueueWriteBuffer( bufferA, CL_FALSE, 0, N_ELEMENTS * sizeof(Atom), A.get() );

		// Read the program source
		std::ifstream sourceFile("vector_add_kernel.cl");
		std::string sourceCode( std::istreambuf_iterator<char>(sourceFile), (std::istreambuf_iterator<char>()));
		cl::Program::Sources source(1, std::make_pair(sourceCode.c_str(), sourceCode.length()));

		// Make program from the source code
		cl::Program program=cl::Program(context, source);

		// Build the program for the devices
		program.build(devices);

		// Make kernel
		cl::Kernel vecadd_kernel(program, "vecadd");

		// Set the kernel arguments
		vecadd_kernel.setArg( 0, bufferA );

		// Execute the kernel
		cl::NDRange global( N_ELEMENTS );
		cl::NDRange local( 1 );
		queue.enqueueNDRangeKernel( vecadd_kernel, cl::NullRange, global, local );

		std::cout<< "Success!\n";
		return( EXIT_SUCCESS );
	}
	catch(cl::Error err) {
		    std::cout << "Error: " << err.what() << "(" << err.err() << ")" << std::endl;
		    return( EXIT_FAILURE );
		
		}
	
}
