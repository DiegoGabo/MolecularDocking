#include "structures_molecule.hpp"
#include <CL/cl.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
void addAtom(Molecule* molecule, Atom atom)
{
	int index = molecule->numberOfAtoms;
	molecule->atoms[index]=atom;
	molecule->links[index].numberOfLinks=0;
	molecule->numberOfAtoms += 1;
}


void setEdge(Molecule * molecule, int src, int dest)
{
	int index = molecule->links[src].numberOfLinks;
	molecule->links[src].atoms[index]=dest;
	molecule->links[src].numberOfLinks+=1;
	
	index = molecule->links[dest].numberOfLinks;
	molecule->links[dest].atoms[index]=src;
	molecule->links[dest].numberOfLinks+=1;
}

int main( int argc, char** argv ) {

    const int N_ELEMENTS=1;
    unsigned int platform_id=0, device_id=0;

	Atom a1, a2, a3, a4, a5, a6, a7, a8;
	a1.x=0.0;  a1.y=1.0;  a1.z=0.0;
	a2.x=0.0;  a2.y=-1.0; a2.z=0.0;
	a3.x=1.0;  a3.y=0.0;  a3.z=0.0;
	a4.x=2.0;  a4.y=0.0;  a4.z=0.0;
	a5.x=3.0;  a5.y=1.0;  a5.z=0.0;
	a6.x=4.0;  a6.y=-1.0; a6.z=0.0;
	a7.x=4.0;  a7.y=0.0;  a7.z=0.0;
	a8.x=-1.0; a8.y=-1.0;  a8.z=0.0;
	
	Molecule m1;
	m1.numberOfAtoms=0;
	addAtom(&m1, a1);
	addAtom(&m1, a2);
	addAtom(&m1, a3);
	addAtom(&m1, a4);
	addAtom(&m1, a5);
	addAtom(&m1, a6);
	addAtom(&m1, a7);
	addAtom(&m1, a8);

	setEdge(&m1, 0, 2);
	setEdge(&m1, 1, 2);
	setEdge(&m1, 1, 7);
	setEdge(&m1, 2, 3);
	setEdge(&m1, 3, 4);
	setEdge(&m1, 3, 5);
	setEdge(&m1, 4, 6);
	setEdge(&m1, 5, 6);

	//std::unique_ptr<Molecule[]> A(new Molecule[N_ELEMENTS]); 
	Molecule* A = new Molecule[N_ELEMENTS];
	A[0] = m1;

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
