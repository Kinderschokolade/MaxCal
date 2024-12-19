#include <iostream>
#include <vector> 
#include "mpi.h" 
#include <math.h> 
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{

	cout<< "hello"<<endl;
   	MPI_Init(&argc, &argv);
        const int np= MPI::COMM_WORLD.Get_size (  );
	const int root=0;
        int rank=MPI::COMM_WORLD.Get_rank ( );
	cout<< "s"<< np << endl;	
	cout<< "p" << rank <<endl;

	MPI_Finalize();

}


