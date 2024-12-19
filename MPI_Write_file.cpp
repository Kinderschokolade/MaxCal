#include "mpi.h"
#include <string.h> 
#include <iostream>
#include <fstream>
#include <sstream> 

using namespace std;


int main(int argc, char *argv[])
{

	MPI::Init(argc, argv);
	const int np= MPI::COMM_WORLD.Get_size (  );
	const int root = 0;
	int rank=MPI::COMM_WORLD.Get_rank ( );

	string link = "test.dat";
	MPI_File Slist;
	char* linkc = new char[link.length() + 1];
	strcpy(linkc,link.c_str());
	MPI_File_open(MPI_COMM_WORLD,linkc, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &Slist);


	double* writev = new double [5];
	double write;
	for (int i=0 ; i < 5 ; i++)
	{
		for (int k =0 ; k <5 ; k++)
		{
			writev[k] = (rank)+k;
//			cout<< writev[k] << "\t";
		}
//		cout<< endl;
//		write = i *(rank+1);

//		MPI_File_write_at(Slist,i * rank * 5 *sizeof(double), writev,5,MPI_DOUBLE, MPI_STATUS_IGNORE);
		MPI_File_write_shared(Slist, &write ,1,MPI_DOUBLE, MPI_STATUS_IGNORE);
		

	}

	MPI_File_close(&Slist);
	MPI::Finalize();
}


