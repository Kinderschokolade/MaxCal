#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <iostream>
#include <fstream>
#include <sstream> 

using namespace std;
int main(int argc, char *argv[])
{

	MPI::Init(argc, argv);
	const int np= MPI::COMM_WORLD.Get_size();
	const int root = 0;
	int rank=MPI::COMM_WORLD.Get_rank();
	printf("Running on rank %d \n", rank);
	string link = "test.dat";
	MPI_File Slist;
	char* linkc = new char[link.length() + 1];
	strcpy(linkc,link.c_str());
	MPI_File_delete(linkc,MPI_INFO_NULL);
	MPI_File_open(MPI_COMM_WORLD,linkc, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &Slist);
	
	double writev[5];
	double write;
	for (int i=0 ; i < (rank+1) ; i++)
	{
		for (int k =0 ; k <5 ; k++)
		{
		    writev[k] = (rank)+k;
		}
		int offset = (5 * rank + np * i * 5) * sizeof(double);
		MPI_File_write_shared(Slist, writev, 5, MPI_DOUBLE, MPI_STATUS_IGNORE);
		cout<< rank << "\t"<< i << endl;
	}
	
	MPI_File_close(&Slist);
	MPI::Finalize();

	int rank =0;
	int root =0;
	if(rank == root)
	{
		double in;
		FILE *f = fopen("test.dat", "r");
		std::cout << "Values : ";
		while(!feof(f) && !ferror(f))
		{
			for (int i = 0; i < 5; i++) {
				fread(&in, sizeof(double), 1, f);
				std::cout << in << "\t";
			}
			std::cout << std::endl;
		}
	}
}



