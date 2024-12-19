#include <iostream>
#include <vector> 
#include "mpi.h" 
#include <math.h> 

using namespace std;

int main(int argc, char *argv[])
{

   	MPI_Init(&argc, &argv);
        const int np= MPI::COMM_WORLD.Get_size (  );
	int procN = sqrt(np); // restriction on number of processors: 1,4,9,16,25,..
	const int root=0;
        int rank=MPI::COMM_WORLD.Get_rank ( );
	
	int seed = 1; //read
	seed = seed *(rank+1);
	int N = 4; // length of vector

	vector < vector <double> > ar(N,vector<double> (N,0));

	for (int i=0 ; i < N ; i++)
	{
		for (int j=0 ; j< N ; j++)
		{
			ar[i][j] = rank;
		}
	}


	vector < vector <double> > ar_sum;
	if (rank ==0)
	{
		ar_sum.resize(N);
		for (int i=0 ; i<N  ; i++)
		{
			ar_sum[i].resize(N,0);
		}
	}

	vector <double> sendvec(N,0);
	vector <double> recvec(N,0);
	for (int k=0 ; k<N; k++)
	{
		sendvec = ar[k];
 		MPI::COMM_WORLD.Reduce(&sendvec[0],&recvec[0], sendvec.size(), MPI_DOUBLE,MPI_SUM, root);
		if (rank ==0)
		{
			ar_sum[k] = recvec;
		}
	}

	if (rank ==root)
	{
		for (int i=0 ; i < N ; i++)
		{
			for (int j=0 ; j < N ; j++)
			{
				cout<< ar_sum[i][j] <<"\t";
			}
			cout<<endl;
		}

	}

	MPI_Finalize();

}


