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
	int N = 2;
	int sizes[2]    = {N*procN, N*procN}; // size after gathering
	int subsizes[2] = {N, N}; // initial size
    	int starts[2]   = {0,0}; // where this one starts
	double* ar_gath = nullptr;
	MPI_Datatype type,subarrtype;
	MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &type);
    	MPI_Type_create_resized(type, 0, N * sizeof(double), &subarrtype);
    	MPI_Type_commit(&subarrtype);

	int sendcounts[procN*procN];
	int displacement[procN*procN];
	double* ar = new double[N  * N];

	if (rank == root)
	{
		ar_gath = new double[N*N*np];
		for (int i=0 ; i < procN*procN ; i++)
		{
			sendcounts[i] =1;
		}
       		int disp = 0;
		for (int i = 0; i < procN; i++) {
            		for (int j = 0; j < procN; j++) {
                		displacement[i * procN + j] = disp;
                		disp += 1;
            		}
            		disp += (N  - 1) * procN;
        	}
	        for (int i=0 ; i < N*N*np ; i++)
        	{
                	ar_gath[i] = 0.;
        	}
	}

	for (int i=0 ; i < N ; i++)
	{
		for (int j=0 ; j < N ; j++)
		{
			ar[i*N+j] = rank;
		}
	}


	for (int p = 0; p < np; p++) {
    		if (rank == p) {
    			printf("ar on rank %d:\n", rank);
            		for (int i = 0; i < N ; i++) {
                		for (int j = 0; j < N ; j++) {
                    			printf("%f ", ar[j * N  + i]);
                		}
                		printf("\n");
            		}
        	}
    		MPI_Barrier(MPI_COMM_WORLD);
    	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gatherv(ar, N*N, MPI_DOUBLE , ar_gath, sendcounts,displacement, subarrtype, root, MPI_COMM_WORLD);

	vector < vector <double> > ar_sum;
	if (rank ==0)
	{
		ar_sum.resize(N);
		for (int i=0 ; i<N  ; i++)
		{
			ar_sum[i].resize(N,0);
		}

		for (int p=0 ; p< np ; p++)
		{
			for (int i =0 ; i< N ; i++)
			{
				for (int j=0 ; j<N ; j++)
				{
					ar_sum[i][j] += ar_gath[displacement[p]+ i*N + j*N*N*procN];
				} // this doesnt work
			}
		}
		

		for (int i=0 ; i < N*procN *N*procN ; i++)
		{
			//for (int j=0 ; j < N*procN ; j++)
			//{
				cout<< ar_gath[i] <<"\t";
			//}
		}
		cout<<endl;	

		for (int i=0 ; i < N ; i++)
		{
			for (int j=0 ; j < N ; j++)
			{
				cout<< ar_sum[i][j] <<"\t";
			}
			cout<<endl;
		}

	}

	delete ar;
	delete ar_gath;
	MPI_Finalize();

}


