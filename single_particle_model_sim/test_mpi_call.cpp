#include "mpi.h"
#include <iostream> 

#define N 4 // grid size
#define procN 2  // size of process grid

using namespace std;

int main(int argc, char **argv) {
    double* gA = nullptr; // pointer to array
    int rank, size;       // rank of current process and no. of processes

    // mpi initialization
    MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // force to use correct number of processes
    if (size != procN * procN) {
		if (rank == 0) fprintf(stderr,"%s: Only works with np = %d.\n", argv[0], procN *  procN);
        MPI_Abort(MPI_COMM_WORLD,1);
    }

    // allocate and print global A at master process
    if (rank == 0) {
        gA = new double[N * N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                gA[j * N + i] = j * N + i;
			}
        }

        printf("A is:\n");
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                printf("%f ", gA[j * N + i]);
			}
            printf("\n");
        }
    }

    // create local A on every process which we'll process
    double* lA = new double[N / procN * N / procN];

    // create a datatype to describe the subarrays of the gA array
    int sizes[2]    = {N, N}; // gA size
    int subsizes[2] = {N / procN, N / procN}; // lA size
    int starts[2]   = {0,0}; // where this one starts
    MPI_Datatype type, subarrtype;
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &type);
    MPI_Type_create_resized(type, 0, N / procN * sizeof(double), &subarrtype);
    MPI_Type_commit(&subarrtype);
	if (rank ==0)
	{
		cout<< "Sizes/ subs "<<endl;
		cout<< sizes[0] << "\t"<< sizes[1] << endl;
		cout<< subsizes[0] << "\t"<< subsizes[1] << endl;
	}
	cout<< endl;
    // compute number of send blocks
    // compute distance between the send blocks
    int sendcounts[procN * procN];
    int displs[procN * procN];

    if (rank == 0) {
        for (int i = 0; i < procN * procN; i++) {
            sendcounts[i] = 1;
        }
        int disp = 0;
	cout << "disp"<<endl;
        for (int i = 0; i < procN; i++) {
            for (int j = 0; j < procN; j++) {
                displs[i * procN + j] = disp;
                disp += 1; 
		cout<< displs[i*procN+j]<< "\t";
            }
		cout<< endl;
            disp += ((N / procN) - 1) * procN;
        }
	cout<< endl;
	cout<< "procN*procN" <<endl;
	cout<< procN*procN << endl;
    }

    // scatter global A to all processes
    MPI_Scatterv(gA, sendcounts, displs, subarrtype, lA,
                 N*N/(procN*procN), MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    // print local A's on every process
    for (int p = 0; p < size; p++) {
    	if (rank == p) {
    		printf("la on rank %d:\n", rank);
            for (int i = 0; i < N / procN; i++) {
                for (int j = 0; j < N / procN; j++) {
                    printf("%f ", lA[j * N / procN + i]);
                }
                printf("\n");
            }
        }
    	MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // write new values in local A's
    for (int i = 0; i < N / procN; i++) {
        for (int j = 0; j < N / procN; j++) {
            lA[j * N / procN + i] = rank;
        }
    }

    // gather all back to master process
    MPI_Gatherv(lA, N*N/(procN*procN), MPI_DOUBLE,
                gA, sendcounts, displs, subarrtype,
                0, MPI_COMM_WORLD);

    // print processed global A of process 0
    if (rank == 0) {
        printf("Processed gA is:\n");
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                printf("%f ", gA[j * N + i]);
            }
            printf("\n");
        }
    }

    MPI_Type_free(&subarrtype);

    if (rank == 0) {
        delete gA;
    }

    delete lA;

    MPI_Finalize();

    return 0;
}
