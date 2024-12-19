#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;


int main(int argc, char *argv[])
{
	int nx = 2;
	double *data = new double [nx]; 
	double *data_read = new double [nx]; 


	FILE *file = fopen("rest.dat", "rb");
	fread(data_read, sizeof(double), nx, file );
	fclose(file);
	
	for ( int i=0 ; i < 10 ; i++)
	{
		cout<< data_read[i] << endl;
	}

}

