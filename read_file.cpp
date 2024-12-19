#include <string> 
#include <iostream>
#include <fstream>
#include <sstream> 

using namespace std;


int main(int argc, char *argv[])
{
	double readv;
	string line;
	ifstream myfile("test.dat",ios::binary);
	if (myfile.is_open())
	{
		myfile.read(&readv,sizeof(readv));
		//cout<< readv[0] << "\t" << readv[1] << endl;
	}
	else cout << "Unable to open file"; 


}
