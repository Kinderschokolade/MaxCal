#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

using namespace std;


void read_flags (int argc, char* argv[0], 
	vector <double> & G , vector <double> & k, 
	vector <vector <double> >  & alpha, string & ident, 
	bool & ierr)
{ 
	ostringstream arg;
	for (int i= 0 ; i< G.size() ; i++) G[i] = 1;
	for (int i= 0 ; i< k.size() ; i++) k[i] = 0.2;
	for (int i= 0 ; i< alpha.size() ; i++)
	{ 
		for (int j=0 ; j<alpha.size() ; j++) {alpha[i][j] = 0;}
	}

	ident = "x";
	ierr =0;

	if(	std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")
	{
		cout<<"-Gi \t free Energy at pos i - default 1" << endl;	
		cout<<"-ki \t equ. exchange rate from i to i+1 - default 0.2" << endl;	
		cout<<"-aij \t driving rate between i and j - default 0" << endl;
		cout<<"-o \t recognition" << endl;
		ierr =1;
	} 
	for(int i=0; i<argc; i++) 
	{
		if (argv[i][1]== 'G' || argv[i][1]== 'k')
		{
			for (int j=0 ; j< G.size() ; j++)
			{
				arg.str("");
				arg.clear();
				arg<<"-G"<<j+1;
     				if(std::string(argv[i]) == arg.str()) 
				{
		      			G[j] = (double)atof(argv[i+1]);
				}
     				arg.str("");
				arg.clear();
				arg<<"-k"<<j+1;
				if(std::string(argv[i]) == arg.str()) 
				{
	      				k[j] = (double)atof(argv[i+1]);
				}
			}
		}	
		else if (argv[i][1]== 'a')
		{
			for (int j=0 ; j< alpha.size() ; j++)
			{
				for (int m=0 ; m< alpha[0].size() ; m++)
				{
					arg.str("");
					arg.clear();
					arg<<"-a"<<j+1<<m+1;
     					if(std::string(argv[i]) == arg.str()) 
					{
		      				alpha[j][m]=(double)atof(argv[i+1]);
					}
				}
			}
		}
	     	else if(std::string(argv[i]) == "-o") 
		{
	      		ident  = std::string(argv[i+1]);
		}
	     	else if(argv[i][0] == '-') 
		{
	      		cout<< argv[i] << " does not exist" << endl;
		}

	}

	if (ident == "x")ierr =1;
}


int main(int argc, char* argv[])
{
	int N = 1000 ; //Number of steps
	int dim = 3;
	int numit = 10;
	double deltaalpha = 0.001;

	int old_pos;
	double norm;
	double S_prod;
	double Q;
	bool ierr;
	string ident;
	double v_length = pow(dim,N);

	vector <double> p(dim,0);
	vector <double> p_old(dim,0);
	vector <double> G(dim,0);
	vector <double> P_stat(dim,0);
	vector <double> prob_norm(dim,0);

	vector < vector < double > > I(dim, vector <double> (dim,0));
	vector < vector < double > > Is(dim, vector <double> (dim,0));
	vector < vector < double > > W(dim, vector <double> (dim,0));
	vector < vector < double > > R(dim, vector <double> (dim,0));
	vector < vector < double > > k(dim, vector <double> (dim,0));
	vector < vector < double > > alpha(dim, vector <double> (dim,0));
	vector < vector < double > > T(dim, vector <double> (dim,0));

	vector <double> kv (dim,0);

	read_flags(argc, argv, G, kv, alpha, ident, ierr);
	//alpha[0][1] = 0.25;
	if (ierr) return EXIT_FAILURE ;	

	k[0][1] = kv[0];
	k[0][2] = kv[1];
	k[0][0] = 1. - k[0][1] - k[0][2] - alpha[0][1] - alpha[0][2];
	k[1][2] = kv[2];
	k[1][0] = exp(G[1] - G[0]) * k[0][1];
	k[1][1] = 1. - k[1][2] - k[1][0] - alpha[1][2] - alpha[1][0];
	k[2][0] = exp(G[2] - G[0]) * k[0][2];
	k[2][1] = exp(G[2] - G[1]) * k[1][2]; 
	k[2][2] = 1. - k[2][1] - k[2][0] - alpha[2][1] - alpha[2][0];

        R[0][0] = 0;
        R[1][0] = exp(G[1]) / k[1][0];
        R[0][1] = R[1][0];
        R[1][1] = 0;
        R[2][0] = exp(G[2]) / k[2][0];
        R[2][1] = exp(G[2]) / k[2][1];
        R[0][2] = R[2][0];
        R[0][1] = R[1][0]; // some mistake here?


	p[0] = 1.;
	p[1] = 0;
	p[2] = 0;	
	p_old = p;	
	ofstream file;
	file.open("data/exact_"+ident+".dat");
       	double R_tot = R[1][0] + R[2][0] + R[2][1];
	double Pnorm;
	file << "#k"<<endl;
	for (int i=0 ; i< dim ; i++)
	{
		for(int j=0; j< dim ; j++) file << k[i][j]<< "\t";
		file << endl;
	}
	file << "#G"<<endl;
	for (int i=0 ; i< dim ; i++)
	{
		file << G[i]<< endl;
	}
	file << "#alpha"<<endl;
	for (int i=0 ; i< dim ; i++)
	{
		for(int j=0 ; j < dim ; j++) {file << alpha[i][j]<< "\t";}
		file<<endl;
	}
	file << endl;
	file << "#R " << endl;
	for (int i=0 ; i< dim ; i++)
	{
		for(int j=0; j< dim ; j++) file << R[i][j]<< "\t";
		file << endl;
	}


	for (int i=0 ; i <dim ; i++)
	{
		for (int j=0 ; j<dim ; j++)
		{
			T[i][j] = (k[i][j] + alpha[i][j]);// / prob_norm[i];
			// T[i][j]   << "\t";
		}
	}

	

	P_stat[0] =  (exp(-G[0] - G[1])*(exp(G[0] + G[1])*
  		  	(exp(G[1])*(alpha[2][0] + alpha[2][1])*k[0][1] + 
      		  	exp(G[0])*(alpha[1][0]*(alpha[2][0] + alpha[2][1]) + 
  	          	alpha[2][0]*(alpha[1][2] + k[1][2]))) + 
  		  	exp(G[2])*(exp(2*G[1])*k[0][1]*k[0][2] + 
       	  	exp(2*G[0])*alpha[1][0]*k[1][2] + 
         		exp(G[0] + G[1])*(k[0][1]*k[1][2] + 
            		k[0][2]*(alpha[1][0] + alpha[1][2] + k[1][2]))
         		)))/(exp(G[1])*k[0][1]*(alpha[0][2] + k[0][2]) + 
    			exp(G[0])*(alpha[1][0]*k[0][2] + 
       		(alpha[0][1] + k[0][1] + k[0][2])*
        		(alpha[1][2] + k[1][2]) + 
       		alpha[0][2]*(alpha[1][0] + alpha[1][2] + k[1][2])));

	P_stat[1] = (exp(G[0] + G[1])*(alpha[0][1]*(alpha[2][0] + alpha[2][1]) + 
       		alpha[2][0]*k[0][1] + 
       		alpha[2][1]*(alpha[0][2] + k[0][1] + k[0][2])) + 
    			exp(G[2])*(exp(G[1])*(alpha[0][1] + k[0][1])*
        		k[0][2] + exp(G[0])*
        		(alpha[0][1] + alpha[0][2] + k[0][1] + k[0][2])*
        		k[1][2]))/
  			(exp(2*G[1])*k[0][1]*(alpha[0][2] + k[0][2]) + 
    			exp(G[0] + G[1])*(alpha[1][0]*k[0][2] + 
       		(alpha[0][1] + k[0][1] + k[0][2])*
        		(alpha[1][2] + k[1][2]) + 
       		alpha[0][2]*(alpha[1][0] + alpha[1][2] + k[1][2])));


	file << "#T"<<endl;
	for (int i=0 ; i< dim ; i++)
	{
		for(int j=0; j< dim ; j++) file << T[i][j]<< "\t";
		file << endl;
	}

	Pnorm = P_stat[0] + P_stat[1] +1 ;	
	P_stat[0] = P_stat[0] / Pnorm;
	P_stat[1] = P_stat[1] / Pnorm;
	P_stat[2] = 1. / Pnorm;

	
	file << "#P_stat"<<endl;;
	for (int i=0 ; i< dim ; i++)
	{
		file << P_stat[i]<< endl;
	}
	vector <double> current (dim,0);
	for (int i=0 ; i< N ; i++)
	{
		p[0] = p_old[0]
		+ p_old[2] * T[2][0] + p_old[1] * T[1][0]
		- p_old[0] * T[0][1] - p_old[0] * T[0][2];
		p[1] = p_old[1]
		+ p_old[2] * T[2][1] + p_old[0] * T[0][1]
		- p_old[1] * T[1][0] - p_old[1] * T[1][2];
		p[2] = p_old[2]
		+ p_old[1] * T[1][2] + p_old[0] * T[0][2]
		- p_old[2] * T[2][0] - p_old[2] * T[2][1];

		current[0] = 1./2. * (p[0] * T[0][1] -  p[1] * T[1][0]);
		current[1] = 1./2. * (p[1] * T[1][2] -  p[2] * T[2][1]);
		current[2] = 1./2. * (p[2] * T[2][0] -  p[0] * T[0][2]);
		
		p_old = p;	
	
	}

	for (int i=0 ; i < dim ; i++)
	{
		for(int j=0 ; j < dim ; j++)
		{
			W[i][j] = T[i][j] * P_stat[i];
		}
	}
	W[0][0] = -T[0][1] * P_stat[0] - T[0][2] * P_stat[0] ;
	W[1][1] = -T[1][0] * P_stat[1] - T[1][2] * P_stat[1] ;
	W[2][2] = -T[2][1] * P_stat[2] - T[2][0] * P_stat[0] ;

	for (int i=0 ; i < dim ; i++)
	{
		for(int j=0 ; j < dim ; j++)
		{
			I[i][j]  = (W[i][j] - W[j][i])/2.;
			Is[i][j] = (W[i][j] + W[j][i])/2.;
		}
	}
	
	S_prod = 0;
	for (int i=0 ; i< dim ; i++)
	{
		for (int j=0 ; j < dim ; j++)
		{
			S_prod += 1./2. * (p[i]*T[i][j]-p[j]*T[j][i])
				* log(T[i][j]/T[j][i]);
		}
	}
	
	file << "#p numerical"<<endl;;
	for (int i=0 ; i< dim ; i++)
	{
		file << p[i]<< endl;
	}

	file << "#physical current"<<endl;;
	for (int i=0 ; i< dim ; i++)
	{
		for (int j=0 ; j< dim ; j++) file << W[i][j]<< "\t";
		file<< endl;
	}

	file << "#Current assymetric"<<endl;;
	for (int i=0 ; i< dim ; i++)
	{
		for (int j=0 ; j< dim ; j++) file << I[i][j]<< "\t";
		file<< endl;
	}

	file << "#Current symmetric"<<endl;;
	for (int i=0 ; i< dim ; i++)
	{
		for (int j=0 ; j< dim ; j++) file << Is[i][j]<< "\t";
		file<< endl;
	}

	file << "#S_prod\t"<< endl << S_prod <<endl;
	file << "#R_tot * I *I\t" <<endl<< R_tot * current[0] * current[0] << endl;

	
	cout<< P_stat[0] << "\t" << P_stat[1] << "\t" << P_stat[2] << endl;

}



