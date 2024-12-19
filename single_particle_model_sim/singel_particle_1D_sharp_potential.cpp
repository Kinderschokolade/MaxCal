#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream> 
#include <random>
#include <string> 

using namespace std;

#include "init.h"
#include "force.h"
#include "force_delta.h"  
#include "find_shift.hpp"

void read_flags (int argc, char* argv[0], 
	vector <double> & U , vector <double> & T,int & lvl, double & gamma, string & ident, 
	double & fullT, double & dt, int & potential, double & k, double & f,
	bool & ierr,int & it, int & ms)
{ 
	ostringstream arg;

	gamma = 1;
	ident = "x";
	fullT = 1;
	dt = 0.01;
	potential = 1;
	k = 50;
	f = 0;	
	ms = 3;

	ierr =0;

	if(	std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")
	{
		cout<<"-lvl \t number of states" << endl;
		cout<<"-Ui \t Potential at pos i - default 0" << endl;	
		cout<<"-Ti \t Temperature at pos i - default 1" << endl;	
		cout<<"-gamma \t fricition coeff - default 1" << endl;
		cout<<"-o \t recognition" << endl;
		cout<<"-fullT \t simulation time " << endl;
		cout<<"-dT \t time steps" << endl;
		cout<<"-pot \t 1 - deltapot" << endl;
		cout<<" \t 2 - numeric delta pot " <<endl;
		cout<<"-k \t sharpness of numerical heaviside fct" <<endl;
		cout<<"-ms \t number of microstates - default 3" <<endl;
		ierr =1;
	} 
	for(int i=0; i<argc; i++) 
	{
		if (std::string(argv[i]) == "-lvl")  lvl = (int) atof(argv[i+1]);
	}
	U.resize(lvl,0);
	T.resize(lvl,1);
	
	for(int i=0; i<argc; i++) 
	{
		if (argv[i][1]== 'T' || argv[i][1]== 'U')
		{
			for (int j=0 ; j< U.size() ; j++)
			{
				arg.str("");
				arg.clear();
				arg<<"-U"<<j+1;
     				if(std::string(argv[i]) == arg.str()) 
				{
		      			U[j] = (double)atof(argv[i+1]);
				}
     				arg.str("");
				arg.clear();
				arg<<"-T"<<j+1;
				if(std::string(argv[i]) == arg.str()) 
				{
	      				T[j] = (double)atof(argv[i+1]);
				}
			}
		}	
		else if (std::string(argv[i]) == "-gamma") 
		{
	      		gamma  = (double)atof(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-ms") 
		{
	      		ms  = (double)atof(argv[i+1]);
		}
	     	else if(std::string(argv[i]) == "-o") 
		{
	      		ident  = std::string(argv[i+1]);
		}
	     	else if(std::string(argv[i]) == "-fullT") 
		{
	      		fullT  = (double)atof(argv[i+1]);
		}
	     	else if(std::string(argv[i]) == "-dT") 
		{
	      		dt = (double)atof(argv[i+1]);
		}
		else if(std::string(argv[i]) == "-pot") 
		{
			if (std::string(argv[i+1]) == "1") potential = 1 ; // delta pot
			if (std::string(argv[i+1]) == "2") potential = 2 ; // tanh approx
			if (std::string(argv[i+1]) == "3") potential = 3 ; // external force
		}
		else if(std::string(argv[i]) == "-k" )
		{
			k = (double) atof(argv[i+1]);
		}
		else if(std::string(argv[i]) == "-it" )
		{
			it = (int) atof(argv[i+1]);
		}
		else if(std::string(argv[i]) == "-f")
		{
			f = (double) atof(argv[i+1]);
		}
		else if(std::string(argv[i]) == "-lvl")
		{ //nothing, only to avoid "does not exist" 
		}
     		else if(argv[i][0] == '-' ) 
		{
	      		cout<< argv[i] << " does not exist" << endl;
		}

	}
	if (ident == "x")ierr =1;
}

int find_ms(double & x, int ms)
{
	return ( floor( ms * x ) );
}


int main(int argc, char *argv[])
{
	const int dim = 1;
	int lvl ; // # of T,U regions
	int ms; // # microstates
	double gamma = 5;
	double dt;
	double fullT;
	double int_steps;
	double Ekin=0, pav=0;
	int blocknum = 10;
	int blocksize ;

	int number_counts=200;

	vector <double> U;	
	vector <double> T;
	
	vector <double> box(dim,1);
	vector <double> v, x, xo;
	vector <double> f(dim,0);

	vector < vector <int> > block_hist;
	
	double E = 0;
	double Eav;
	int pos;
	bool err = false;
	bool bc = false;
	bool ierr;
	string ident;
	int block;
	double k;
	int pos_old;
	int potential;
	double ext_f;

	int countfac =3000; 
	int countfac_v[number_counts];  
	for (int i =0 ; i < number_counts ; i++)
	{
		countfac_v[i] = i*5; // adjust...
	}	
	countfac_v[0] = 1;


	int write_freq = 1000;

	int it;
	double int_cut =0.1; // integration cut off for shift.

	read_flags(argc, argv, U ,T ,lvl, gamma, ident,fullT, dt,potential,k,ext_f, ierr,it,ms);
	vector <double> p(dim,0);
	vector <double> Perr(ms,0);
	vector < vector < vector <double> > > Tcount (blocknum, vector<vector <double> > 
		 (ms,vector <double> (ms,0)));

	vector <int> hist(ms,0);
	vector < vector <vector <int> > > Count(number_counts, 
		vector <vector <int> > (ms, vector <int>(ms,0)));

	vector < vector <int> > histo(number_counts, vector <int> (ms,0));
	vector < vector <double > > transition (ms, vector<double> (ms,0));
	vector < vector <double > > tr_err (ms, vector<double> (ms,0));
	vector < vector <double > > W (ms, vector<double> (ms,0));
	vector < vector <double > > Werr(ms, vector<double>(ms,0));
	vector < vector <double > > S (ms, vector<double> (ms,0));
	vector < vector <double > > A (ms, vector<double> (ms,0));

	double v_cut = 1./(double)lvl;

	
	vector <double> xc(lvl+1,0);

	int_steps = fullT / dt;
	blocksize = int_steps / blocknum ;
	Tcount.resize(blocknum);
	block_hist.resize(blocknum);
	for(int i=0 ; i< Tcount.size(); i++) 
	{
		block_hist[i].resize(ms,0);
		Tcount[i].resize(ms);
		for(int j=0 ; j < lvl ; j++) Tcount[i][j].resize(ms,0);
	}
	
	
	default_random_engine generator;

	init_single_particle( x, xo, box, v, dim);
	vector <double> beta(lvl,0); 	
	beta[0] = 1./T[0];
	for (int i=1 ; i< lvl ; i++)
	{
		beta[i] = 1./T[i];
		xc[i] = double(i)/lvl+findroot(beta[i-1],beta[i],U[i-1],U[i],k,k,int_cut);
	}
	xc[0] = 0.+findroot(beta[lvl-1],beta[0],U[lvl-1],U[0],k,k,int_cut);
	xc[lvl] = 1.+xc[0];
		
	ofstream file;
	file.open("/data/isilon/bause/single_particle/trajectory_"+ident+".dat");


	ofstream param_file;
	param_file.open("/data/isilon/bause/single_particle/param_"+ident+".dat");

	for (int i=0 ; i<lvl ; i++)
	{
		param_file <<"U"<<i<< "\t\t" << U[i] << endl;
	}
	for (int i=0 ; i<lvl ; i++)
	{
		param_file <<"T"<<i<< "\t\t" << T[i] << endl;
	}	
	for (int i=0 ; i<lvl+1 ; i++)
	{
		param_file <<"xc"<<i<< "\t\t" << xc[i] << endl;
	}

	param_file <<"gamma\t\t" << gamma <<endl;
	param_file <<"potential\t" << potential << endl;
	param_file <<"k \t\t" << k << endl;
	param_file <<"extf\t" << ext_f << endl;
	param_file <<"dT \t\t" << dt << endl;
	param_file <<"microstates \t" << ms << endl;
	double sigmaf;

	x[0] = 0.5;
	xo[0] = 0.5;
	pos = find_ms(x[0],ms);
	pos_old= pos;

	vector <int> survival_hist;
	vector <int> prev_pos(number_counts,pos_old);
	int survival = 0;

	for (unsigned long long i=0 ; i < int_steps; i++)
	{
		force_nstate(f, x, xo, box, U, T, 
			gamma, dt, generator,p,k,bc, err, sigmaf, xc, ext_f);
		pos =find_ms(x[0],ms);
		survival++;
		if (i% countfac == 0){
			hist[pos]++;
			block = i / blocksize;
			block_hist[block][pos]++;
			Tcount[block][pos_old][pos]++;
			pos_old = pos;
		}
		if (pos != pos_old)
		{
			if (survival_hist.size()<survival) 
				survival_hist.resize(survival+1);
			survival_hist[survival]++;
			survival_hist.push_back(survival);
			survival =0;  
		}

		for (int c =0 ; c < number_counts ; c++ )
		{
			if (i % countfac_v[c] == 0)
			{
				histo[c][pos]++;
				Count[c][prev_pos[c]][pos]++;
				prev_pos[c] = pos;
			}
		}
		

		
		Ekin = (Ekin *i + 1./2. * p[0] * p[0]         )/(double)(i+1);
		Eav =  (Eav    *i + 1./2. * p[0] * p[0] + U[pos])/(double)(i+1);
		pav = (pav*i + p[0])/(double)(i+1);
		if (i % write_freq == 0)
		{
			file << i << "\t" << x[0] << "\t" << p[0] << "\t" 
				<< f[0] << "\t" << sigmaf << "\t" << Ekin << "\t" 
				<< Eav << "\t" << pav << endl;
		}
	}
	ofstream count_file;
	count_file.open("/data/isilon/bause/single_particle/count_"+ident+".dat");
	
	count_file << "# ";	
	for (int i=0 ; i< number_counts ; i++)
	{
		count_file << countfac_v[i] << "\t";
	}
	count_file << endl;

	for (int k=0 ; k< ms ; k++)
	{
		for (int j=0 ; j< ms ; j++)
		{
			count_file<<"C"<<k+1<<j+1<<"\t";
			for (int i=0 ; i< number_counts ; i++)
			{
				count_file << Count[i][j][k] << "\t";
			}
			count_file << endl;
		}
	}
	for (int k=0 ; k< ms ; k++)
	{
		count_file << "H"<<k+1 <<"\t";
		for (int i=0 ; i< number_counts ; i++)
		{
			count_file << histo[i][k] << "\t";
		}
		count_file << endl;
	}

	ofstream survival_file;
	survival_file.open("/data/isilon/bause/single_particle/survival_"+ident+".dat");

	for (int i=0 ; i< survival_hist.size() ; i++)
	{
		survival_file <<i << "\t" << survival_hist[i] << endl;
	}
	
	ofstream res_file;
	res_file.open("/data/isilon/bause/single_particle/res_"+ident+".dat");
	
	//  rewrite if neccessary
	res_file << "#log(hist[i]/hist[j])\tUi-Uj"<<endl;
	res_file<<"logl\t";
	for (int i=0 ; i< ms ; i++)
	{
		for (int j=i+1 ; j< ms ; j++)
		{
		res_file<<log(hist[i]/ (double) hist[j]) << "\t" << U[i] - U[j] << "\t";
		}
	}
	res_file << endl;
	double tr, Wtr, Pval;
	vector <double> Pstat(ms,0);

	for (int b=0 ; b< blocknum ; b++)
	{
		for (int i=0 ; i < ms ; i++)
		{
			//res_file <<"T \t";
			for (int j=0 ; j< ms ; j++) 
			{
				tr = Tcount[b][i][j]  / block_hist[b][i];
				Wtr = Tcount[b][i][j] / blocksize;
				tr_err[i][j] = (tr *tr + b*tr_err[i][j] ) / (double)(b+1);
				transition[i][j] = (tr+b* transition[i][j])/(double) (b+1);
				W[i][j] = (Wtr +b* W[i][j])/(double)(b+1);
				Werr[i][j] =  (Wtr * Wtr + b* Werr[i][j])/(double)(b+1);
			}
			Pval = block_hist[b][i] /(double) blocksize;
			Pstat[i] = (Pval + b* Pstat[i]) / (double) (b+1);
			Perr[i] = (Pval*Pval + b* Perr[i]) / (double)(b+1);
		}
	}
	
	ofstream block_file;
	block_file.open("/data/isilon/bause/single_particle/block_"+ident+".dat");

	for (int i=0 ; i < ms ; i++)
	{
		for (int j=0 ; j< ms ; j++) 
		{
			block_file<<"T" << i << j << endl;
			for (int b=0 ; b<blocknum ; b++)
			{
				block_file <<"\t" <<  Tcount[b][i][j] 
					 / block_hist[b][i]<<endl;
			}
			block_file<<endl;
		}
	}

	for (int i=0 ; i < ms ; i++)
	{
		res_file <<"Tval \t";
		for (int j=0 ; j < ms ; j++) 
		{
			res_file << transition[i][j] << "\t";
			tr_err[i][j]=sqrt(tr_err[i][j]-transition[i][j]*transition[i][j])
					/sqrt(blocknum);
			Werr[i][j] = sqrt(Werr[i][j] - W[i][j] * W[i][j])
					/sqrt(blocknum);
		}
		Perr[i] = sqrt(Perr[i] - Pstat[i] *Pstat[i]) / sqrt(blocknum);
		res_file << endl;
	}
	
	for (int i=0 ; i < ms ; i++)
	{
		res_file <<"Terr\t";
		for (int j=0 ; j< ms ; j++) 
		{
			res_file << tr_err[i][j] << "\t";
		}
		res_file << endl; 
	}
	
	for (int i=0 ; i < ms ; i++)
	{
		res_file <<"Count\t";
		for (int j=0 ; j < ms ; j++) 
		{
			res_file << W[i][j] << "\t";
		}
		res_file << endl;
	}
	
	for (int i=0 ; i < ms ; i++)
	{
		res_file <<"Aval \t";
		for (int j=0 ; j< ms ; j++) 
		{
			S[i][j] = (W[i][j] + W[j][i])/2.;
			A[i][j] = (W[i][j] - W[j][i])/2. ;
			res_file << A[i][j] << "\t";
		}
		res_file << endl; 
	}
	
	for (int i=0 ; i < ms ; i++)
	{
		res_file <<"Aerr \t";
		for (int j=0 ; j< ms ; j++) 
		{
			res_file<<1./2 *( fabs(fabs(Werr[i][j])-fabs(Werr[j][i])) )<< "\t";
		}
		res_file << endl; 
	}

	
	for (int i=0 ; i < ms ; i++)
	{
		res_file <<"Sval \t";
		for (int j=0 ; j< ms ; j++) 
		{
			res_file << S[i][j] << "\t";
		}
		res_file << endl;
	}
	
	double Sprod =0;
	for (int i=0 ; i < ms ; i++)
	{
		res_file <<"Serr \t";
		for (int j=0 ; j< ms ; j++)
		{ 
			res_file<<1./2 *( fabs(fabs(Werr[i][j])+fabs(Werr[j][i])) )<< "\t";
			if (S[i][j] > 0 || A[i][j] > 0)
			{
				Sprod += 1./2. * (W[i][j]-W[j][i]) 
					* log((S[i][j]+A[i][j])/(S[j][i]-A[i][j]));
			}
		}
		res_file<< endl;
	}
	res_file << "Sprod" << "\t"<< Sprod<< endl;
	
	for (int i=0 ; i< ms ; i++)
	{
		res_file <<"Pval\t" << Pstat[i] << endl;
	}
	for (int i=0 ; i< ms ; i++)
	{
		res_file <<"Perr\t" << Perr[i] << endl;
	}


}









 
