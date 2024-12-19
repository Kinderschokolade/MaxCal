#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream> 
//#include <boost/random.hpp>
#include <string> 
#include <unistd.h>

using namespace std;

//typedef boost::random::mt19937 generator;
#include <random>

#include "init.h"
#include "find_shift.hpp"
//#include "force_delta.h"  
#include "force.h"



void read_flags (int argc, char* argv[0], 
	vector <double> & U , double & T1, double & T2,int & lvl, double & gamma, string & ident, 
	double & fullT, double & dt,  double & sigma, vector < double > & x, double & slope,
	bool & ierr, int & ms, int & write_freq, double & seed)
{ 
	ostringstream arg;

	gamma = 1;
	ident = "x";
	fullT = 1;
	dt = 0.01;
	ms = 3;
	ierr =0;
	write_freq = 500;
	seed = 100;

	if(	std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")
	{
		cout<<"-lvl \t number of basins" << endl;
		cout<<"-Ui \t Potential at basin i - default 0" << endl;	
		cout<<"-sigma \t sigma of basin - default 0" << endl;	
		cout<<"-xi \t position of basin - default 0" << endl;	
		cout<<"-T1 \t Temperature at beginning - default 1" << endl;	
		cout<<"-T2 \t Temperature at end - default 1" << endl;	
		cout<<"-gamma \t fricition coeff - default 1" << endl;
		cout<<"-o \t recognition" << endl;
		cout<<"-fullT \t simulation time " << endl;
		cout<<"-dT \t time steps" << endl;
		cout<<"-ms \t number of microstates - default 3" <<endl;
		cout<<"-wf \t writing freq - default 500" <<endl;
		cout<<"-seed \t seed - default 100" << endl;
		ierr =1;
	} 
	for(int i=0; i<argc; i++) 
	{
		if (std::string(argv[i]) == "-lvl")  lvl = (int) atof(argv[i+1]);
	}
	U.resize(lvl,0);
	x.resize(lvl,0);	

	for(int i=0; i<argc; i++) 
	{
		if (argv[i][1]== 'x' || argv[i][1]== 'U')
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
				arg<<"-x"<<j+1;
				if(std::string(argv[i]) == arg.str()) 
				{
	      				x[j] = (double)atof(argv[i+1]);
				}
			}
		}	
		else if (std::string(argv[i]) == "-sigma") 
		{
	      		sigma  = (double)atof(argv[i+1]);
		}		
		else if (std::string(argv[i]) == "-T1") 
		{
	      		T1  = (double)atof(argv[i+1]);
		}	
		else if (std::string(argv[i]) == "-T2") 
		{
	      		T2  = (double)atof(argv[i+1]);
		}	
	
		else if (std::string(argv[i]) == "-slope") 
		{
	      		slope  = (double)atof(argv[i+1]);
		}	
		else if (std::string(argv[i]) == "-gamma") 
		{
	      		gamma  = (double)atof(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-wf")
		{
	      		write_freq  = (int)atof(argv[i+1]);
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
		else if(std::string(argv[i]) == "-seed")
		{
			seed = (double) atof(argv[i+1]);
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

int find_ms(double & x, int & ms, double & cutl, double & cutr, double & box)  // redefine such that boundary corresponds to large ms
{
	int pos;
	if (x < cutl) 
	{
		pos = 0;
	}
	else if (x>cutr)
	{
		pos = ms -1;
	}
	else
	{
		pos = round( (x-cutl)/  (box - 2*cutl) *(ms-1)) ;
	}
	return pos;
}



int main(int argc, char *argv[])
{
	char hostname[HOST_NAME_MAX];
	char username[LOGIN_NAME_MAX];
	char keyname[] = "mbause";
	string path0, path1;
	gethostname(hostname, HOST_NAME_MAX);
	getlogin_r(username, LOGIN_NAME_MAX);
	if (strcmp(username,keyname) == 0 || strcmp(username,"") ==  0)
	{	
		path0 = "/u/mbause/data/";
		path1 = "/u/mbause/data/";
	}
	else
	{
		path0 = "/data/isilon/bause/";
		path1 = "/data/pckr194/bause/";
	}

	//cout<< path0 << endl;
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

	int number_counts=1000;

	vector <double> U;	
	
	vector <double> box(dim,1);
	vector <double> v, x, xo;
	vector <double> f(dim,0);

	vector < vector <int> > block_hist;
	
	double E = 0;
	double Eav;
	int pos;
	bool err = false;
	int bc;
	bool ierr;
	string ident;
	int block;
	double k;
	int pos_old;
	int potential;
	double ext_f;
	double padd=0.;

	int countfac =100; 
	int countfac_v[number_counts];  
	for (int i =0 ; i < number_counts ; i++)
	{
		countfac_v[i] = i; // adjust...maybe.
	}	
	countfac_v[0] = 1;
	vector <int> bc_part (number_counts,0);

	int write_freq;

	double int_cut =0.1; // integration cut off for shift.
	double seed;

	vector <double> xbas;
	double sigma,slope;
	double T1, T2;
	read_flags(argc, argv, U ,T1,T2 ,lvl, gamma, ident,fullT, dt,sigma,xbas,slope, ierr,ms, write_freq,seed);
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
	
	//generator rng(seed);	

	mt19937 rng(seed);
	init_single_particle( x, xo, box, v, dim);
	x[0] = xbas[0];
	xo[0] = xbas[0];
		
	ofstream file;
	string link; // use to build paths
//	cout<< it << "\t" << xc[0] << "\t" <<xc[1]<<"\t"<<xc[2]<<"\t"<<xc[3]<< endl;
	if (write_freq >0)
	{	
		link = path1+"Tgrad/trajectory_"+ident+".dat";
		file.open(link.c_str());
	}
	else
	{
		file.open("/data/pckr194/bause/Tgrad/mumpitz.dat");
	}



	link = path0+"Tgrad/param_"+ident+".dat";

	ofstream param_file;
	param_file.open(link.c_str());

	for (int i=0 ; i<lvl ; i++)
	{
		param_file <<"U"<<i<< "\t\t" << U[i] << endl;
	}
	for (int i=0 ; i<lvl+1 ; i++)
	{
		param_file <<"xbas"<<i<< "\t\t" << xbas[i] << endl;
	}

	param_file <<"sigma\t\t" << sigma <<endl;
	param_file <<"slope\t\t" << slope <<endl;
	param_file <<"T\t\t" << T1<< "\t" << T2 <<endl;
	param_file <<"gamma\t\t" << gamma <<endl;
	param_file <<"potential\t" << potential << endl;
	param_file <<"dT \t\t" << dt << endl;
	param_file <<"microstates \t" << ms << endl;
	param_file <<"seed \t" << seed << endl;
	param_file <<"Time \t" << fullT << endl;
	double sigmaf;
	double fav, f2av;

	double cutl = xbas[0] - sigma;
	double cutr = box[0] - cutl;

	pos = find_ms(x[0],ms,cutl,cutr,box[0]);
	pos_old= pos;

	
	vector <int> prev_pos(number_counts,pos_old);

	double prec_F = 10000;
	double xmul = box[0] / prec_F;

	ofstream pot_file;
	ofstream potms_file;
	link = path0+"Tgrad/potential_"+ident+".dat";
	pot_file.open(link.c_str());

	link = path0+"Tgrad/potential_ms_"+ident+".dat";
	potms_file.open(link.c_str());
	

	vector < double> Uvec = get_U(prec_F,U,xbas,box,sigma,slope);
	vector < double> Fvec = get_F(prec_F,U,xbas,box,sigma,slope);

	double start_ms, end_ms;
	start_ms =0;
	vector <double > ms_x(ms,0);
	int ms_count =0,pp;
	int poscu= 0,posold =0;
	double xcu=0. ,xold=0.;

	for (int i=1 ; i < prec_F ; i++)
	{
		xold = xcu;
		posold = poscu;
		xcu = i / (double)(prec_F);
		poscu = find_ms(xcu,ms,cutl, cutr,box[0]);
		pot_file << xcu << "\t" << Uvec[i] << "\t" <<Fvec[i]  <<endl;
		//cout << xcu<< "\t"<< poscu << "\t" << Uvec[i] << "\t" <<Fvec[i]  <<endl;
		if (posold != poscu)
		{	
			end_ms = xold;
			ms_x[ms_count] = (end_ms + start_ms) / 2.;
			start_ms = xcu;
			pp = round( ms_x[ms_count] * box[0] *prec_F);
			potms_file <<ms_x[ms_count] << "\t" << Uvec[pp]  <<endl;
			ms_count +=1;
		}
	}
	ms_x[ms-1] = (end_ms + box[0]) / 2.;
	pp = round( ms_x[ms-1] * box[0] *prec_F);
	potms_file <<ms_x[ms-1] << "\t" << Uvec[pp]  <<endl;

	potms_file.close();
	int progress = int_steps / 10.;

	for (unsigned long long i=0 ; i < int_steps; i++)
	{
		//cout<< i << endl;
		force_Tgrad(f, x, xo, box, T1,T2, gamma, dt, rng,p, err, sigmaf, Fvec,xmul);
		if (write_freq != 0)
		{
			if (i % write_freq == 0)
			{
				file <<x[0]<<"\t"<< p[0] <<"\t"<<f[0]/gamma*dt << "\t"<< sigmaf << endl;
			}
		}
		if ( x[0] > 0 and x[0] < 1.)
		{
			pos =find_ms(x[0],ms,cutl,cutr,box[0]);
		}
		else if(x[0] < 0.)
		{
			cout<<x[0] << endl; 
			x[0] = 0.001;
			pos = 0;	
	}
		else
		{
			cout<<x[0] << endl; 
			x[0] = 0.999;
			pos =ms -1; // hoffe ich
		}
		if (i% countfac == 0){
			hist[pos]++;
			block = i / blocksize;
			block_hist[block][pos]++;
			Tcount[block][pos_old][pos]++;
		}
/*		if (x[0] > 0.95*box[0] || x[0] < 0 || x[0] > 1.)
		{
			//cout<< x[0] << endl;
			xo[0] = 0.1;
			x[0] = 0.1;
			pos =find_ms(x[0],ms);
			for (int c =0 ; c < number_counts ; c++ )
			{
				prev_pos[c] = pos;
			}
			file << endl; // empty line to indicate restart
		}*/
	//	else{
			for (int c =0 ; c < number_counts ; c++ )
			{
				//v_part[c] = v_part[c] + p[0]; // this does not work.. why?
				if (i % countfac_v[c] == 0)
				{
					Count[c][prev_pos[c]][pos]++;
					histo[c][pos]++;
					prev_pos[c] = pos;
				}
			}
			//file<< "av" <<endl;
			//if (i % 100000==0) cout<< i <<"\t" << int_steps <<  endl;
	
//			Ekin = (Ekin *i + 1./2. * p[0] * p[0]         )/(double)(i+1);
//			Eav =  (Eav    *i + 1./2. * p[0] * p[0] + U[pos])/(double)(i+1);
//			pav = (pav*i + p[0])/(double)(i+1);
//			fav = (fav*i + f[0])/(double)(i+1);
//			f2av =  (f2av*i + f[0]*f[0])/(double)(i+1);
//	}// of else
		//file<< "endl" <<endl;
		if (i % progress == 0)
		{
			cout<< i / (double) int_steps << endl;
		}
	}
	cout<<"Sim finished"<<endl;
	
	ofstream count_file;
	link = path0+"Tgrad/count_"+ident+".dat";
	count_file.open(link.c_str());
	

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
	
}









 
