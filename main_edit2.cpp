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
	vector <double> & U , vector <double> & T,int & lvl, double & gamma, string & ident, 
	double & fullT, double & dt, int & potential, double & f,
	bool & ierr,int & it, int & ms, int & write_freq, double & seed)
{ 
	ostringstream arg;

	gamma = 1;
	ident = "x";
	fullT = 1;
	dt = 0.01;
	potential = 1;
	f = 0;	
	ms = 3;
	ierr =0;
	write_freq = 500;
	seed = 100;

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
		else if(std::string(argv[i]) == "-pot") 
		{
			if (std::string(argv[i+1]) == "1") potential = 1 ; // cos approx.
			if (std::string(argv[i+1]) == "2") potential = 2 ; // tanh approx
			if (std::string(argv[i+1]) == "3") potential = 3 ; // external force
		}
		else if(std::string(argv[i]) == "-it" )
		{
			it = (int) atof(argv[i+1]);
		}
		else if(std::string(argv[i]) == "-f")
		{
			f = (double) atof(argv[i+1]);
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

int find_ms(double & x, int ms)
{
	return ( floor( ms * x ) );
}

int force_nstate_overdamped(vector <double> & f, vector <double> & x, vector <double> & xo, 
	vector <double > & box, vector <double > & U, vector <double> & T, 
	double & gamma, double & dt, mt19937 & generator,
	vector <double> & p ,double & k, int & bc, bool & err, double & sigmap, 
	vector <double> & xc, double & ext_f, int & potential
	)
{
	int pos;
	int dim = f.size();
	double mean = 0;
	double sigma;	

	int lvl = xc.size();


	
//	const auto& F = F_tanh;	//stanni!
	
//	if (potential == 1)
//	{
//		const auto& F =F_sin;
//	}
	vector <double> xnew(dim,0);

	//normal_distribution < double > gauss(mean, sigma);
	//uniform_real_distribution <double> lin(-0.5,0.5);
	normal_distribution<double> gauss(0.0,1.0);

	for (int i=0; i < dim ; i++)
	{
		pos = find_pos(x[i],xc);
		//pos = floor(3*x[i]);
		sigma = sqrt(2. * T[pos] * dt / gamma);
		sigmap = sigma * gauss(generator);
		//xnew[i] = 2*x[i] - xo[i] + f[i]*dt*dt; 
		xnew[i] = x[i] + (F_tanh(x[i],k,U,xc) + ext_f)/gamma *dt + sigmap;
		p[i] = (xnew[i] - x[i]) / dt;
		periodic_boundary(xnew[i],box[i],bc);
	}	

	xo[0] = x[0] ;
	x[0]  = xnew[0] ;
	
	if (x[0] > box[0] || x[0] < 0) 
	{
		cout<< "out of bound err" <<endl;
		err = 1;
	}
	return 0;


}



int main(int argc, char *argv[])
{
	char hostname[HOST_NAME_MAX];
	char username[LOGIN_NAME_MAX];
	char keyname[] = "mbause";
	string path0, path1;
	gethostname(hostname, HOST_NAME_MAX);
	getlogin_r(username, LOGIN_NAME_MAX);
	if (strcmp(username,keyname) == 0 )
	{	
		path0 = "/u/mbause/data/";
		path1 = "/u/mbause/data/";
	}
	else
	{
		path0 = "/data/isilon/bause/";
		path1 = "/data/pckr194/bause/";
	}
	cout<<username << "\t" << keyname << "\t "<<  strcmp(username,keyname) << endl;
	cout<< path0 << endl;	

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
	vector <double> T;
	
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

	int it;
	double int_cut =0.1; // integration cut off for shift.
	double seed;

	read_flags(argc, argv, U ,T ,lvl, gamma, ident,fullT, dt,potential,ext_f, ierr,it,ms, write_freq,seed);
	
	vector <double> p(dim,0);
	vector <double> Perr(ms,0);
	vector < vector < vector <double> > > Tcount (blocknum, vector<vector <double> > 
		 (ms,vector <double> (ms,0)));

	vector <int> hist(ms,0);
	vector < vector <vector <int> > > Count_pl(number_counts, 
		vector <vector <int> > (ms, vector <int>(ms,0)));
	
	vector < vector <vector <int> > > Count_mi(number_counts, 
		vector <vector <int> > (ms, vector <int>(ms,0)));

	vector < vector <int> > histo(number_counts, vector <int> (ms,0));
	vector < vector <double > > transition (ms, vector<double> (ms,0));
	vector < vector <double > > tr_err (ms, vector<double> (ms,0));
	vector < vector <double > > W (ms, vector<double> (ms,0));
	vector < vector <double > > Werr(ms, vector<double>(ms,0));
	vector < vector <double > > S (ms, vector<double> (ms,0));
	vector < vector <double > > A (ms, vector<double> (ms,0));

	double v_cut = 1./(double)lvl;


	double conf = 0.999; //How good is U approx
	k = - lvl * log(1./conf -1);

	//if (potential ==1 ) ierr = 1;	// not defined
	
	vector <double> xc(lvl+1,0);
	//if (ierr) return EXIT_FAILURE;

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
	vector <double> beta(lvl,0); 	
	beta[0] = 1./T[0];
	for (int i=1 ; i< lvl ; i++)
	{
		beta[i] = 1./T[i];
		xc[i] = double(i)/lvl+findroot(beta[i-1],beta[i],U[i-1],U[i],k,k,int_cut);
		if (potential == 1)
		{
			xc[i] = double(i) / lvl;
		}
	}
	xc[0] = 0.+findroot(beta[lvl-1],beta[0],U[lvl-1],U[0],k,k,int_cut);
	xc[lvl] = 1.+xc[0];
	if (potential == 1)
	{
		xc[0] = 0.;
		xc[lvl] = 1.;
	}
		
	ofstream file;
	string link; // use to build paths
//	cout<< it << "\t" << xc[0] << "\t" <<xc[1]<<"\t"<<xc[2]<<"\t"<<xc[3]<< endl;
	if (write_freq ==0)
	{	
		link = path1+"single_particle/trajectory_"+ident+".dat";
		file.open(link.c_str());
	}
	else
	{
		file.open("/data/pckr194/bause/single_particle/mumpitz.dat");
	}
	ofstream param_file;
	link = path0+"single_particle/param_"+ident+".dat";
	param_file.open(link.c_str());

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
	param_file <<"seed \t" << seed << endl;
	double sigmaf;
	double fav, f2av;

	x[0] = 0.5;
	xo[0] = 0.5;
	pos = find_ms(x[0],ms);
	pos_old= pos;

	vector <int> survival_hist;
	vector <int> prev_pos(number_counts,pos_old);
	int survival = 0;
	int surv_max = 500;
	survival_hist.resize(surv_max);

	ofstream pot_file;
	ofstream potms_file;
	link = path0+"single_particle/potential_"+ident+".dat";
	pot_file.open(link.c_str());

	link = path0+"single_particle/potential_ms_"+ident+".dat";
	potms_file.open(link.c_str());



	vector <double> Fvec = get_F(200,  k, U, xc,potential);
	vector <double> Uvec = get_U(200,  k, U, xc,potential);

	vector < double> Uvec_ms = get_U(ms,k,U,xc,potential);
	for (int i=0 ; i < 200 ; i++)
	{
		pot_file <<i/200. << "\t" << Uvec[i] << "\t" <<Fvec[i]  <<endl;
	}

	for (int i=0 ; i < ms ; i++)
	{
		potms_file <<i/(double)ms << "\t" << Uvec_ms[i]  <<endl;
	}

	potms_file.close();

//	exit(EXIT_FAILURE);
	for (unsigned long long i=0 ; i < int_steps; i++)
	{
		//file<< "force" <<endl;
		force_nstate_overdamped(f, x, xo, box, U, T, gamma, dt, rng,p,k,bc, err, sigmaf, xc, ext_f,potential);
		//file<< "pos" <<endl;
		pos =find_ms(x[0],ms);
		//pos = floor(x[0] * 3);
		//cout<< x[0] << "\t" << pos << endl;
		//file<< "cfac" <<endl;
		if (i% countfac == 0){
			survival++;
			hist[pos]++;
			block = i / blocksize;
			block_hist[block][pos]++;
			Tcount[block][pos_old][pos]++;
			//file << "surv"<< endl;	
			if (pos != pos_old)
			{
				if (survival<surv_max) 
				{
					survival_hist.resize(survival+1);
					survival_hist[survival]++;
					//survival_hist.push_back(survival);
				}
				survival =0;  
			}
			pos_old = pos;
		}
		//file<< "hist" <<endl;
/////////////////
//I think this version is ancient

//////////////


		for (int c =0 ; c < number_counts ; c++ )
		{
			//v_part[c] = v_part[c] + p[0]; // this does not work.. why?
			bc_part[c] += bc;
			if (i % countfac_v[c] == 0)
			{
				histo[c][pos]++;
				if(bc_part[c] > 0)
				{
					Count_pl[c][prev_pos[c]][pos]++;
				}
				else if (bc_part[c] < 0){
					Count_mi[c][prev_pos[c]][pos]++;
				}
				else if (prev_pos[c] < pos) 
				{
					Count_pl[c][prev_pos[c]][pos]++;
				}
				else
				{		
					Count_mi[c][prev_pos[c]][pos]++;
				}
				bc_part[c] =0;
				prev_pos[c] = pos;
			}
		}
		
		//file<< "av" <<endl;
		//if (i % 100000==0) cout<< i <<"\t" << int_steps <<  endl;

		Ekin = (Ekin *i + 1./2. * p[0] * p[0]         )/(double)(i+1);
		Eav =  (Eav    *i + 1./2. * p[0] * p[0] + U[pos])/(double)(i+1);
		pav = (pav*i + p[0])/(double)(i+1);
		fav = (fav*i + f[0])/(double)(i+1);
		f2av =  (f2av*i + f[0]*f[0])/(double)(i+1);
		if (write_freq != 0)
		{
			if (i % write_freq == 0)
			{
				file <<x[0]<<"\t"<< p[0] << endl;
			}
		}
		//file<< "endl" <<endl;
	}
	cout<<"Sim finished"<<endl;
	
	ofstream count_file;
	link = path0+"single_particle/count_"+ident+".dat";
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
			count_file<<"Cpl"<<k+1<<j+1<<"\t";
			for (int i=0 ; i< number_counts ; i++)
			{
				count_file << Count_pl[i][j][k] << "\t";
			}
			count_file << endl;
		}
	}
	for (int k=0 ; k< ms ; k++)
	{
		for (int j=0 ; j< ms ; j++)
		{
			count_file<<"Cmi"<<k+1<<j+1<<"\t";
			for (int i=0 ; i< number_counts ; i++)
			{
				count_file << Count_mi[i][j][k] << "\t";
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
	link = path0+"single_particle/survival_"+ident+".dat";
	survival_file.open(link.c_str());

	for (int i=0 ; i< survival_hist.size() ; i++)
	{
		survival_file <<i << "\t" << survival_hist[i] << endl;
	}
	
	ofstream res_file;
	link = path0+"single_particle/res_"+ident+".dat";
	res_file.open(link.c_str());
	
	//  rewrite if neccessary
//	res_file << "#log(hist[i]/hist[j])\tUi-Uj"<<endl;
//	res_file<<"logl\t";
//	for (int i=0 ; i< ms ; i++)
//	{
//		for (int j=i+1 ; j< ms ; j++)
//		{
//		res_file<<log(hist[i]/ (double) hist[j]) << "\t" << U[i] - U[j] << "\t";
//		}
//	}
//	res_file << endl;
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
	link = path0+"single_particle/block_"+ident+".dat";
	block_file.open(link.c_str());

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

	res_file<< "Ekinav\t"<< Ekin<<  endl;
	res_file<< "Eav\t"<< Eav<<endl ;
	res_file<< "pav\t"<< pav<<endl ;
	res_file<< "fav\t"<< fav<<endl ;
	res_file<< "ferr\t"<< sqrt(f2av - fav * fav)<<endl ;

}









 
