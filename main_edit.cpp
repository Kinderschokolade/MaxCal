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
	vector <double> & U , vector <double> & T,int & lvl, vector <int> & lag, int &lagmin, int & lagmax,  double & gamma, string & ident, 
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
		cout<<"-pot \t 1 - potential.. better choose 1" << endl;
		cout<<" \t 2 - numeric delta pot " <<endl;
		cout<<"-ms \t number of microstates - default 3" <<endl;
		cout<<"-wf \t writing freq - default 500" <<endl;
		cout<<"-seed \t seed - default 100" << endl;
		cout<<"-lagmin \t minimum steps of tra" << endl;
		cout<<"-lagmax \t maximum steps of tra" << endl;
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
		else if (std::string(argv[i]) == "-lag") 
		{
	      		int c = i+1;
			while(std::string(argv[c])[0] != '-')
			{
				cout<<( std::string(argv[c]) ) <<endl;
				lag.push_back((int)atof(argv[c]));
				c++;
			}
		}
		else if (std::string(argv[i]) == "-lagmin") 
		{
			lagmin = (int)atof(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-lagmax") 
		{
			lagmax = (int)atof(argv[i+1]);
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

int find_ms(double & x, int & ms)
{
	return ( floor( ms * x ) );
}

int find_ms(double & x, int & ms, double & shift)
{
	int pos = (int)((x+shift)*ms*2);
	if(pos%2 == 0)
	{
		pos = pos/2;
	}
	else
	{
		pos = ms;
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


	cout<< path0 << endl;
	cout<< path1 << endl;
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



	vector <int> lag; // read dynamically
	int lagmin;
	int lagmax; // needed for trajectory analysis

	int write_freq;

	int it;
	double int_cut =0.1; // integration cut off for shift.
	double seed;
	read_flags(argc, argv, U ,T , lvl,lag, lagmin,lagmax, gamma, ident,fullT, dt,potential,ext_f, ierr,it,ms, write_freq,seed);
	int number_lags= lag.size();
	int lagm = max(lagmax,lag[number_lags-1]); // defines the max length of saved trajectory
	vector <int> bc_part (number_lags,0);

	vector <double> p(dim,0);
	vector <double> Perr(ms,0);
	vector < vector < vector <double> > > Tcount (blocknum, vector<vector <double> > 
		 (ms,vector <double> (ms,0)));

	vector <int> hist(ms,0);
	
	vector < vector <vector <int> > > Count(number_lags, 
		vector <vector <int> > (ms, vector <int>(ms,0)));

	vector < vector <int> > histo(number_lags, vector <int> (ms,0));


	double v_cut = 1./(double)lvl;

	double conf = 0.999; //How good is U approx
	k = - lvl * log(1./conf -1);

	double shift = -1./(ms*4); // 50% core position 
	//if (potential ==1 ) ierr = 1;	// not defined
	
	double fac2 = dt/(gamma*T[0]);
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

		if (potential == 1)
		{
			xc[i] = double(i) / lvl;
		}
		else
		{
			xc[i] = double(i)/lvl+findroot(beta[i-1],beta[i],U[i-1],U[i],k,k,int_cut);
		}
	}
	if (potential == 2)
	{
		xc[0] = 0.+findroot(beta[lvl-1],beta[0],U[lvl-1],U[0],k,k,int_cut);
		xc[lvl] = 1.+xc[0];
	}
	if (potential == 1)
	{
		xc[0] = 0.;
		xc[lvl] = 1.;
	}
		
	ofstream traj_file;
	string link; // use to build paths
//	cout<< it << "\t" << xc[0] << "\t" <<xc[1]<<"\t"<<xc[2]<<"\t"<<xc[3]<< endl;
	if (write_freq >0)
	{	
		link = path1+"single_particle/trajectory_"+ident+".dat";
		traj_file.open(link.c_str(),ios::out | ios::binary);
	}
	else
	{
		traj_file.open("/data/pckr194/bause/single_particle/mumpitz.dat");
	}
	ofstream param_file;
	link = path0+"single_particle/param_"+ident+".dat";
	param_file.open(link.c_str());

	link = path1+"single_particle/Sprod_"+ident+".dat";
	ofstream Sprod_file;
	Sprod_file.open(link.c_str());

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
	param_file <<"fullT \t" << fullT << endl;
	param_file <<"lagmin \t" << lagmin << endl;
	param_file <<"lagmax \t" << lagmax << endl;
	param_file <<"lag\t";
	for (int i=0 ; i<number_lags ; i++)
	{
		param_file<< "\t"<< lag[i];
	}
	param_file << endl;
	double sigmaf;
	double fav, f2av;

	x[0] = 0.52;
	xo[0] = 0.52;
	pos = find_ms(x[0],ms,shift);
	pos_old= pos;

	vector <int> survival_hist;
	vector <int> prev_pos(number_lags,pos_old);
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
	int progress = int_steps / 10.;
	int therm_steps = 100;
	double ran; // random number 
	
	vector <double> traj (lagm+1,0);
	vector <double> F_l (lagm+1,0);
	vector <double> Sprod_l (lagm+1,0);
	vector <double> pos_l (lagm+1,0);

	double Sprod=0;
	double structT=0., structTav=0.;	
	double structTav2=0. , f4av = 0., deviation =0.;
	
	force_nstate_overdamped(f, x, xo, box, U, T, gamma, dt, rng,p,k,bc, err, sigmaf, xc, ext_f,potential);
	
	vector < vector <int > > steps (ms, vector<int> (ms,0));
	vector < vector <int > > checkpos (ms, vector<int> (ms,0));
	
	for (int i=0 ; i < lagm; i++)
	{
		force_nstate_overdamped(f, x, xo, box, U, T, gamma, dt, rng,p,k,bc, err, sigmaf, xc, ext_f,potential);
		F_l.erase(F_l.begin());
		F_l.push_back(f[0]);
		pos =find_ms(x[0],ms,shift);
		pos_l.push_back(pos);
		//traj.erase(traj.begin());
		//traj.push_back(x[0]);
		pos_l.erase(pos_l.begin());
		Sprod_l.push_back((traj[lagm] - traj[lagm-1])*(F_l[lagm]+F_l[lagm-1])/(2.*T[0]) );
		Sprod_l.erase(Sprod_l.begin());
		traj.erase(traj.begin());
		traj.push_back(x[0]);
	}
	cout<< "thermed"<<endl;
//	exit(EXIT_FAILURE);
	int start,length;
	double posi;

	for (unsigned long long i=0 ; i < int_steps; i++)
	{
		force_nstate_overdamped(f, x, xo, box, U, T, gamma, dt, rng,p,k,bc, err, ran, xc, ext_f,potential);
		//file<< "pos" <<endl;
		pos =find_ms(x[0],ms,shift);
		pos_l.push_back(pos); // int to double casting is not useful
		traj.erase(traj.begin());
		traj.push_back(x[0]);
		pos_l.erase(pos_l.begin());
		F_l.erase(F_l.begin());
		F_l.push_back(f[0]);
		Sprod_l.push_back((traj[lagm-1] - traj[lagm-2])*(F_l[lagm]+F_l[lagm-1])/(2.*T[0]) );
		Sprod_l.erase(Sprod_l.begin());
		for (int u =0 ; u<ms ; u++){for(int l =0 ; l<ms ; l++){steps[u][l]+=1;}}


		if(pos != pos_old and pos != ms)
		{
			pos_old = pos;
		//	fill(steps[pos].begin(), steps[pos].end(), 1);	
		//	fill(checkpos[pos].begin(), checkpos[pos].end(), 1);	
			for(int l =0 ; l < ms ; l++)
			{
				if (checkpos[pos][l] ==0)
				{
					steps[pos][l] =1;
					checkpos[pos][l] =  1;
				}
			}
			

			checkpos[pos][pos] = 0;
			for(int l =0 ; l < ms ; l++)
			{
				length = steps[l][pos];
				if (checkpos[l][pos] == 1 and length > lagmin)
				{
					Sprod= 0;
					if (length < lagm)
					{
						for(int p =lagm - length+1 ; p<= lagm ; p++)
						{
							Sprod += Sprod_l[p];
						} // any way of saving previous calculations? 
						double lendoub = length;
					//	cout<< pos_l[lagm-length+1]<< "\t"<< pos_l[lagm-length]<< "\t"<< pos_l[lagm-length-1]  << "\t"<< l << endl;
						Sprod_file.write( (char*)&pos_l[lagm-length+1], sizeof(double));
						Sprod_file.write( (char*)&pos_l[lagm], sizeof(double));
						Sprod_file.write( (char*)&Sprod, sizeof(double));
						Sprod_file.write( (char*)&lendoub, sizeof(double));
					} // lagm and lagmax is not needed twice
					
					steps[l][pos]=0 ;
					checkpos[l][pos] = 0;
				}
			}
		}

		for (int c =0 ; c < number_lags ; c++ )
		{
			pos = pos_l[lag[c]+1];
			if (pos != ms) // only core states matter
			{

				if (pos != ms && i% lag[c] ==0) // only core sates matter
				{
					histo[c][pos]++;
					Count[c][prev_pos[c]][pos]++;
					prev_pos[c] = pos;
				}
			}

		}
		
		//file<< "av" <<endl;
		//if (i % 100000==0) cout<< i <<"\t" << int_steps <<  endl;

		//Ekin = (Ekin *i + 1./2. * p[0] * p[0]         )/(double)(i+1);
		//Eav =  (Eav    *i + 1./2. * p[0] * p[0] + U[pos])/(double)(i+1);
		//pav = (pav*i + p[0])/(double)(i+1);
		//fav = (fav*i + f[0])/(double)(i+1);
		f2av =  (f2av*i + f[0]*f[0])/(double)(i+1);
		f4av =  (f4av*i + f[0]*f[0]*f[0]*f[0])/(double)(i+1);
		structT = struct_T(x[0], k, U, xc,f[0]) ;
		structTav = (structTav *i + structT) / (double)(i+1); // is the double needed here?
		structTav2 = (structTav2 *i + structT*structT) / (double)(i+1); // is the double needed here?

		if (write_freq != 0)
		{
			if (i % write_freq == 0)
			{
				deviation = fabs(f4av - f2av *f2av) / (structTav) + f2av * fabs(structTav2 - structTav*structTav )/ (structTav *structTav) ;
				traj_file <<x[0]<<"\t"<< p[0] <<"\t"<< ran <<"\t"<< f2av/structTav <<"\t"<< sqrt(deviation/i) << endl;
			}
		}
		//file<< "endl" <<endl;
		if (i % progress == 0)
		{
			cout<< i / (double) int_steps << endl;
		}

	}

	cout<<"Sim finished"<<endl;
	
	ofstream count_file;
	link = path0+"single_particle/count_"+ident+".dat";
	count_file.open(link.c_str());
	

	count_file << "# ";	
	for (int i=0 ; i< number_lags ; i++)
	{
		count_file << lag[i] << "\t";
	}
	count_file << endl;

	for (int k=0 ; k< ms ; k++)
	{
		for (int j=0 ; j< ms ; j++)
		{
			count_file<<"C"<<k+1<<j+1<<"\t";
			for (int i=0 ; i< number_lags ; i++)
			{
				count_file << Count[i][j][k] << "\t";
			}
			count_file << endl;
		}
	}
	


	for (int k=0 ; k< ms ; k++)
	{
		count_file << "H"<<k+1 <<"\t";
		for (int i=0 ; i< number_lags ; i++)
		{
			count_file << histo[i][k] << "\t";
		}
		count_file << endl;
	}
	cout<< "count written"<< endl;
/*	
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
	
	Sprod =0;
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
	*/

}







 