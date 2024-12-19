#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream> 
#include <string> 
#include <unistd.h>
#include <random> 
#include "mpi.h"
using namespace std;

//#include "init.h"
#include "find_shift.hpp"
//#include "force_delta.h"  
#include "force.h"

void ReadNumbers( const string & s, vector <double> & v) {
        istringstream is( s );
        string stuff;
        double n;
        is >> n;
        v.push_back(n);
        while(!is.eof())
        {
                is >> n;
                if(is.fail()) break;
                v.push_back(n);
        }
}


void import_matrix(string filename_X, vector < vector < double > >& v)
{
        ifstream file_X;
        string line;
        file_X.open(filename_X.c_str());
        int rows =1 ;
        if (file_X.is_open())
        {
                getline(file_X, line);
                v.resize(1);
                ReadNumbers( line, v[0]);
                while (getline(file_X, line))
                {
                  	v.resize(rows + 1);
                  	ReadNumbers( line, v[rows]);
                    	rows++;
                }
                file_X.close();
        }
        else{
                cout << "file open failed"<<endl;
        }
    
}


//function to read flags of call
void read_flags (int argc, char* argv[0], 
	vector <double> & U , double & T,int & lvl, vector <int> & lag, double & gamma, vector<double> & sigma, string & ident, 
	double & fullT, double & dt, int & potential, vector <double> & f,
	bool & ierr,int & it, vector <int> & MS, int & write_freq, double & seed, double & Ux)
{ 
	ostringstream arg;

	gamma = 1;
	ident = "x";
	fullT = 1;
	dt = 0.01;
	potential = 1;
	ierr =0;
	write_freq = 500;
	seed = 100;
	sigma[0] = 0.1;
	sigma[1] = 0.1;
	MS[0] = 3;
	MS[1] = 3;
	f[0] = 0.;
	f[1] = 0.;
	Ux = 0.;

	if(	std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")
	{
		cout<<"-lvl \t number of states" << endl;
		cout<<"-Ui \t Potential at pos i - default 0" << endl;	
		cout<<"-Ux \t strength x^2 Potential in y direction - default 0" << endl;	
		cout<<"-T \t Temperature - default 1" << endl;	
		cout<<"-sigma[1,2] \t deviation of gauss potential" << endl;	
		cout<<"-gamma \t fricition coeff - default 1" << endl;
		cout<<"-o \t recognition" << endl;
		cout<<"-fullT \t simulation time " << endl;
		cout<<"-dT \t time steps" << endl;
		cout<<"-pot \t 1 - deltapot" << endl;
		cout<<" \t 2 - numeric delta pot " <<endl;
		cout<<"-msi \t number of microstates in direction i - default 3" <<endl;
		cout<<"-fi \t exteral force in direction i - default 0" <<endl;
		cout<<"-wf \t writing freq - default 500" <<endl;
		cout<<"-seed \t seed - default 100" << endl;
		ierr =1;
	} 
	for(int i=0; i<argc; i++) 
	{
		if (std::string(argv[i]) == "-lvl")  lvl = (int) atof(argv[i+1]);
	}
	U.resize(lvl,0);
	T=1;
	
	for(int i=0; i<argc; i++) 
	{
		if ( argv[i][1]== 'U')
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
			}
				
			if(std::string(argv[i]) == "-Ux") 
			{	
	      			Ux  = (double)atof(argv[i+1]);
			}
		}
		else if (std::string(argv[i]) == "-lag") 
		{
	      		int c = i+1;
			while(std::string(argv[c])[0] != '-')
			{
				lag.push_back((int)atof(argv[c]));
				c++;
			}
		}
		else if(std::string(argv[i]) == "-sigma1") 
		{
	      		sigma[0]  = (double)atof(argv[i+1]);
		}
		else if(std::string(argv[i]) == "-sigma2") 
		{
	      		sigma[1]  = (double)atof(argv[i+1]);
		}	
		else if (std::string(argv[i]) == "-T") 
		{
	      		T  = (double)atof(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-gamma") 
		{
	      		gamma  = (double)atof(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-wf")
		{
	      		write_freq  = (int)atof(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-MS1") 
		{
	      		MS[0]  = (double)atof(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-MS2") 
		{
	      		MS[1]  = (double)atof(argv[i+1]);
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
		else if(std::string(argv[i]) == "-f1")
		{
			f[0] = (double) atof(argv[i+1]);
		}
		else if(std::string(argv[i]) == "-f2")
		{
			f[1] = (double) atof(argv[i+1]);
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


// highly inefficient, check for more efficient alg
int find_ms( vector <double> & x, vector< vector <double> > & vor, vector <double> & box, int & ms) 
{
	int pos;
	double d = box[0] * box[0] + box[1] * box[1];
	double dnew,p,q;
	for(int i=0; i< ms ; i++)
	{
		p = min( abs(x[0] - vor[i][0]), min(x[0],vor[i][0]) +box[0] - max(x[0],vor[i][0]));
		q = min( abs(x[1] - vor[i][1]), min(x[1],vor[i][1]) +box[1] - max(x[1],vor[i][1]));
		dnew = p*p + q*q;
		if (dnew < d)
		{
			pos = i;
			d = dnew;
		}
	}
	return pos;
}


int find_ms(vector <double> & x, vector <int> & MS, vector <double> & box)
{
	return ( floor( MS[1] * x[1]/box[1] ) * MS[0] + floor(MS[0] * x[0]/box[0] ) );
}



int main(int argc, char *argv[])
{
	MPI::Init(argc, argv);
       	const int np= MPI::COMM_WORLD.Get_size (  );
        const int root = 0;
        int rank=MPI::COMM_WORLD.Get_rank ( );

	char hostname[HOST_NAME_MAX];
	char username[LOGIN_NAME_MAX];
	char keyname[] = "mbause";
	string path0, path1;
	gethostname(hostname, HOST_NAME_MAX);
	getlogin_r(username, LOGIN_NAME_MAX);
	path0 = "/data/isilon/bause/";
	path1 = "/data/pckr194/bause/";
	
	if (rank == root)
	{
		cout<< "np "<< np << endl;
	}
	const int dim = 2; // controls number of dims 
	int lvl ; // # of T,U regions
	vector < int > MS(dim,0);
	int ms; // # microstates
	double gamma = 5;
	double dt;
	double fullT;
	double int_steps;
	double Ekin=0, pav=0, Ekin_mom;
	int blocksize ;

	int number_steps = 10;

	vector <double> U;	
	double T;
	
	vector <double> box(dim,1);
	box[0] = 3;
	box[1] = 1;
	vector <double> v, x, acc;
	vector <double> f(dim,0);
	
	double E = 0;
	double Eav;
	int pos;
	bool err = false;
	vector <int> bc(dim,0);
	bool ierr;
	string ident;
	int block;
	int blocknum = np;

	int pos_old;
	int potential;
	vector <double> ext_f(dim,0);
	double padd=0.;

//	int countfac = 100; 
	vector<int>  lag;

	int write_freq;
	double Uxsq;
	int it;
	double seed;
	vector <double> sigma(2,0);

	read_flags(argc, argv, U,T,lvl, lag, gamma,sigma, ident,fullT, dt,potential,ext_f, ierr,it,MS, write_freq,seed,Uxsq);
	// think of a better way to define potential form
	ms = MS[0] * MS[1];
	int number_lags= lag.size();
	int lagm = lag[number_lags-1];
	
	blocksize = fullT/dt  ;

	vector <double> p(dim,0);
	vector <double> Perr(ms,0);
//	vector < vector <vector <int> > > Counnumber_lags, vector< vector <int> > (ms,vector< int > (ms, 0));
	int* Count = new int[number_lags*ms*ms];
//	vector < vector <vector <int> > > Tcount(blocknum, vector< vector <int> > (ms,vector< int > (ms, 0)));
//	vector <vector <int> > block_hist (blocknum, vector <int>(ms,0));
//	vector < vector <int> > histo(number_lags, vector <int> (ms,0));
	int* histo = new int[number_lags*ms];
//	vector < vector <double > > transition (ms, vector<double> (ms,0));
//	vector < vector <double > > tr_err (ms, vector<double> (ms,0));
//	vector < vector <double > > W (ms, vector<double> (ms,0));
//	vector < vector <double > > Werr(ms, vector<double>(ms,0));
//	vector < vector <double > > S (ms, vector<double> (ms,0));
//	vector < vector <double > > A (ms, vector<double> (ms,0));
	for (int i=0 ; i < number_lags*ms*ms ; i++)
	{
		Count[i] =0;
	}
	for (int i=0 ; i < number_lags*ms ; i++)
	{
		histo[i] =0;
	}



	cout<< "alloc"<< endl;

	double pot_mid = box[0]/(double)lvl;

	vector <double> xc(lvl+1,0);

	int_steps = fullT / (dt); // /np
	
	mt19937 rng(seed*(rank+1)); // seed each rank differently	

	double beta = 1. / T; 	
	for (int i=1 ; i< lvl ; i++)
	{

		if (potential == 1)
		{
			xc[i] = box[0] *  double(i) / lvl;	
		}
		
	}
	if (potential == 1)
	{
		xc[0] = 0.;
		xc[lvl] = 3.;
	}
	else
	{
		xc[lvl] = 1.+xc[0];
	}

	// Center of potential fixed here
	vector < vector <double> > xm(lvl, vector <double> (2,0.5));
	xm[0][0] =  0.5;
	xm[1][0] =  1.5;
	xm[2][0] =  2.5;

	string link; // use to build paths
	ofstream traj_file;

	if (write_freq !=0)
	{
		link = path1+"single_particle/trajectory_"+ident+"_"+to_string(rank)+".dat";
		traj_file.open(link.c_str());
	}
	else
	{
		traj_file.open("/data/pckr194/bause/single_particle/mumpitz.dat");
	}

	vector <MPI_File> Slist(number_lags);
	for (int i=0 ; i < number_lags ; i++)
	{
		link = path1+"single_particle/Sprod_"+ident+"_"+to_string(lag[i])+".dat";
		MPI_File_open(MPI_COMM_WORLD, link.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &Slist[i]);
	}

	if (rank == root){
	cout<< "write param"<< endl;
	ofstream param_file;
	link = path0+"single_particle/param_"+ident+".dat";

	param_file.open(link.c_str());
	for (int i=0 ; i<lvl ; i++)
	{
		param_file <<"U"<<i<< "\t\t" << U[i] << endl;
	}
	param_file <<"Uxsq"<< "\t\t" << Uxsq << endl;
	param_file <<"T"<< "\t\t" << T << endl;
	for (int i=0 ; i<lvl+1 ; i++)
	{
		param_file <<"xc"<<i<< "\t\t" << xc[i] << endl;
	}

	param_file <<"fullT\t\t" << fullT <<endl;
	param_file <<"gamma\t\t" << gamma <<endl;
	param_file <<"potential\t" << potential << endl;
	param_file <<"sigma1\t \t" << sigma[0] << endl;
	param_file <<"sigma2\t \t" << sigma[1] << endl;
	for (int i=0 ; i < dim ; i++)
	{
		param_file <<"extf"<< i <<"\t\t" << ext_f[i] << endl;
	}
	param_file <<"dT \t\t" << dt << endl;
	param_file <<"microstates \t" << MS[0] << " x " << MS[1]  << endl;
	param_file <<"seed \t\t" << seed << endl;
	param_file <<"fullT \t\t" << fullT  << endl;
	} // end rank ==0
//	vector < vector <double> > voronoi;
//	string filename = path0+"single_particle/cluster_8000.dat";
//	import_matrix(filename,voronoi);
//	int m= voronoi.size();
//	int mm = voronoi[0].size();
//	if (rank ==0){
//	for (int i=0 ; i< m ; i++)
//	{
//		for (int j=0 ; j< mm ; j++)
//		{
//			cout<< voronoi[i][j]<< " \t";
//		}
//		cout<< endl;
//	}}

	vector <double> ran(dim,0);
	double fav, f2av;
	
	x.resize(2);
	acc.resize(2);
	x[0] = 0.0;
	x[1] = 0.5;
	acc[0] = 0.;
	acc[1] = 0.;

	pos = find_ms(x,MS,box);
	pos_old= pos;

	vector <int> prev_pos(number_lags,pos_old);

	ofstream pot_file;
	ofstream potms_file;
	if (rank == root)
	{
		cout<< "write potential" << endl;
		link = path0+"single_particle/potential_"+ident+".dat";
		pot_file.open(link.c_str());
	
		link = path0+"single_particle/potential_ms_"+ident+".dat";
		potms_file.open(link.c_str());
		vector<int> steps(2,100);
		vector < vector < double> > Uvec = get_U_Gauss(steps,sigma,  U, Uxsq, xm,box);
	
		vector < vector <double> > Uvec_ms = get_U_Gauss(MS,sigma,U, Uxsq, xm,box);
	
		for (int i=0 ; i < steps[1] ; i++)
		{
			for (int j =0 ; j< steps[0] ; j++)
			{
				pot_file <<i/100. * box[1] << "\t" << j / 100. * box[0] << "\t"<< Uvec[i][j]  <<endl;
			}
		}
		for (int i=0 ; i < MS[1] ; i++)
		{
			for (int j =0 ; j < MS[0] ; j++)
			{
				potms_file <<i/(double)MS[1] *box[1] << "\t"<< j/(double)MS[0] * box[0] << "\t" << Uvec_ms[i][j]  <<endl;
			}
		}
	
	
		potms_file.close();
	
	} // end rank ==0
	vector<double> trajy(lagm,0);
	vector<double> trajx(lagm,0);
	vector<double> Fx_l(lagm,0);
	vector<double> Fy_l(lagm,0);
	vector<double> pos_l(lagm,0);
	vector<double> SprodSx_l(lagm,0);
	vector<double> SprodSy_l(lagm,0);

	int therm_steps = 1000;	
	double SprodSx = 0;
	double SprodSy =0;
		
	double fac2 = dt/(gamma*T);
	for (int i=0 ; i < therm_steps; i++)
	{
		force_nstate_overdamped_2D(f, x, acc, box, U,Uxsq, T, 
			gamma, dt, rng,p,sigma,bc, err,ran, xm, ext_f,potential);
		pos =find_ms(x,MS,box);
		pos_l.push_back(pos); // int to double casting is not useful
		trajx.erase(trajx.begin());
		trajx.push_back(x[0]);
		trajy.erase(trajy.begin());
		trajy.push_back(x[1]);
		pos_l.erase(pos_l.begin());
		Fx_l.erase(Fx_l.begin());
		Fx_l.push_back(f[0]);
		Fy_l.erase(Fy_l.begin());
		Fy_l.push_back(f[1]);	
		SprodSx_l.push_back((trajx[lagm-1] - trajx[lagm-2])*(Fx_l[lagm]+Fx_l[lagm-1])/(2.*T) );
		SprodSx_l.erase(SprodSx_l.begin());
		SprodSy_l.push_back((trajy[lagm-1] - trajy[lagm-2])*(Fy_l[lagm]+Fy_l[lagm-1])/(2.*T) );
		SprodSy_l.erase(SprodSy_l.begin());
	}
	if (rank == 0)
	{
		cout<<  "start"<<endl;
	}
	MPI_Request request;
	double writev[4];
	int kk = int_steps /10;
	int ms2 = ms*ms;
	int start;
	for (unsigned long long i=0 ; i < int_steps; i++)
	{
		if ( i % kk == 0)
		{
			cout<< i << endl;
		}
		force_nstate_overdamped_2D(f, x, acc, box, U, Uxsq, T, 
			gamma, dt, rng,p,sigma,bc, err, ran, xm, ext_f,potential);
		pos =find_ms(x,MS,box);
		pos_l.push_back(pos); // int to double casting is not useful
		trajx.erase(trajx.begin());
		trajx.push_back(x[0]);
		trajy.erase(trajy.begin());
		trajy.push_back(x[1]);
		pos_l.erase(pos_l.begin());
		Fx_l.erase(Fx_l.begin());
		Fx_l.push_back(f[0]);
		Fy_l.erase(Fy_l.begin());
		Fy_l.push_back(f[1]);	
		SprodSx_l.push_back((trajx[lagm-1] - trajx[lagm-2])*(Fx_l[lagm]+Fx_l[lagm-1])/(2.*T) );
		SprodSx_l.erase(SprodSx_l.begin());
		SprodSy_l.push_back((trajy[lagm-1] - trajy[lagm-2])*(Fy_l[lagm]+Fy_l[lagm-1])/(2.*T) );
		SprodSy_l.erase(SprodSy_l.begin());
		SprodSx =0;
		SprodSy =0;

		for (int c =0 ; c < number_lags ; c++ )
		{
			start = c >0 ? lag[c-1]+2 : 2;
			for(int p = start ; p < lag[c]+2 ; p++) // speed up by making use of previous integration
			{
				SprodSx += SprodSx_l[p];
				SprodSy += SprodSy_l[p];
			}
			writev[0] = pos_l[1];
			writev[1] = pos_l[lag[c]+1];
			writev[2] = SprodSx;
			writev[3] = SprodSy;
			//MPI_File_iwrite_shared(Slist[c], writev, 5, MPI_DOUBLE, &request);	
			MPI_File_write_shared(Slist[c], writev, 4, MPI_DOUBLE, MPI_STATUS_IGNORE);
			if (i % lag[c] == 0)
			{
				writev[0] = prev_pos[c];
				histo[c*ms+pos]++;
				Count[c*ms2+prev_pos[c]*ms+pos]++; 
				prev_pos[c] = pos;
			}

		}
		if (write_freq != 0)
		{
			if (i % write_freq == 0)
			{	
				Ekin_mom = 1./2. * (p[0] * p[0] + p[1] * p[1]);
				Ekin = (Ekin *i + Ekin_mom)/(double)(i+1);
				Eav =  (Eav    *i + Ekin_mom + U[pos])/(double)(i+1);
				//pav = (pav*i + p[0])/(double)(i+1);
				//fav = (fav*i + f[0])/(double)(i+1);
				f2av =  (f2av*i + f[0]*f[0]+ f[1]*f[1])/(double)(i+1);

				if(bc[0] != 0 || bc[1] != 0) {traj_file << endl;}
				traj_file <<i<< "\t"<<x[0]<<"\t"<< x[1]<< "\t" << ran[0] << "\t"<< ran[1] << "\t"<< f[0] << "\t"<< f[1]<< "\t" << Ekin << "\t"<< Eav << endl;
			}
		}
	}
	int* Count_sum;

	int* histo_sum;
	if (rank ==0)
	{
		Count_sum = new int[number_lags*ms*ms];
		histo_sum = new int[number_lags*ms];
	}
//  collecting only works for smaller matrices

	MPI::COMM_WORLD.Reduce(Count,Count_sum, number_lags*ms*ms, MPI_INT,MPI_SUM, root);
	MPI::COMM_WORLD.Reduce(histo,histo_sum, number_lags*ms, MPI_INT,MPI_SUM, root);

        stringstream ss;
        ss << rank;
        string rank_str = ss.str();
        ofstream count_file;
        link = path0+"single_particle/count_"+ident+"_"+rank_str+".dat";

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
                                        count_file << Count[i*ms*ms+j*ms+k] << "\t";
                                }
                                count_file << endl;
                }
        }

        for (int k=0 ; k< ms ; k++)
        {
                count_file << "H"<<k+1 <<"\t";
                for (int i=0 ; i< number_lags ; i++)
                {
                        count_file << histo[i*ms+k] << "\t";
                }
                count_file << endl;
        }


        MPI::COMM_WORLD.Barrier();

	if(rank ==root)
	{
	
		Count = Count_sum;
		histo = histo_sum;
		
		cout<<"full " << Count[0] << endl;

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
					count_file << Count[i*ms*ms+j*pos+k] << "\t";
				}
				count_file << endl;
			}
		}
		
		for (int k=0 ; k< ms ; k++)
		{
			count_file << "H"<<k+1 <<"\t";
			for (int i=0 ; i< number_lags ; i++)
			{
				count_file << histo[i*ms+k] << "\t";
			}
			count_file << endl;
		}
	}
	for (int i=0 ; i < number_lags; i++)
	{
		MPI_File_close(&Slist[i]);
	}

/*	
does not work here, dont know why	
		ofstream res_file;
		link = path0+"single_particle/res_"+ident+".dat";
		res_file.open(link.c_str());
		
		double tr, Wtr, Pval;
		vector <double> Pstat(ms,0);
	
		cout<< blocknum << "\t"<< block_hist.size() <<"\t"<< Tcount.size()  << endl;	

		cout<< Pstat.size() << "\t"<< Perr.size() << "\t"<< ms << endl;
		cout<< W.size()<< "\t" << W[0].size() << "\t"<< ms << endl;
		cout<< tr_err.size()<< "\t" << tr_err[0].size() << "\t"<< ms << endl;
		cout<< transition.size()<< "\t" << transition[0].size() << "\t"<< ms << endl;
		cout<< Werr.size()<< "\t" << Werr[0].size() << "\t"<< ms << endl;
		cout<< block_hist[0].size()  << "\t"<< ms << endl;

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
*/
//	} // end of root part

	MPI::Finalize();
}









 