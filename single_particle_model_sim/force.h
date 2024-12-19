void periodic_boundary(double & xnew, double & box, int & bc)
{
	
	bc = 0;
	if(xnew < 0) {
		xnew = box + xnew;
		bc = -1;
	}
	if(xnew > box) {
		xnew = xnew -box;
		bc = 1;
	}
}
double dH(double & x, double & k, double shift)
{
	return(-2.*exp(-2.*k*(-shift+x))*k)/
		((1.+exp(-2*k*(-shift+x)))*(1.+exp(-2*k*(-shift+x))));
}

double ddH(double & x, double & k, double shift)
{
	double e = exp(-2.*k*(x-shift));
	return (4.*k*k*e*e - 4. *k*k*e )/((1+e)*(1+e)*(1+e));
}

double H(double & x, double & k, double shift)
{
	return(1./(1.+exp(-2*k*(x-shift))));
}

double Gauss(vector <double> &x , vector <double> & sigma)
{
	return  exp(-x[0]*x[0]/(2*sigma[0]*sigma[0])-x[1]*x[1]/(2*sigma[1]*sigma[1]));
}


double Gauss(double x , double sigma)
{
	return exp(-x*x/(2*sigma*sigma));
}


//1D Sin force
double F_sin(double & x, double & k, vector <double> & U, vector <double> & xc)
{
	
	double Uup = U[0];
	double Udown = U[1];
	double lvl = U.size();

	double pos = fmod(x, 2*xc[1]);
	
	double	f = (Uup- Udown) * sin(pos * lvl * M_PI);

	return f;
}



//1D tanh force
double F_tanh(double & x, double & k, vector <double> & U, vector <double> & xc)
{
	int lvl = U.size();
	double f = (dH(x,k,xc[0]) - dH(x,k,xc[1])+ dH(x,k,xc[lvl]))* U[0] 
		+(dH(x,k,xc[lvl-1]) - dH(x,k,xc[lvl]) -dH(x,k,xc[0]))* U[lvl-1] ;

	for (int i= 1 ; i < lvl -1 ; i++)
	{
		f += (dH(x,k,xc[i]) - dH(x,k,xc[i+1]))* U[i] ;
	}
	return f;
}


//1D gaussian force
double F_gauss(double & x, double & sigma, vector <double> & U, vector<double>  & x0, double & box) 
{
	int lvl = x0.size()-1;
	double f = 0;
 	double xshift;
	double xmax = 0.7 *box;
	for (int i=0 ; i< lvl ; i++)
	{
		xshift = x - x0[i];

		if (xshift< xmax)
		{
			f+= U[i] * -xshift/(sigma*sigma) * Gauss(xshift,sigma);
		}
	}

	xshift = x - x0[lvl];
	if (xshift< xmax)
	{
		f+= U[0] * -xshift/(sigma*sigma) * Gauss(xshift,sigma);
	}

	return f;
}


double struct_T(double & x, double & k, vector <double> & U, vector <double> & xc,double & f) //1D
{
	int lvl = U.size();
	double fp = (ddH(x,k,xc[0]) - ddH(x,k,xc[1])+ ddH(x,k,xc[lvl]))* U[0] 
		+(ddH(x,k,xc[lvl-1]) - ddH(x,k,xc[lvl]) -ddH(x,k,xc[0]))* U[lvl-1] ;

	for (int i= 1 ; i < lvl -1 ; i++)
	{
		fp += (ddH(x,k,xc[i]) - ddH(x,k,xc[i+1]))* U[i] ;
	}
	
	return fp;
}

double struct_T_gauss(vector <double> & x, vector <double> & sigma, vector <double> & U, vector <vector<double> > & x0) 
{
	int lvl = x0.size();
	double fp = 0;
 	vector <double> xshift(2,0);
	double sigmax2 = sigma[0]*sigma[0];
	double sigmay2 = sigma[1]*sigma[1];
	double sigmax4 = sigmax2*sigmax2;
	double sigmay4 = sigmay2*sigmay2;

	for (int i=0 ; i< lvl ; i++)
	{
		xshift[0] = x[0]-x0[i][0];
		xshift[1] = x[1]-x0[i][1];
		fp+= U[i] * (xshift[0]*xshift[0]/sigmax4 + xshift[1]*xshift[1]/ sigmay4 - 1./sigmax2 - 1./sigmay2) * Gauss(xshift,sigma);
	}
	return fp;
}



double F_tanh(double & x, double & k, vector <double> & U, vector <double> & xc, int & dir) // 2D version 
{
	int lvl = U.size();
	double f;
	if (dir == 0)
	{
		f = (dH(x,k,xc[0]) - dH(x,k,xc[1])+ dH(x,k,xc[lvl]))* U[0] 
			+(dH(x,k,xc[lvl-1]) - dH(x,k,xc[lvl]) -dH(x,k,xc[0]))* U[lvl-1] ;

		for (int i= 1 ; i < lvl -1 ; i++)
		{
			f += (dH(x,k,xc[i]) - dH(x,k,xc[i+1]))* U[i] ;
		}
	}
	else{
		f = 0;
	}

	return f;
}

double F_gauss(vector <double> & x, vector <double> & sigma, vector <double> & U, double & Uscale, vector <vector<double> > & x0, int & dir) 
{
	int lvl = x0.size();
	double f = 0;
 	vector <double> xshift(2,0);
	for (int i=0 ; i< lvl ; i++)
	{
		xshift[0] = x[0]-x0[i][0];
		xshift[1] = x[1]-x0[i][1];
		f+= U[i] * -xshift[dir]/(sigma[dir]*sigma[dir]) * Gauss(xshift,sigma);
	}
	if (dir == 1)
	{
		f-= Uscale * (x[1] -x0[0][1]) ; // careful: the middle is defined at the y-center of the first well.. usually in the middle anyway
	}
	return f;
}

//1D tanh potential
vector <double> get_F(int steps,  double & k, vector <double> & U, vector <double> & xc, int potential)
{
	vector <double> Fprof;
	double pos;
	//if (potential == 2)
	//{
		for (int i=0 ; i < steps ; i++)
		{
			pos = double(i+ 1./2.) / double(steps);
			Fprof.push_back(F_tanh(pos,k,U,xc));
		}
	//}
	//if (potential == 1)
	//{
	//	for (int i=0 ; i < steps ; i++)
	//	{
	//		pos = double(i+1./2.) / double(steps);
	//		Fprof.push_back(F_sin(pos,k,U,xc));
	//		//cout<<pos << "\t" << Fprof[i] << endl;
	//	}
	
//	}
	return Fprof;
}

//1D Gaussian all forces
vector <double> get_Fgauss(int steps,  double & k, vector <double> & U, vector <double> & xc, double & box)
{
	vector <double> Fprof;
	double pos;
	for (int i=0 ; i < steps ; i++)
	{
		pos = double(i+ 1./2.) / double(steps);
		Fprof.push_back(F_gauss(pos,k,U,xc,box));
	}
	return Fprof;
}

//Gauss 2D
vector < vector<double> > get_U_Gauss(vector<int> steps,  vector <double>  sigma, vector <double>  U, double Uxsq, vector < vector<double> >  xm, vector <double> box)
{
	vector < vector<double> > Uvec;
	double E;
	double Exsq;
	int lvl = U.size();
	vector <double> xsh(2,0);
	vector <double> pos(2,0);

	Uvec.resize(steps[1]);
	double scale = 4.* Uxsq / (box[1]*box[1]) ; // such that Uxsq is the max of the potential at boundary
	for (int i= 0; i< steps[1] ; i++)
	{
		pos[1] = double(i+ 1./2.) / double(steps[1]) *box[1];
		Exsq = (pos[1] - box[1] / 2.) * (pos[1] - box[1] / 2. ) * scale;
		for (int z = 0 ; z< steps[0]; z++)
		{
			pos[0] = double(z+ 1./2.) / double(steps[0])* box[0];
			E =0;
			for (int j=0 ; j< lvl ; j++)
			{
				xsh[0] = pos[0] - xm[j][0];
				xsh[1] = pos[1] - xm[j][1];
				E-= U[j] * Gauss(xsh,sigma);
			}
			E+=Exsq;
			Uvec[i].push_back(E);
		}
	}
	return Uvec; 
}

//potential with repelling sides and Gaussian wells included. Linear slope included.
vector <double> get_U(int steps, vector <double> & Umax , vector <double> xbas, vector <double>  box, double sigma, double slope)
{
	vector <double> U(steps,0);
	double dx = box[0]/ (double)steps;	
	double x =0.;
	int lvl = xbas.size();
	double boundary  = xbas[0] -sigma;
	double bound_mul = 20;
	int bound_steps = boundary / dx;	

	vector <double> norm(lvl,0);
	for (int l =0 ; l<lvl ; l++)
	{
		norm[l] = Umax[l]/ Gauss(0,sigma);
	}
	for(int i= 0 ; i <steps ; i++)
	{
		U[i] = -slope * x;
		for (int l=0 ; l < lvl ; l++)
		{
			U[i] -= norm[l]*Gauss(x-xbas[l],sigma);
		}
		if (i < bound_steps)//left booundary repelling
		{
			U[i] += exp((boundary -x)*(slope*bound_mul)) -1;
		}
		if (i > steps-bound_steps )//left booundary repelling
		{
			U[i] += exp((-box[0] +boundary + x)*(slope*bound_mul)) -1;
		}
		
		x = x + dx;
	}
	return U;
}
//force for gaussian/slope/repelling boundaries
vector <double> get_F(int steps, vector <double> & Umax , vector <double> xbas, vector <double>  box, double sigma, double slope)
{
	vector <double> F(steps,0);
	double dx = box[0]/ (double)steps;	
	double x =0.;
	int lvl = xbas.size();
	double xsh;
	double boundary  = xbas[0]-sigma;
	int bound_steps = boundary / dx;	
	int bound_mul=20;

	cout<< bound_steps << endl;
	vector <double> norm(lvl,0);
	for (int l =0 ; l<lvl ; l++)
	{
		norm[l] = Umax[l]/ (Gauss(0,sigma)*sigma*sigma);
	}
	
	for(int i= 0 ; i <steps ; i++)
	{
		F[i] = +slope;
		for (int l=0 ; l < lvl; l++)
		{
			xsh = x - xbas[l];	
			F[i] -= norm[l] * xsh *  Gauss(xsh, sigma);
		}
		if (i < bound_steps)//left boondary repelling
		{
			F[i] += (slope*bound_mul)*exp((boundary -x)*(slope*bound_mul));
		}
		if (i > steps-bound_steps )//right boundary repelling
		{
			F[i] -= (slope*bound_mul)*exp((-box[0] +boundary + x)*(slope*bound_mul));
		}
	
		x = x + dx;
	}
	return F;
}

//1D tanh potential
vector < double > get_U(int steps,  double & k, vector <double> & U, vector <double> & xc, int potential)
{
	vector < double > Uvec;
	double pos,E;
	int lvl = U.size();

	for (int i=0 ; i < steps ; i++)
	{
		pos = double(i+ 1./2.) / double(steps);
		E = (H(pos,k,xc[0]) - H(pos,k,xc[1])+ H(pos,k,xc[lvl]))* U[0] 
				+(H(pos,k,xc[lvl-1]) - H(pos,k,xc[lvl]) - H(pos,k,xc[0]))* U[lvl-1] ;
	
		for (int j = 1 ; j< lvl -1 ; j++)
		{
			E += (H(pos,k,xc[j]) - H(pos,k,xc[j+1]))* U[j] ;
		}
		
		Uvec.push_back(E);

	}
	return Uvec;
}
//1D Gauss potential
vector < double > get_Ugauss(int steps,  double & k, vector <double> & U, vector <double> & x0, double & box) 
{
	vector < double > Uvec;
	double pos,E,xshift;
	int lvl = U.size();
	double xmax = 0.5 *box;

	for (int j=0 ; j < steps ; j++)
	{
		pos = double(j+ 1./2.) / double(steps);

		E=0;
		for (int i=0 ; i< lvl ; i++)
		{
			xshift = pos - x0[i];
			if (xshift< xmax)
			{
				E-= U[i] *  Gauss(xshift,k);
			}
		}
		for (int i=0 ; i< lvl ; i++) // periodic boundary conditions enforced
		{
			xshift = x0[i] + (box - pos);
			if (xshift < xmax) // ommit very small forces
			{
				E-= U[i] * Gauss(xshift,k);
			}
		}

		Uvec.push_back(E);

	}
	return Uvec;
}


//2D tanh potential
vector < vector <double> > get_U2D(int steps,  double & k, vector <double> & U, vector <double> & xc, int potential)
{
	vector < vector <double> > Uvec;
	double pos,posy,E;
	int lvl = U.size();

	Uvec.resize(steps);

	for (int i=0 ; i < steps ; i++)
	{
		posy = double(i+ 1./2.) / double(steps);
		for (int z = 0 ; z< steps; z++)
		{
			pos = double(z+ 1./2.) / double(steps);
			E = (H(pos,k,xc[0]) - H(pos,k,xc[1])+ H(pos,k,xc[lvl]))* U[0] 
				+(H(pos,k,xc[lvl-1]) - H(pos,k,xc[lvl]) - H(pos,k,xc[0]))* U[lvl-1] ;
	
			for (int j = 1 ; j< lvl -1 ; j++)
			{
				E += (H(pos,k,xc[j]) - H(pos,k,xc[j+1]))* U[j] ;
			}
			Uvec[i].push_back(E);
		}
	}
	return Uvec;
}


int find_pos(double & x,vector <double> & xc, vector <double> & box)
{
	int pos;
	int lvl = xc.size();

	if (x < xc[0]) 
	{
		pos =lvl-1;
	}
	else 
	{ 
		for (int i=1 ; i < lvl ;i++)
		{
			if (x < xc[i])
			{
				pos = i -1;
				return pos ;
			}
		}
	}
	pos =0;
	return pos;
}

double Temp(double & x , double & T1, double & T2, double & box)
{
	return T1 + (-T1+T2)*x/box;
}


int force_Tgrad(vector <double> & f, vector <double> & x, vector <double> & xo, 
	vector <double > & box, double & T1, double & T2,
	double & gamma, double & dt, mt19937 & generator,
	vector <double> & p , bool & err, double & sigmap, 
	vector <double> & Fvec, double & xmul 
	)
{
	int pos;
	int dim = x.size();
	double mean = 0;
	double sigma;	

	
	vector <double> xnew(dim,0);
	normal_distribution<double> gauss(0.0,1.0);

	pos = 0;
	for (int i=0; i < dim ; i++)
	{
		sigma = sqrt(2. * Temp(x[i],T1,T2,box[i]) * dt / gamma);
		sigmap = sigma * gauss(generator);
		pos = floor(x[i]/xmul);
		xnew[i] = x[i] + Fvec[pos] /gamma *dt + sigmap;
		p[i] = (xnew[i] - x[i]) / dt;
	}	


	xo[0] = x[0];
	x[0]  = xnew[0] ;

	return 0;


}


int force_nstate_overdamped(vector <double> & f, vector <double> & x, vector <double> & xo, 
	vector <double > & box, vector <double > & U, double & T, 
	double & gamma, double & dt, mt19937 & generator,
	vector <double> & p ,double & k, int & bc, bool & err, double & sigmap, 
	vector <double> & xc, double & ext_f,double & fparam2, double & ftype, int & potential
	)
{
	int pos;
	int dim = f.size();
	double mean = 0;
	double sigma;	

	int lvl = xc.size();


	vector <double> xnew(dim,0);

	normal_distribution<double> gauss(0.0,1.0);

	pos = 0;
	sigma = sqrt(2. * T * dt / gamma);
	for (int i=0; i < dim ; i++)
	{
		sigmap = sigma * gauss(generator);
		f[i] = F_tanh(x[i],k,U,xc,i) + ext_f;
		
		xnew[i] = x[i] + f[i]/gamma *dt + sigmap;
		p[i] = (xnew[i] - x[i]) / dt;
		periodic_boundary(xnew[i],box[i],bc);
	}	

	x[0]  = xnew[0] ;
	xo[0] = x[0] ;

	
	if (x[0] > box[0] || x[0] < 0) 
	{
		cout<< "out of bound err" <<endl;
		err = 1;
	}
	return 0;

}

int force_nstate_overdamped_floc(vector <double> & f, vector <double> & x, vector <double> & xo, 
	vector <double > & box, vector <double > & U, double & T, 
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

	double s_squ = 2*0.02*0.02;
	vector <double> xnew(dim,0);

	normal_distribution<double> gauss(0.0,1.0);

	pos = 0;
	sigma = sqrt(2. * T * dt / gamma);
	for (int i=0; i < dim ; i++)
	{
		sigmap = sigma * gauss(generator);
		
		f[i] = F_tanh(x[i],k,U,xc,i) + ext_f * exp(-(x[i]-xc[1])*(x[i]-xc[1])/s_squ );  // not generalised  
		xnew[i] = x[i] + f[i]/gamma *dt + sigmap;
		
		p[i] = (xnew[i] - x[i]) / dt;
		periodic_boundary(xnew[i],box[i],bc);
	}	

	x[0]  = xnew[0] ;
	xo[0] = x[0] ;

	
	if (x[0] > box[0] || x[0] < 0) 
	{
		cout<< "out of bound err" <<endl;
		err = 1;
	}
	return 0;

}

int force_nstate(vector <double> & f, vector <double> & x, vector <double> & xo, 
	vector <double > & box, vector <double > & U, double & T, 
	double & gamma, double & dt, mt19937 & generator,
	vector <double> & p ,double & k, int & bc, bool & err, double & sigmaf, 
	vector <double> & xc, double & ext_f
	)
{
	int pos;
	int dim = f.size();
	double mean = 0;
	double sigma;	
	
	vector <double> xnew(dim,0);

	uniform_real_distribution <double> lin(-0.5,0.5);
	
	sigma = sqrt(24. * gamma * T / dt);
	for (int i=0; i < dim ; i++)
	{
		sigmaf = sigma * lin(generator);
		f[i] = (sigmaf - gamma * p[i] + F_tanh(x[i],k,U,xc,i) + ext_f);
		if (fabs(x[i] -xo[i] ) < 0.5)
		{
			xnew[i] = 2*x[i] - xo[i] + f[i]*dt*dt; 
		}
		else if (xo[i] < x[i])
		{		
			xo[i] +=1.;
			xnew[i] = 2*x[i] - xo[i] + f[i]*dt*dt; 
		}
		else
		{
			xo[i] -=1.;
			xnew[i] = 2*x[i] - xo[i] + f[i]*dt*dt; 
		}
		p[i] = (xnew[i] - x[i]) / dt;//watch out for xo[i] can be other side of bdry
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

//gaussian force function
int force_nstate_overdamped(vector <double> & f, vector <double> & x, vector <double> & xo, 
	vector <double > & box, vector <double > & U, double & T, 
	double & gamma, double & dt, mt19937 & generator,
	vector <double> & p ,double & k, int & bc, bool & err, double & sigmap, 
	vector <double> & xc, double & fmax, double & fsig, double & fx, int & potential
	)
{
	int pos;
	int dim = f.size();
	double mean = 0;
	double sigma;	
	double frel;
	int lvl = xc.size();


	vector <double> xnew(dim,0);

	normal_distribution<double> gauss(0.0,1.0);


	sigma = sqrt(2. * T * dt / gamma);
	for (int i=0; i < dim ; i++)
	{
		sigmap = sigma * gauss(generator);
		frel = fx - x[i] ;
		if ( frel > box[i]/2.) {frel -= box[i];}
		else if(frel < -box[i]/2.) {frel += box[i];}  // boundary conditions
		f[i] = F_tanh(x[i],k,U,xc,i) + fmax * Gauss(frel,fsig )  ;
		
		xnew[i] = x[i] + f[i]/gamma *dt + sigmap;
		p[i] = (xnew[i] - x[i]) / dt;
		periodic_boundary(xnew[i],box[i],bc);
	}	

	x[0]  = xnew[0] ;
	xo[0] = x[0] ;

	
	if (x[0] > box[0] || x[0] < 0) 
	{
		cout<< "out of bound err" <<endl;
		err = 1;
	}
	return 0;
}





int force_nstate_overdamped_2D(vector <double> & f, vector <double> & x, vector <double> & xo, 
	vector <double > & box, vector <double > & U, double & Uxsq, double & T, 
	double & gamma, double & dt, mt19937 & generator,
	vector <double> & p , vector <double> & sigma, vector <int> & bc, bool & err, vector <double> & ran, 
	vector < vector <double> > & xm, vector <double> & ext_f, int & potential 
	)
{
	int pos;
	int dim = f.size();
	double mean = 0;
	double fac = sqrt(2. * T * dt / gamma);
	int lvl = xm.size();
	
	vector <double> xnew(dim,0);

	normal_distribution<double> gauss(0.0,1.0);

	for (int i=0; i < dim ; i++)
	{
		ran[i] = gauss(generator);
		f[i] = F_gauss(x,sigma,U,Uxsq,xm,i);// + ext_f[i];
		xnew[i] = x[i] + (f[i]+ext_f[i])/gamma *dt + fac*ran[i];
		p[i] = (xnew[i] - x[i]) / dt;
		periodic_boundary(xnew[i],box[i],bc[i]);
	}	
	int dir = 0;

	xo = x ;
	x  = xnew ;
	return 0;


}


