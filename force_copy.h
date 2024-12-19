void periodic_boundary(double & xnew, double & box, int & bc)
{
	
	bc = 0;
	if(xnew < 0) {
		xnew = xnew - floor(xnew) ;// wrong for changing box
		bc = -1;
	}
	if(xnew > box) {
		xnew = xnew - floor(xnew) ;
		bc = 1;
	}
}
double dH(double & x, double & k, double shift)
{
	return(-2.*exp(-2.*k*(-shift+x))*k)/
		((1.+exp(-2*k*(-shift+x)))*(1.+exp(-2*k*(-shift+x))));
}

double H(double & x, double & k, double shift)
{
	return(1./(1.+exp(-2*k*(x-shift))));
}


double F_sin(double & x, double & k, vector <double> & U, vector <double> & xc)
{
	
	double Uup = U[0];
	double Udown = U[1];
	double lvl = U.size();

	double pos = fmod(x, 2*xc[1]);
	
	double	f = (Uup- Udown) * sin(pos * lvl * M_PI);

	return f;
}




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

vector <double> get_F(int steps,  double & k, vector <double> & U, vector <double> & xc, int potential)
{
	vector <double> Fprof;
	double pos;
	if (potential == 2)
	{
		for (int i=0 ; i < steps ; i++)
		{
			pos = double(i+ 1./2.) / double(steps);
			Fprof.push_back(F_tanh(pos,k,U,xc));
			//cout<<pos << "\t" << Fprof[i] << endl;
		}
	}
	if (potential == 1)
	{
		for (int i=0 ; i < steps ; i++)
		{
			pos = double(i+1./2.) / double(steps);
			Fprof.push_back(F_sin(pos,k,U,xc));
			//cout<<pos << "\t" << Fprof[i] << endl;
		}
	
	}
	return Fprof;
}

vector <double> get_U(int steps,  double & k, vector <double> & U, vector <double> & xc, int potential)
{
	vector <double> Uvec;
	double pos,E;
	int lvl = U.size();

	if (potential == 2)
	{
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
	}
	if (potential == 1)
	{
		for (int i=0 ; i < steps ; i++)
		{
			double x = double(i+ 1./2.) / double(steps);
			pos = fmod(x, 2. * xc[1]);
			Uvec.push_back((U[0]- U[1]) * cos(pos * lvl * M_PI));
		}	
	}
	return Uvec;
}



int find_pos(double & x,vector <double> & xc)
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

/*
int force_3state(vector <double> & f, vector <double> & x, vector <double> & xo, 
	vector <double > & box, vector <double > & U, vector <double> & T, 
	double & gamma, double & dt, default_random_engine & generator,
	vector <double> & p ,double & k, int & bc, bool & err, double & sigmaf, 
	vector <double> & xc, double & ext_f
	)
{
	int pos;
	int dim = f.size();
	double mean = 0;
	double sigma;	
	
	vector <double> xnew(dim,0);

	//normal_distribution < double > gauss(mean, sigma);
	uniform_real_distribution <double> lin(-0.5,0.5);
	
	for (int i=0; i < dim ; i++)
	{
		pos = find_pos(x[i],xc);
		//pos = floor(3*x[i]);
		sigma = sqrt(24. * gamma * T[pos] / dt);
		sigmaf = sigma * lin(generator);
		f[i] = (sigmaf - gamma * p[i] + F(x[i],k,U,xc) + ext_f);
		//xnew[i] = 2*x[i] - xo[i] + f[i]*dt*dt; 
		xnew[i] = x[i] + p[i]*dt + f[i] * dt *dt; //improve to verlet!!
		p[i] = (xnew[i] - x[i]) / dt;//watch out for xo[i] can be other side of bdry
		periodic_boundary(xnew,box,bc);
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

int force_nstate(vector <double> & f, vector <double> & x, vector <double> & xo, 
	vector <double > & box, vector <double > & U, vector <double> & T, 
	double & gamma, double & dt, default_random_engine & generator,
	vector <double> & p ,double & k, int & bc, bool & err, double & sigmaf, 
	vector <double> & xc, double & ext_f
	)
{
	int pos;
	int dim = f.size();
	double mean = 0;
	double sigma;	
	
	vector <double> xnew(dim,0);

	//normal_distribution < double > gauss(mean, sigma);
	uniform_real_distribution <double> lin(-0.5,0.5);
	
	for (int i=0; i < dim ; i++)
	{
		pos = find_pos(x[i],xc);
		//pos = floor(3*x[i]);
		sigma = sqrt(24. * gamma * T[pos] / dt);
		sigmaf = sigma * lin(generator);
		f[i] = (sigmaf - gamma * p[i] + F(x[i],k,U,xc) + ext_f);
		//xnew[i] = 2*x[i] - xo[i] + f[i]*dt*dt; 
		xnew[i] = x[i] + p[i]*dt + f[i] * dt *dt; //improve to verlet!!
		p[i] = (xnew[i] - x[i]) / dt;//watch out for xo[i] can be other side of bdry
		periodic_boundary(xnew,box,bc);
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
*/
/*
int force_nstate_overdamped(vector <double> & f, vector <double> & x, vector <double> & xo, 
	vector <double > & box, vector <double > & U, vector <double> & T, 
	double & gamma, double & dt, default_random_engine & generator,
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
		xnew[i] = x[i] + (F(x[i],k,U,xc) + ext_f)/gamma *dt + sigmap;
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



int force_nstate_overdamped(vector <double> & f, vector <double> & x, vector <double> & xo,
        vector <double > & box, vector <double > & U, vector <double> & T,
        double & gamma, double & dt, generator & rng,
        vector <double> & p ,double & k, int & bc, bool & err, double & sigmap,
        vector <double> & xc, double & ext_f, int & potential
        )
{
        int pos;
        int dim = f.size();
        double mean = 0;
        double sigma;

        int lvl = xc.size();

//      const auto F& = F_tanh; //stanni!

//      if (potential == 1)
//      {
//              const auto& F =F_sin;
//      }
        vector <double> xnew(dim,0);

        //normal_distribution < double > gauss(mean, sigma);
        //uniform_real_distribution <double> lin(-0.5,0.5);
        boost::random::normal_distribution<double> gauss(0.0,1.0);

        for (int i=0; i < dim ; i++)
        {
                pos = find_pos(x[i],xc);
                //pos = floor(3*x[i]);
                sigma = sqrt(2. * T[pos] * dt / gamma);
                sigmap = sigma * gauss(rng);
                //xnew[i] = 2*x[i] - xo[i] + f[i]*dt*dt; 
                xnew[i] = x[i] + (F_tanh(x[i],k,U,xc) + ext_f)/gamma *dt + sigmap;
                p[i] = (xnew[i] - x[i]) / dt;
                periodic_boundary(xnew[i],box[i],bc);
        }

        xo[0] = x[0] ;
        x[0]  = xnew[0];

        if (x[0] > box[0] || x[0] < 0)
        {
                cout<< "out of bound err" <<endl;
                err = 1;
        }
        return 0;
}

*/

