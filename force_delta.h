
void periodic_boundary(vector <double> & xnew, vector <double> & box)
{
	int dim = xnew.size();
	for (int i =0 ; i < dim ; i++)
 	{	
		if(xnew[i] < 0) {
			xnew[i] = xnew[i] - floor(xnew[i]) ;// wrong for changing box
		}
		if(xnew[i] > box[i]) {
			xnew[i] = xnew[i] - floor(xnew[i]) ;
		}
	}
}

int pos (double x)
{
	int y = floor(3*x);
	if(y < 0) y = 2;
	if(y > 2) y =0;
	return y;
}

void border(vector <double>  & x, vector < double> & xnew, vector <double> & p, 
	vector <double> & Fold, vector <double> & box,
	vector <double> & F, double & dt
	)
{
	int pos1 = pos(x[0]);
	int pos2 = pos(xnew[0]);
	int border = round(3 * x[0]);
	int direction;
	if (border > 2) border = 0;
	double p_loc = (x[0] - xnew[0])/ dt;
	double x0 = border * 1./3.;
	double ddx = fabs(x[0]-x0);
	if (ddx > 0.5 *box[0]) ddx=fabs(box[0]-ddx); // not elegant
	double ddt = fabs(ddx/(p_loc));
	double rt = dt -ddt;
	//if(rt < 0.01 * dt) rt = 0.1 * rt; //because causes problems..
	double Fnew = Fold[0] + F[pos1] - F[pos2];
	//cout<< xnew[0]  << endl;
	if (fabs(Fold[0]) > (-F[pos1]+F[pos2]) )
	{
		Fold[0] =Fnew;	
		periodic_boundary(xnew,box);
		if ((pos2 -pos1)>1) // 0 -> 2
			x0 = 1;
		if (p[0] > 0) {
			p[0] = fabs(xnew[0] - x0) / rt;
			xnew[0] = x0 + p[0] *rt + Fnew * rt *rt ;
		}
		else {
			p[0] = -fabs(xnew[0] - x0 )/ rt;
			xnew[0] = x0 + p[0] *dt + Fnew * rt *rt ;
		}
		periodic_boundary(xnew, box);
	}
	else 
	{
		p[0] = -p_loc;
		xnew[0] = x0 + p_loc*rt;
		//cout<<"mirror " <<  xnew[0] << "\t" << p[0] <<endl;
		periodic_boundary(xnew,box);
	}
	
	
}


int delta_3state(vector <double> & f, vector <double> & x, vector <double> & xo, 
	vector <double > & box, vector <double > & U, vector <double> & T, 
	double & gamma, double & dt, default_random_engine & generator,
	vector <double> & p ,double & k, bool & bc, bool & err
	)
{
	int dim = f.size();
	double mean = 0;
	double sigma;	
	int position;
	
	vector <double> xnew(dim,0);
	double v_cut = box[0]/3. ; //only in x direction

	//normal_distribution < double > gauss(mean, sigma);
	uniform_real_distribution <double> lin(-0.5,0.5);
	

	for (int i=0; i < dim ; i++)
	{
		position = pos(x[i]);
		sigma = sqrt(24. * gamma * T[position] / dt);
		f[i] = sigma * lin(generator) - gamma * p[i];
		xnew[i] = x[i] + p[i]*dt + f[i] * dt *dt;
		if (pos(xnew[i]) != pos(x[i])) 
			border(x, xnew, p, f, box, U, dt);
		else p[i] = (xnew[i] - x[i]) / dt; 
		//periodic_boundary(xnew,box,bc);
		//cout<< x[i] << "\t" << xnew[i] << "\t" << endl;
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
