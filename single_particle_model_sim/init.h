
void init_single_particle_acc(vector <double> & x, vector <double> & v, vector <double> & a, vector <double> box, int dim)
{
	x.resize(dim);
	a.resize(dim);
	v.resize(dim);

	for(int i=0 ; i < dim ; i++)
	{
		x[i] = box[i]*0.5;
		a[i] = 0.;
		v[i] = 0.;
	}
}

void init_single_particle(vector <double> & x, vector <double> & x0, vector <double> & box, vector <double> v, int dim)
{
	x.resize(dim);
	x0.resize(dim);
	v.resize(dim);

	for(int i=0 ; i < dim ; i++)
	{
		x[i] = box[i]*0.5;
		x0[i] = box[i]*0.5;
		v[i] = 0.;
	}
}
