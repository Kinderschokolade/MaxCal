void periodic_boundary(vector<double>& xnew, vector<double>& box) {
    for (size_t i = 0; i < xnew.size(); ++i) {
        if (xnew[i] < 0 || xnew[i] > box[i]) {
            xnew[i] -= floor(xnew[i]);
        }
    }
}

int pos(double x) {
    int y = floor(3 * x);
    return (y < 0) ? 2 : (y > 2) ? 0 : y;
}

void border(vector<double>& x, vector<double>& xnew, vector<double>& p, 
            vector<double>& Fold, vector<double>& box, 
            vector<double>& F, double& dt) {
    int pos1 = pos(x[0]);
    int pos2 = pos(xnew[0]);
    int border = (round(3 * x[0]) > 2) ? 0 : round(3 * x[0]);
    double p_loc = (x[0] - xnew[0]) / dt;
    double x0 = border / 3.0;
    double ddx = fabs(x[0] - x0);
    ddx = (ddx > 0.5 * box[0]) ? fabs(box[0] - ddx) : ddx;
    double ddt = fabs(ddx / p_loc);
    double rt = dt - ddt;
    double Fnew = Fold[0] + F[pos1] - F[pos2];
    
    if (fabs(Fold[0]) > fabs(F[pos1] - F[pos2])) {
        Fold[0] = Fnew;
        periodic_boundary(xnew, box);
        if ((pos2 - pos1) > 1) x0 = 1;
        p[0] = (p[0] > 0) 
             ? fabs(xnew[0] - x0) / rt 
             : -fabs(xnew[0] - x0) / rt;
        xnew[0] = x0 + p[0] * dt + Fnew * rt * rt;
        periodic_boundary(xnew, box);
    } else {
        p[0] = -p_loc;
        xnew[0] = x0 + p_loc * rt;
        periodic_boundary(xnew, box);
    }
}

int delta_3state(vector<double>& f, vector<double>& x, vector<double>& xo, 
                 vector<double>& box, vector<double>& U, vector<double>& T, 
                 double& gamma, double& dt, default_random_engine& generator, 
                 vector<double>& p, double& k, bool& bc, bool& err) {
    size_t dim = f.size();
    vector<double> xnew(dim, 0);
    double v_cut = box[0] / 3.0;
    uniform_real_distribution<double> lin(-0.5, 0.5);

    for (size_t i = 0; i < dim; ++i) {
        int position = pos(x[i]);
        double sigma = sqrt(24.0 * gamma * T[position] / dt);
        f[i] = sigma * lin(generator) - gamma * p[i];
        xnew[i] = x[i] + p[i] * dt + f[i] * dt * dt;
        if (pos(xnew[i]) != pos(x[i])) {
            border(x, xnew, p, f, box, U, dt);
        } else {
            p[i] = (xnew[i] - x[i]) / dt;
        }
    }

    xo[0] = x[0];
    x[0] = xnew[0];

    if (x[0] > box[0] || x[0] < 0) {
        cout << "out of bound err" << endl;
        err = true;
    }

    return 0;
} 
