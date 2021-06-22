double vmag(double vec[3]);
void cross(double a[3], double b[3], double *c);
double dot(double a[3], double b[3]);
double fAlt(double e, double i, double argp, double nu, double h, 
            double Re, double Rp, double mu);
double f(double e, double nu, double argp, double c0);
double fp(double e, double nu, double argp, double c0);
double fpp(double e, double nu, double argp, double c0);
int isAngBetween(double theta, double lb, double ub);
double acosr(double x);
void keplerMinMaxSphere(double r0[3], double v0[3], double rf[3], double vf[3],
                        double tf, double mu, double R, double *minAlt, 
                        double *maxAlt);
void keplerMinMax(double r0[3], double v0[3], double rf[3], double vf[3], 
                  double tf, double mu, double Re, double Rp, int Npts, 
                  double *minAlt, double *maxAlt);