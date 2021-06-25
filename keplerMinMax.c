// Purpose:
//
// Minimum/Maximum Altitude values corresponding to a two-body trajectory
// over an oblate spheroid.  This routine will also support min/max
// altitude calculations over a perfect sphere with fewer computations.
//
//
// Revision History:
// Darin C. Koblick, PhD                                     (c) 06-14-2021
// ======================= Begin Code Sequence ============================
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double vmag(double vec[3])
{
    // vector magnitude (e.g. L2 norm)
    return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

void cross(double a[3], double b[3], double *c)
{
    // cross product a x b = c
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

double dot(double a[3], double b[3])
{
    // dot product (a dot b)
    return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

double fAlt(double e, double i, double argp, double nu, double h, 
            double Re, double Rp, double mu)
{
  // Altitude function wrt nu
  return (h*h/(mu*(1+e*cos(nu))) - 
                     (Re-((Re-Rp)*pow(sin(i),2)*pow(sin(nu+argp),2))));
}
    
double f(double e, double nu, double argp, double c0)
{
    // Root function in which we're trying to solve for a zero, f(nu)
    return (e*sin(nu) + c0*sin(2*(nu+argp))*pow(e*cos(nu)+1,2));
}

double fp(double e, double nu, double argp, double c0)
{
    // First derivative of the root function df/dnu with respect to nu
    double u;
    u = nu+argp;
    return (e*cos(nu)+ 2*c0*cos(2*u)*pow(e*cos(nu)+1,2) -
           2*c0*sin(2*u)*(e*cos(nu)+1)*(e*sin(nu)));
}

double fpp(double e, double nu, double argp, double c0)
{
    // Second derivative of the root function (ddf/ddnu) with respect to nu
    double u;
      u = nu+argp;
    return (-e*sin(nu) - 4*c0*sin(2*u)*pow(e*cos(nu)+1,2) +
            4*c0*cos(2*u)*(e*cos(nu)+1)*(-e*sin(nu)) -
           (4*c0*cos(2*u)*(e*cos(nu)+1)*(e*sin(nu)) +
            2*c0*sin(2*u)*(-e*sin(nu))*(e*sin(nu)) + 
            2*c0*sin(2*u)*(e*cos(nu)+1)*(-e*cos(nu))));
}

int isAngBetween(double theta, double lb, double ub)
{
    double twopi;
    twopi = 2*acos(-1);
    // Check if an angle,theta, is between lb and ub
    if (theta < 0)
    {theta = theta + twopi;}
    else if (theta > twopi)
    {theta  = theta - twopi;}
    return (((lb <= theta) & (theta<=ub)) | 
           ((((0 <= theta) & (theta <= ub)) | (lb <= theta)) & (lb > ub)));
}

double acosr(double x)
{
    //  real part of theta = arccos(x)
    if (x > 1){x = 1;}else if (x <-1){x = -1;}
    return acos(x);
}


void keplerMinMaxSphere(double r0[3], double v0[3], double rf[3], double vf[3],
                        double tf, double mu, double R, double *minAlt, 
                        double *maxAlt)
{
// When Re = Rp, the planet is a perfect sphere and computations may be
// simplified
double e[3],h[3];
double r0Mag, v0Mag, rfMag, eMag, hMag, a, nu0, nuf, P, rp, ra, N;
double pi;
int i;
r0Mag = vmag(r0);    
cross(r0,v0,h);                                    // momentum vector
cross(v0,h,e);      
for (i=0; i<3; i++){e[i] = e[i]/mu - r0[i]/r0Mag;} // eccentricity vector
eMag = vmag(e);
// Return Solution when the orbit is circular as minAlt = maxAlt
if (eMag < 1e-8)
{
    minAlt[0] = r0Mag-R;
    maxAlt[0] = minAlt[0];
    return;
}
   pi = acos(-1);
v0Mag = vmag(v0);
rfMag = vmag(rf);
 hMag = vmag(h);
    a = -mu/(v0Mag*v0Mag-2*mu/r0Mag);
  nu0 = acosr(dot(e,r0)/(eMag*r0Mag));
  if (dot(r0,v0) < 0){nu0 = 2*pi - nu0;}
  nuf = acosr(dot(e,rf)/(eMag*rfMag));
  if (dot(rf,vf) < 0){nuf = 2*pi - nuf;}
    P = 2*pi*sqrt(a*a*a/mu);
   rp = a*(1-eMag);
  if (fabs(eMag-1)<1e-8){rp = hMag*hMag/(2*mu);} // parabolic case
   ra = a*(1+eMag);
  minAlt[0] = rp-R;
  maxAlt[0] = ra-R;
          N = tf/P;
  if ((eMag >= 1) | ((isAngBetween(pi,nu0,nuf) == 0) & (N < 1)))
  {
      maxAlt[0] = fmax(r0Mag,rfMag)-R;
  }
  if ((isAngBetween(0,nu0,nuf) == 0) & ((N < 1) | (eMag >= 1)))
  {
      minAlt[0] = fmin(r0Mag,rfMag)-R;
  }
}



void keplerMinMax(double r0[3], double v0[3], double rf[3], double vf[3], 
                  double tf, double mu, double Re, double Rp, int Npts, 
                  double *minAlt, double *maxAlt)
{
    // Purpose:
    //
    // Main Routine to ingest the inital state vector (r0,v0), final
    // state vector (rf,vf), flight time, tf, and central
    // body parameters (mu,Re,Rp).  It will then return the altitude
    // extrema of the orbit segment.
    //
    // Inputs:
    // --------
    // r0                       [1 x 3]             Position Vector in
    //                                              Cartesian Coordiantes
    //                                              corresponding to t=0
    //
    //  v0                      [1 x 3]             Velocity Vector in
    //                                              Cartesian Coordinates
    //                                              corresponding to t=0
    //
    //  rf                      [1 x 3]             Position Vector in
    //                                              Cartesian Coordinates
    //                                              corresponding to t=tf 
    //
    //  vf                      [1 x 3]             Velocity Vector in
    //                                              Cartesian Coordinates
    //                                              corresponding to t=tf 
    //
    //  tf                      double              Duration of Segment
    //
    //  mu                      double              Standard Gravitational
    //                                              Parameter of Central 
    //                                              Body
    //
    //  Re                      double              Equatorial Radius of
    //                                              the Central Body
    //
    //  Rp                      double              Polar Radius of the 
    //                                              Central Body
    //
    //  Npts                    integer             Number of Iniital
    //                                              Points to use for the
    //                                              true anomaly grid
    //
    // Outputs:
    // ---------
    //
    // minAlt                   double              Minimum Surface 
    //                                              Altitude of the orbit
    //                                              segment
    //
    // maxAlt                   double              Maximum Surface 
    //                                              Altitude of the orbit 
    //                                              segment
    //
    // Revision History:
    // Darin C. Koblick, PhD                                (c) 06-15-2021
    // --------------------- Begin Code Sequence --------------------------
    double r0Mag, rfMag, eMag, hMag, v0Mag, nMag;
    double h[3],e[3],kCrossH[3],k[3]={0,0,1},n[3];
    double a, P, nu0, inc, nuf, kxhMag, argp, c0, pi;
    double f_nu, fp_nu, fpp_nu, nu_pg, nu, thisAlt;
    int i, j, maxIter = 20;
    // Perform a quick check on the equatorial/polar radius, if they are
    // equal, then call the spherical implmentation of kepler min/max
    if (fabs(Re-Rp) < 1e-8)
    {
        keplerMinMaxSphere(r0,v0,rf,vf,tf,mu,Re,minAlt,maxAlt);
        return;
    }
    // Otherwise, continue with the oblate spheroid implmentation
           pi = acos(-1);
        r0Mag = vmag(r0);
        rfMag = vmag(rf);
        v0Mag = vmag(v0);
        cross(r0,v0,h);                             // Momentum vector
         hMag = vmag(h);                            // |Momentum|
            a = -mu/(v0Mag*v0Mag - 2*mu/r0Mag);     // SMA
            P = 2*pi*sqrt(a*a*a/mu);                // Period
         cross(v0,h,e);                             // Eccentricity vector
         for (i=0; i<3; i++){e[i] = e[i]/mu - r0[i]/r0Mag;}
         eMag = vmag(e);                            // |Eccentricty|
          nu0 = acosr(dot(e,r0)/(eMag*r0Mag));      
         if (dot(r0,v0) < 0){nu0 = 2*pi - nu0;}     // True Anomaly @ t=0
          nuf = acosr(dot(e,rf)/(eMag*rfMag));      
         if  (dot(rf,vf) < 0){nuf = 2*pi-nuf;}      // True Anomaly @ t=tf
          inc = acosr(h[2]/hMag);                   // Inclination [rad]
         cross(k,h,kCrossH);
         kxhMag = vmag(kCrossH);
         argp = acosr(dot(kCrossH,e)/(kxhMag*eMag));
         if (e[2] < 0){argp = 2*pi - argp;}         // Argument of Perigee
         // Recompute the argument of perigee for equatorial orbits as
         // it's undefined when k x h = [0 0 0]
         if (inc < 1E-8)
         {
          // See https://en.wikipedia.org/wiki/Argument_of_periapsis
             argp = atan2(e[1],e[0]);
             if (h[2] < 0){argp = 2*pi-argp;}
         }
         // For Circular inclined orbits, compute the argument of latitude,
         // set argp = 0, and compute the effective true anomaly, 
         // u = argp + nu see Eq. 2-89 Vallado
         if ((eMag < 1E-8) & (inc >= 1E-8))
         {
           argp = 0;
             cross(k,h,n); // n = k x h
           nMag = vmag(n); 
            nu0 = acosr(dot(n,r0)/(nMag*r0Mag));
            if (r0[2] < 0){nu0 = 2*pi - nu0;} // nu0 = u0-argp = u0
            nuf = acosr(dot(n,rf)/(nMag*rfMag));
            if (rf[2] < 0){nuf = 2*pi - nuf;} // nuf = uf-argp = uf 
         }
         else if (eMag < 1E-8)
         {
         // For Circular Equatorial Orbits, compute the true longitude
         // set RAAN = 0, compute the argp, and back out the effective
         // true anomaly, lambda_true = argp + RAAN + nu, Eq. 2-92 Vallado
            cross(k,h,n);             // n = k x h 
            nu0 = acosr(r0[0]/r0Mag);
            if (r0[1] < 0){nu0 = 2*pi-nu0;}
            nuf = acosr(rf[0]/rfMag);
            if (rf[1] < 0){nuf = 2*pi-nuf;}
            nu0 = nu0-argp;
            nuf = nuf-argp;
         }
         // Ensure that nu0 and nu0 are betwen 0 and 2*pi
         if (nu0 < 0){nu0 = nu0 + 2*pi;}
         if (nuf < 0){nuf = nuf + 2*pi;}
         // Constant for Root function:
                   c0 = (mu*(Re-Rp)*sin(inc)*sin(inc))/(hMag*hMag);
         // Determine the min/max altitudes of the orbit segements:  
            minAlt[0] = fAlt(eMag,inc,argp,nu0,hMag,Re,Rp,mu);
            maxAlt[0] = fAlt(eMag,inc,argp,nuf,hMag,Re,Rp,mu);
         if (maxAlt[0] < minAlt[0])
         {    thisAlt = minAlt[0];
            minAlt[0] = maxAlt[0];
            maxAlt[0] = thisAlt;
         }
         for (i=0; i<Npts; i++)
         {
             // Initial guess [0: 2*pi/Npts : 2*pi-2*pi/Npts]
                   nu = i*(2*pi/Npts);
             // Halley's method for each initial guess:
             for (j=0; j<maxIter; j++)
             {
                nu_pg = nu;
                 f_nu = f(eMag,nu,argp,c0);
                fp_nu = fp(eMag,nu,argp,c0);
               fpp_nu = fpp(eMag,nu,argp,c0);
                   nu = nu-2*f_nu*fp_nu/(2*fp_nu*fp_nu-f_nu*fpp_nu);
              if (fabs(nu-nu_pg) < 1e-12){break;}
             }
             // Check if true anomaly is contained within the segment:
             if ((isAngBetween(nu,nu0,nuf) == 1) | (P <= tf+1e-8))
             {
                  thisAlt = fAlt(eMag,inc,argp,nu,hMag,Re,Rp,mu);
                minAlt[0] = fmin(minAlt[0],thisAlt);
                maxAlt[0] = fmax(maxAlt[0],thisAlt);
             } 
         }
}
