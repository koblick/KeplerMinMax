/*
 * keplerMinMaxMex.c - determine the minimum/maximum altitude of a two
 * body trajectory over an oblate spheroid.
 *
 * The calling syntax is:
 * [minAlt,maxAlt] = keplerMinMaxMex(r0,v0,rf,vf,tf,mu,Re,Rp);
 * 
 * [minAlt,maxAlt] = keplerMinMaxMex([0 0 0],[0 0 0],[0 0 0],[0 0 0],0,0,0,0);
 * To compile this program from a MATLAB terminal simply:
 *
 * mex keplerMinMaxMex.c keplerMinMax.c
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2021 Darin C. Koblick
 */

#include "stdio.h"
#include "stdlib.h"
#include "mex.h"
#include "matrix.h"
#include "keplerMinMax.h"

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double thisInitPos[3],thisInitVel[3],thisFinalPos[3],thisFinalVel[3];
  double thisAltMin,thisAltMax;
  double *r0,*v0,*rf,*vf,*tf,*minAlt,*maxAlt,mu,Re,Rp;
  int N;
  size_t mrows,ncols,i,j;
  /* create a pointer to the geodetic coordinates array */
           r0 = mxGetPr(prhs[0]);
           v0 = mxGetPr(prhs[1]);
           rf = mxGetPr(prhs[2]);
           vf = mxGetPr(prhs[3]);
           tf = mxGetPr(prhs[4]);
           mu = mxGetScalar(prhs[5]);
           Re = mxGetScalar(prhs[6]);
           Rp = mxGetScalar(prhs[7]);
  /*  get the dimensions of the r0,v0 arrays */
          mrows = mxGetM(prhs[0]);
          ncols = mxGetN(prhs[0]);   
  /*  set the output pointer to the output matrix */
        plhs[0] = mxCreateDoubleMatrix( (mwSize)mrows, (mwSize)1, mxREAL);
        plhs[1] = mxCreateDoubleMatrix( (mwSize)mrows, (mwSize)1, mxREAL);
  /*  create a C pointer to a copy of the output matrix */
           minAlt = mxGetPr(plhs[0]);  
           maxAlt = mxGetPr(plhs[1]);
  // Write this in a for-loop for increased speed!
  for (i=0; i<mrows; i++)
  {   
      /* Extract the PCI position vectors */
      thisInitPos[0] = *(r0+i);
      thisInitPos[1] = *(r0+i+mrows);
      thisInitPos[2] = *(r0+i+mrows*2);
     thisFinalPos[0] = *(rf+i);
     thisFinalPos[1] = *(rf+i+mrows);
     thisFinalPos[2] = *(rf+i+mrows*2); 
      /* Extract the PCI velocity vectors */
      thisInitVel[0] = *(v0+i);
      thisInitVel[1] = *(v0+i+mrows);
      thisInitVel[2] = *(v0+i+mrows*2);
     thisFinalVel[0] = *(vf+i);
     thisFinalVel[1] = *(vf+i+mrows);
     thisFinalVel[2] = *(vf+i+mrows*2);   
     /*  call the C subroutine */
     keplerMinMax(thisInitPos,thisInitVel,thisFinalPos,thisFinalVel, 
                  *(tf+i), mu, Re, Rp, 20, &thisAltMin, &thisAltMax);
     /* Assign Min/Max Altitudes to output arrays */
        *(minAlt+i) = thisAltMin;
        *(maxAlt+i) = thisAltMax;
  }
}