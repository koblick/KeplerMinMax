function [minAlt,maxAlt] = keplerMinMaxAdapter(r0,v0,rf,vf,tf,params,dim3)
%% Purpose:
%
%  This routine will compute the minimum and maximum altitude of a two-body
%  trajectory over an oblate spheroid as described in the 2021 AMOS
%  conference paper titled "Novel Closed Form Solution for Orbit Segment
%  Altitude Extrema Over Spherical and Oblate Central Bodies".
%
%% Inputs:
%
%  r0                           [N x 3]             Position Vector in
%                                                   Cartesian Coordiantes
%                                                   corresponding to t=0
%
%  v0                           [N x 3]             Velocity Vector in
%                                                   Cartesian Coordinates
%                                                   corresponding to t=0
%
%  rf                           [N x 3]             Position Vector in
%                                                   Cartesian Coordinates
%                                                   corresponding to t=tf 
%
%  vf                           [N x 3]             Velocity Vector in
%                                                   Cartesian Coordinates
%                                                   corresponding to t=tf 
%
%  tf                           [N x 1]             Duration of Segment
%
%
%  params                       struct             Structure Containing
%                                                  Properties of the
%                                                  Primary Body/Planet
%                                                  Re = Equatorial radius
%                                                  Rp = Polar Radius
%                                                  mu = Standard Grav Param
%
% dim3                          integer            Singleton Specifier for
%                                                  [x,y,z] dimension of 
%                                                  the cartesian vectors
%
%% Outputs:
%
% minAlt                       [N x 1]             Minimum Altitude of the
%                                                  satellite over the
%                                                  planet surface at any
%                                                  point in the orbit
%                                                  segment
%
% maxAlt                       [N x 1]             Maximum Altitude of the
%                                                  satellite over the
%                                                  planet surface at any
%                                                  point in the orbit
%                                                  segment
%
%% Revision History:
%  Darin C. Koblick                                         (c) 06-15-2021
%% ------------------------- Begin Code Sequence --------------------------

%% Flatten the Inputs:
       r0 = fDim(r0,dim3);
       v0 = fDim(v0,dim3);
       rf = fDim(rf,dim3);
       vf = fDim(vf,dim3);
[tf,fSeq] = fDim(tf,dim3);

%% Ensure the vector lengths are all uniform:
if size(r0,1) ~= size(v0,1) || ...
   size(r0,1) ~= size(rf,1) ||  ...
   size(rf,1) ~= size(vf,1) || ...
   size(vf,1) ~= size(tf,1)
   error([mfilename,':: Input Vector Dimensions Not Uniform']); 
end

%% Call the C MATLAB Executable:
[minAlt,maxAlt] = keplerMinMaxMex(r0,v0,rf,vf,tf,params.mu,params.Re,params.Rp);

%% Reshape the Outputs:
minAlt = eDim(minAlt,fSeq);
maxAlt = eDim(maxAlt,fSeq);
end