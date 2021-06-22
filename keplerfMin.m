function [minAlt,maxAlt,alt] = keplerfMin(r0,v0,tf,params)
%% Purpose:
%
%  This routine will compute the min/max altitude of a two-body trajectory
%  using fminbnd optimization. This will be much slower than a pure
%  analytic solution to the problem, but helps with verification and
%  validation of runtime performance and accuracy.
%
%% Inputs:
%
%  r0                           [1 x 3]             Position Vector in
%                                                   Cartesian Coordiantes
%                                                   corresponding to t=0
%
%  v0                           [1 x 3]             Velocity Vector in
%                                                   Cartesian Coordinates
%                                                   corresponding to t=0
%
%  rf                           [1 x 3]             Position Vector in
%                                                   Cartesian Coordinates
%                                                   corresponding to t=tf 
%
%  params                       struct             Structure Containing
%                                                  Properties of the
%                                                  Primary Body/Planet
%                                                  Re = Equatorial radius
%                                                  Rp = Polar Radius
%                                                  mu = Standard Grav Param
%
%% Outputs:
%
% minAlt                       [1 x 1]             Minimum Altitude of the
%                                                  satellite over the
%                                                  planet surface at any
%                                                  point in the orbit
%                                                  segment
%
% maxAlt                       [1 x 1]             Maximum Altitude of the
%                                                  satellite over the
%                                                  planet surface at any
%                                                  point in the orbit
%                                                  segment
%
%% Revision History:
%  Darin C. Koblick, PhD
%% --------------------------- Begin Code Sequence ------------------------
if nargin == 0
       params.mu = 398600.4418;
       params.Re = 6378.137;
       params.Rp = 6356.75231424518;
              r0 = [-10052.2019157853         -21361.4927313261                       0];
              v0 = [4.9063933729238           3.1136594128655                         0];
              tf = 60*60*2;
 [minAlt,maxAlt,alt] = keplerfMin(r0,v0,tf,params);
 figure('color',[1 1 1]);
 plot(alt); hold on;
 plot(1:numel(alt),repmat(minAlt,size(alt)),'r');
 plot(1:numel(alt),repmat(maxAlt,size(alt)),'b');
 grid on;
 ylabel('Surface Altitude [km]');
 return;
end
 stepSize = 5;
%% Step One: Propagate the trajectory from t=0 to t=tf:
     opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
   [t,rv] = ode113(@twoBody,[0 tf],[r0,v0],opts,params.mu);

%% Step Two: Compute the surface altitude:
      alt = rv2alt(rv,params);

%% Step Three: Guess the min/max altitude solutions:
               t_i = linspace(0,tf,ceil(tf)*10)';
             alt_i = interp1(t,alt,t_i,'spline');
  [minAltG,idxMin] = min(alt_i);
  [maxAltG,idxMax] = max(alt_i);
              tMin = t_i(idxMin); 
              tMax = t_i(idxMax);

              
%% Step Four: Numerically refine guesses within bounds of the system:
             fOpts = optimset('TolX',1e-12);
              opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
           tmin_lb = max([0,tMin-stepSize]);
           tmin_ub = min([tf,tMin+stepSize]);
%Check that minimum altitude does not occur at the boundary:
if tMin > 0 && tMin < tf
        [~,minAlt] = fminbnd(@(t)propAlt(t,r0,v0,params,opts),tmin_lb,tmin_ub,fOpts); 
else
           minAlt = minAltG;
end
           tmax_lb = max([0,tMax-stepSize]);
           tmax_ub = min([tf,tMax+stepSize]);
%Check that the maximum altitude does not occur at the boundary:
if tMax > 0 && tMax < tf
        [~,maxAlt] = fminbnd(@(t)-propAlt(t,r0,v0,params,opts),tmax_lb,tmax_ub,fOpts);
            maxAlt = -maxAlt;
else
          maxAlt = maxAltG;
end
end

function alt = propAlt(tGuess,r0,v0,params,opts)
    [~,rv] = ode113(@twoBody,[0 tGuess],[r0,v0],opts,params.mu);
       alt = rv2alt(rv(end,:),params);
end

function etaDot = twoBody(t,eta,mu)
%Two-body force model:
etaDot = NaN(size(eta));
etaDot(1:3) = eta(4:6);                         %velocity
etaDot(4:6) = -mu.*eta(1:3)./norm(eta(1:3)).^3; %acceleration
end

function alt = rv2alt(rv,params)
% [r,v] ==> surface altitude
  rMag = sqrt(sum(rv(:,1:3).^2,2));
  hVec = cross(rv(:,1:3),rv(:,4:6),2);
  hMag = sqrt(sum(hVec.^2,2));
   inc = real(acos(hVec(:,3)./hMag));
%Argument of Latitude, u:
  nVec = cross(repmat([0 0 1],[size(hVec,1) 1]),hVec,2);
  nMag = sqrt(sum(nVec.^2,2));
     u = real(acos(dot(nVec,rv(:,1:3),2)./(nMag.*rMag)));
   idx = rv(:,3) < 0;
u(idx) = 2*pi - u(idx); 
%Geocentric latitude:
   phi = asin(sin(inc).*sin(u));
%For equatorial case, assume phi = 0;
phi(inc == 0) = 0;
%Planet Surface Radius:
  R_sl = 0.5.*((params.Re + params.Rp) + (params.Re-params.Rp).*cos(2.*phi));
   alt = rMag-R_sl;
end

