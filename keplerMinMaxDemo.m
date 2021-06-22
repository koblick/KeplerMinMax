function keplerMinMaxDemo()
%% Purpose:
%
%  MATLAB Demo Routine Demonstrating the accuracy and speed of the
%  keplerMinMax C implementation as reported in the 2021 AMOS conference
%  paper titled "Novel Closed Form Solution for Orbit Segment Altitude 
%  Extrema Over Spherical and Oblate Central Bodies"
%
%% Revision History:
%  Darin C. Koblick                                         (c) 06-22-2021
%% ------------------------- Begin Code Sequence --------------------------

%% Oblate Spheroid (WGS-84) Constants:
 params.mu = 398600.4418;
 params.Re = 6378.137;  
 params.Rp = params.Re*(1-1/298.257223563);

 
%% Load pre-genterated interception data:
 data = load('GEO2LEOIntData.mat');
   r0 = data.r0;
   v0 = data.V1;
   rf = data.rf;
   vf = data.V2;
   tf = data.tf;

%% Compute the Extrema w/ Novel Method:
           timer = tic;
 [minAlt,maxAlt] = keplerMinMaxAdapter(r0,v0,rf,vf,tf,params,2);
 fprintf(1,'%s\n',[mfilename,':: keplerMinMax Execution Time = ', ...
                  num2str(toc(timer)),' [s]']);
              
%% Compute the Extrema w/ fminbnd:  
minAltN = NaN(size(tf));
maxAltN = NaN(size(tf));
timer2 = tic;
for te=1:size(r0,1)   
    [minAltN(te,1),maxAltN(te,1)] = keplerfMin(r0(te,:),v0(te,:),tf(te),params);
end
fprintf(1,'%s\n',[mfilename,':: fminbnd Execution Time = ', ...
                  num2str(toc(timer2)),' [s]']);
       
%% Show Relative Errors Between Optimization Methods:              
figure('color',[1 1 1]);
minAltErr = minAltN-minAlt;
maxAltErr = maxAltN-maxAlt;
semilogy(abs(minAltErr)); hold on;
semilogy(abs(maxAltErr),'r');
ylabel('Relative Error [km]');
xlabel('Transfer Trajectory');
legend('alt_{min}','alt_{max}');
grid on;
              
%% Show Trajectories:              
figure('color',[1 1 1]);
[x,y,z] = sphere(30);
surf(x.*params.Re,y.*params.Re,z.*params.Re,'FaceColor','b','EdgeColor','k'); hold on;
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
for tt=1:size(r0,1)
    [~,rv] = ode113(@twoBody,[0 tf(tt)],[r0(tt,:),v0(tt,:)],opts,params.mu);
    plot3(rv(:,1),rv(:,2),rv(:,3),'-','color',[1 0 0 0.2]);
end
axis equal;
axis off;

end

function etaDot = twoBody(t,eta,mu)
%Two-body force model:
etaDot = NaN(size(eta));
etaDot(1:3) = eta(4:6);                         %velocity
etaDot(4:6) = -mu.*eta(1:3)./norm(eta(1:3)).^3; %acceleration
end