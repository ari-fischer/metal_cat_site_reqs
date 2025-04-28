clc
clear

global gam_bz TOF

% import kinetic data from .csv file
% columns: {Temp (K), Bz pressure (bar), TOF (s-1)}
data_path = 'INPUT DATA PATH THERE'
data = csvread(data_path)

% exctract TOFs from data
TOF = data(:,3);

% specify integer Bz site requirements
gam_bz=6;%5.59;

% initialize the enthalpy and entropy of parameters alpha and beta
x0 = [60.3781827987836	11.0928140014420	-0.0500306109261013	0.0572592003077920]

% upper and lower bound
lb = [-100 -100 -.1 -.1]
ub = [100 100 .1 .1]

%run the regression with matlab fn lsqcurvefit
[x,resnorm,residual,exitflag,output,lambda,J] = lsqcurvefit(@eqs_system,x0,data(:,1:2),data(:,3).*0)

%calculate the error bars
ci = nlparci(x,residual,'jacobian',J)

%report outputs
outs=eqs_system(x,data(:,1:2))