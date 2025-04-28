function out = eqs_system(xs,data)
global gam_bz TOF

% universal parameters
kb = 1.380649E-23;
h = 6.62607015E-34;
R_ig = 8.3144E-3;

%extract the ethalpy and entropy from regression parameters
H_alp = xs(1);
H_beta = xs(2);

S_alp = xs(3);
S_beta = xs(4);

%extract the temperature and benzene pressures from data passed to solver
Ts = data(:,1);
PBzs = data(:,2);

%initialize the predicted TOF vector
ys = [];

%specify alpha and beta in the T-range
alphas = kb.*Ts./h.*exp(-(H_alp-Ts.*S_alp)./Ts./R_ig)
betas = exp(-(H_beta-Ts.*S_beta)./Ts./R_ig)

%loop over each experimental same and use fsolve to get coverages
for i=1:length(data(:,1))

F = @(v) v+gam_bz.*betas(i).*PBzs(i).*v.^gam_bz-1;
v_out =  fsolve(F,.1);

y = alphas(i).*PBzs(i).*v_out.^(gam_bz+1);
ys = [ys,y];
end

%report output (normalizing residuals by experiment TOF 
out = (TOF-ys')./TOF;

