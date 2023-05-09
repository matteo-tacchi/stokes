% Approximation of the volume of a ball with p-norm bounds
% with moments and Stokes constraints
% D. Henrion, M. Tacchi, 1 Feb 22

% uses GloptiPoly for modeling the moment problem
% and MOSEK interfaced through YALMIP for solving the SDP problem

r = 3/4; % ball radius
n = 4; % dimension
p = 4; % norm (must be even)
d = 4; % relaxation degree

mset clear
mpol('xmu',n,1);
mu = meas(xmu);
mpol('xmuhat',n,1);
muhat = meas(xmuhat);
mpol('xnu',n,1);
nu = meas(xnu);
gmuhat = 1-sum(xmuhat.^p);
gmu = r^2-xmu'*xmu;
gnu = r^2-xnu'*xnu;
pows = genpow(n+1,d); pows = pows(:,2:end);
vmu = mmon(xmu,d); vmuhat = mmon(xmuhat,d); vnu = mmon(xnu,d);
y = momball(pows,p);
ME = [mom(vmu)+mom(vmuhat)==y]; % moment equations
SE = [];
for i = 1:n
 % u_i(x) = v(x), u_j(x) = 0 for j \neq i
 % Stokes equation \int div(u) \mu + \int grad(g).u \nu = 0
 SE = [SE; mom(diff(vmu,xmu(i)))+mom(diff(gnu,xnu(i))*vnu)==0];
end
mset('yalmip',true);
mset(sdpsettings('solver','mosek'))
P = msdp(max(mass(mu)),ME,SE,gmuhat>=0,gmu>=0,gnu==0);
msol(P);

disp(['bound = ' num2str(double(mass(mu)))]);
disp(['volume = ' num2str(pi^(n/2)*r^n/gamma(n/2+1))]);

