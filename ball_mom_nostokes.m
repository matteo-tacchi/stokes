% Approximation of the volume of a 3d ball
% with moments and without Stokes constraints
% D. Henrion, M. Tacchi, 1 Feb 22

% uses GloptiPolyfor modeling the moment problem
% and MOSEK for solving the SDP problem

r = 3/4; % ball radius
n = 3; % ball dimension
d = 4; % relaxation degree

mset clear
mpol('x',n,1);
mu = meas(x);
mpol('xh',n,1);
muh = meas(xh);
b = 1-xh'*xh;
g = r^2-x'*x;
pows = genpow(n+1,d); pows = pows(:,2:end);
v = mmon(x,d); vh = mmon(xh,d);
y = momball(pows);
ME = [mom(v)+mom(vh)==y];
mset('yalmip',true);
mset(sdpsettings('solver','mosek'))
P = msdp(max(mass(mu)),ME,b>=0,g>=0);
msol(P);

disp(['bound = ' num2str(double(mass(mu)))]);
disp(['volume = ' num2str(pi^(n/2)*r^n/gamma(n/2+1))]);
