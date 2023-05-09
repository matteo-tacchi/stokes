% Approximation of the area of an l4 disk
% with SOS and without Stokes constraints
% D. Henrion, M. Tacchi, 1 Feb 22

% uses YALMIP for modeling the SOS problem
% and MOSEK for solving the SDP problem

d = 16;  % relaxation degree

x = sdpvar(2,1);
b = 1-x'*x;
g = (25/72)^4-x(1)^4-x(2)^4;
[w,wc] = polynomial(x,d);
[sb,sbc] = polynomial(x,d-2);
[sk,skc] = polynomial(x,d-2);
pows = genpow(3,d); pows = pows(:,2:end);
y = momball(pows);
ops = sdpsettings('solver','mosek');
solvesos([sos(w-sb*b),sos(sb),sos(w-1-sk*g),sos(sk)],y'*wc,ops,[wc;sbc;skc]);

double(y'*wc)

P = vectorize(sdisplay(w));
P = strrep(P,'x(1)','X1');
P = strrep(P,'x(2)','X2');
[X1,X2] = meshgrid(linspace(-1,1,1e2));
P = double(eval(P));
P(X1.^2+X2.^2>=0.95)=0;
close all
surf(X1,X2,P,'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
axis([-1 1 -1 1 0 1.5]);
view(-60,50)
axis vis3d
camlight left
colormap spring
xlabel x_1
ylabel x_2
zlabel w
