% Approximation of the area of a planar disk
% with SOS and with Stokes constraints
% D. Henrion, M. Tacchi, 1 Feb 22

% uses YALMIP for modeling the SOS problem
% and MOSEK for solving the SDP problem

d = 16;  % relaxation degree

x = sdpvar(2,1);
b = 1-x'*x;
g = (1/16-(x(1)-1/2)^2-x(2)^2)*((x(1)+1/2)^2+x(2)^2-1/16);
[w,wc] = polynomial(x,d);
[sb,sbc] = polynomial(x,d-2);
[sk,skc] = polynomial(x,d-4);
[pk,pkc] = polynomial(x,d-4);
[v1,vc1] = polynomial(x,d-3);
[v2,vc2] = polynomial(x,d-3);
pows = genpow(3,d); pows = pows(:,2:end);
y = momball(pows);
ops = sdpsettings('solver','mosek');
solvesos([sos(w-sb*b),sos(sb),sos(w-1-jacobian(v1,x(1))-jacobian(v2,x(2))-sk*g),...
    sos(sk),sos(-v1*jacobian(g,x(1))-v2*jacobian(g,x(2))-pk*g)],y'*wc,ops,[wc;sbc;skc;pkc;vc1;vc2]);

double(y'*wc)

P = vectorize(sdisplay(w));
P = strrep(P,'x(1)','X1');
P = strrep(P,'x(2)','X2');
[X1,X2] = meshgrid(linspace(-1,1,1e2));
P = double(eval(P));
P(X1.^2+X2.^2>=0.95)=0;
close all
surf(X1,X2,P,'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
axis([-1 1 -1 1 0 2.2]);
view(-60,50)
axis vis3d
camlight left
colormap spring
xlabel x_1
ylabel x_2
zlabel w
