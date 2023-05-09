function y = momball(a,p)
% Returns the moment of monomial x^a on the unit p-norm ball
% y = \int_{x : \sum_i |x|^p_i <= 1} x^a dx
%
% For p=2 cf. Theorem 3.1 in J. B. Lasserre, E. S. Zeron. Solving a class
% of multivariate integration problems via Laplace techniques.
% Applicationes Mathematicae 28(4):391-405, 2001.
% Note that there is an incorrect factor 2^{âˆ’n} in the right handside
% of equation (3.3) in this reference.
%
% For other values of p, see Lemma A.3 in M. Tacchi, J. B. Lasserre, D. Henrion. 
% Stokes, Gibbs and volume computation of semi-algebraic sets. HAL 02947268, 2022.
%
% D. Henrion, M. Tacchi, 1 Feb 22

[m,n] = size(a);
y = zeros(m,1);
if nargin < 2, p=2; end;

for k = 1:m

 if all(~rem(a(k,:),2))
  y(k) = (2/p)^n*prod(gamma((a(k,:)+1)/p))/ gamma(1+(n+sum(a(k,:)))/p);
 end

end
