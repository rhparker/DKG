% Sine-Gordon eq

function [F, J] = SGeq(x,u,D2,par)
% returns the right-hand side of our equation

%% operator

eps = par.eps;

% identity operator
N = length(x);
Id = speye(N);

% linear part
LN = eps*D2;

% linear and nonlinear part
F  = LN*u' + sin(u');
F  = F';

%% Jacobian
if nargout > 1     

% nonlinear terms
NL = -sparse(1:N,[1:N],-cos(u'),N,N); 
J = LN + NL;
 
end