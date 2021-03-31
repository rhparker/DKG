clear; format long;

% number of sites;
n=100;

% initial coupling; typically at AC-limit eps=0; 
eps=0;

par.p = 1
par.c = 1

% convenience
c = par.c;
ht = 2*c;

% field initialization
u1=zeros(1,n);
u1(n/2)=ht;
u1(n/2-5)=ht;
% u1(n/2+5)=ht;
% u1(n/2+4)=ht;

% iteration   
it=1;

% lattice index
x=linspace(1,n,n)-n/2;

% paramaters for continuation
delta_eps = 0.02;
% end_eps = 0.05;
end_eps = 0.5;

% threshold for Newton's Method
threshold = 1e-10;

% difference operators
D1 = fddiffeasy(n, 1, 1, 'none');
D2 = fddiffeasy(n, 2, 1, 'none');

% continuation in coupling epsilon

while ( eps < end_eps )
    u=zeros(1,n);
    par.eps = eps;
    while ( norm(u-u1) > threshold )
        u = u1;
        % evaluation of function and computation of Jacobian
        [F,J] = DKG(x,u1,D2,par);
        % Newton correction step
        cor = ( J \ F' )'; 
        u1 = u - cor;
        % convergence indicator: should converge quadratically
    end
    disp(norm(cor));
    
    % construction of linear stability matrix for eigenvalue problem
    [F,J] = DKdV3_int(x,u1,D2,par);
    
    % eigenvalues d and eigenvectors v of stability matrix 
    [v,d]=eig(full(J));
    d1=diag(d);
    
    % store solution and stability
    u_store(:,it)=u1;
    d_store(:,it)=d1;
    e_store(it)=eps;

    % visualize the continuation profiles and stability on the fly 
    subplot(2,1,1)
    plot(x,u,'-o')
    drawnow;
    subplot(2,1,2)
    plot(d1, zeros(length(d1)),'o')
    axis([-3 10 -1 1]);
    drawnow;
    
    % increment indices and epsilon
    it=it+1;
    eps=eps+delta_eps;
    
end;