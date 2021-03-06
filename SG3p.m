% clear all;
format long;

% number of sites;
n=60;
center=n/2;

% initial coupling; typically at AC-limit eps=0; 
eps=0;

% field initialization
uk=ones(1,n) * -pi;

% kink
uk(n/2:end)= pi; 
% % site-centered
% uk(n/2) = 0;

% % kink-antikink
% u1 = ones(1,n) * -pi;
% Loffset = 4;
% Roffset = 4;
% u1(n/2-Loffset:n/2+Roffset-1) = pi;

% kink-antikink-kink
width = 10;
u1 = ones(1,n) * -pi;
u1(n/2:end) = pi;
u1(n/2-width+1:n/2) = pi;
u1(n/2+1:n/2+width) = -pi;

% plot(u1, '.', 'MarkerSize',40);

% iteration   
it=2;

% lattice index
x=linspace(1,n,n)-n/2;

% paramaters for continuation
delta_eps = 0.05;
end_eps = 0.25;
% threshold for Newton's Method
threshold = 1e-08;

% difference operators
D1 = fddiffeasy(n, 1, 1, 'none');
D2 = fddiffeasy(n, 2, 1, 'none');
D2(1,2) = 2;
D2(n,n-1) = 2;

% continuation in coupling epsilon

while ( eps <= end_eps )
    par.eps = eps;
    
    % kink
    u=zeros(1,n);
    while ( norm(u-uk) > threshold )
        u = uk;
        % evaluation of function and computation of Jacobian
        % DNLS eq
        [F,J] = SGeq(x,uk,D2,par);
        % Newton correction step
        cor = ( J \ F' )'; 
        uk = u - cor;
        % convergence indicator: should converge quadratically
        % disp( norm(cor) );
    end;
    
    % kink-antikink
    u=zeros(1,n);
    while ( norm(u-u1) > threshold )
        u = u1;
        % evaluation of function and computation of Jacobian
        % DNLS eq
        [F,J] = SGeq(x,u1,D2,par);
        % Newton correction step
        cor = ( J \ F' )'; 
        u1 = u - cor;
        % convergence indicator: should converge quadratically
        % disp( norm(cor) );
    end;
        
    % eigenvalues d and eigenvectors v of stability matrix 
    [F,J] = SGeq(x,u1,D2,par);
    N = length(D2);
    A2 = eye(N);
    A1 = 0 * A2;
    % EVP is A0 + lambda A1 + lambda^2 A2
    % [V,lambda] = polyeig(A0,A1,A2);
    [V, lambda]  = quadeig(A2, A1, -J);
    [V2,lambda2] = eig(full(J));
    lambda2 = diag(lambda2);
    
    [F,Jk] = SGeq(x,uk,D2,par);
    [Vk, lambdak] = quadeig(A2, A1, -Jk);
    
    l0 = lambdak(1);
    w0 = l0^2;
    v0 = Vk(:,1);
    M = v0'*v0;
    offset = width/2;
    a1 = v0(center+offset)^2 - v0(center+offset-1)^2;
    deltaw = eps*sqrt(2)*a1/M;
    wpred = [ w0 - deltaw ; w0; w0 + deltaw ];
    wact  = lambda2(1:3);
    
    
    % store solution and stability
    u_store(:,it)=u1;
      
    % visualize the continuation profiles and stability on the fly 
    subplot(3,2,1)
    plot(x,u,'-o','LineWidth',1)
    drawnow;
    
    subplot(3,2,2)
    plot(lambda,'.','MarkerSize',20)
    axis([-2 2 -2 2]);
    drawnow;
    
    subplot(3,2,3)
    plot(x,V(:,1),'-o','LineWidth',1)
%     axis([-10,10,-0.1, 0.7]);
    drawnow;
    
    subplot(3,2,4)
    plot(x,V(:,3),'-o','LineWidth',1)
    drawnow;
    
%     subplot(3,2,5)
%     plot(x,V(:,141),'-o','LineWidth',1)
%     drawnow;
%     
%     subplot(3,2,6)
%     plot(x,V(:,143),'-o','LineWidth',1)
%     drawnow;
    
    target = 0.67;
    threshold = 0.2;
    inteigs = lambda( find( abs( imag(lambda) - target) < threshold ) );
    
%     
%     sm = find(abs(eigenv) < 1e-6 );
%     big = find(abs(eigenv) > 1e-6 );
    
%     % eigenvalues near 0
%     subplot(3,2,3)
%     plot(x,vi(:,sm));
%     axis([x(1) x(end) -1 1]);
%     drawnow;
%     
%     subplot(3,2,4)
%     plot(x,wi(:,sm));
%     axis([x(1) x(end) -1 1]);
%     drawnow;
    
%     % other eigenvalues
%     if (length(big) > 0)
%         subplot(3,2,5)
%         plot(x,vi(:,big));
%         axis([x(1) x(end) -1 1]);
%         drawnow;
% 
%         subplot(3,2,6)
%         plot(x,wi(:,big));
%         axis([x(1) x(end) -1 1]);
%         drawnow;
%     end
    
    
    % increment indices and epsilon
    it=it+1;
    eps=eps+delta_eps;
   
end;

% disp(inteigs);


