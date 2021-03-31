%% 1p

load intersited05;

x = [-n/2:n/2-1];

figure('DefaultAxesFontSize',20);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on;
plot(x,uk,'ob','MarkerSize',10, 'LineWidth',1.5)
plot(x,uk,'-b','LineWidth',1.5)
xlabel('$n$', 'Interpreter','latex');
ylabel('$u_n$', 'Interpreter','latex');

figure('DefaultAxesFontSize',20);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(lambdak, '.','MarkerSize',30);
axis( [-1e-10,1e-10,-1.6,1.6] );
xlabel('Re $\lambda$', 'Interpreter','latex');
ylabel('Im $\lambda$', 'Interpreter','latex');

figure('DefaultAxesFontSize',20);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on;
plot(x,Vk(:,1),'ob','MarkerSize',10, 'LineWidth',1.5)
plot(x,Vk(:,1),'-b','LineWidth',1.5)
xlabel('$n$', 'Interpreter','latex');
ylabel('$v_n$', 'Interpreter','latex');

figure('DefaultAxesFontSize',20);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on;
plot(x,Vk(:,41),'ob','MarkerSize',10, 'LineWidth',1.5)
plot(x,Vk(:,41),'-b','LineWidth',1.5)
xlabel('$n$', 'Interpreter','latex');
ylabel('$v_n$', 'Interpreter','latex');

%% 2p

load 2kink;

x = [-n/2:n/2-1];

figure('DefaultAxesFontSize',20);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on;
plot(x,u1,'ob','MarkerSize',10, 'LineWidth',1.5)
plot(x,u1,'-b','LineWidth',1.5)
xlabel('$n$', 'Interpreter','latex');
ylabel('$u_n$', 'Interpreter','latex');
% axis([-20,20,-4,4]);

%% 3p

load 3kink;

x = [-n/2:n/2-1];

figure('DefaultAxesFontSize',20);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on;
plot(x,u1,'ob','MarkerSize',10, 'LineWidth',1.5)
plot(x,u1,'-b','LineWidth',1.5)
xlabel('$n$', 'Interpreter','latex');
ylabel('$u_n$', 'Interpreter','latex');
axis([-20,20,-4,4]);

%% plot of turning point in d vs N

N = [2 3 4 5 6 7 8 9 10];
dvals = [0.599465 0.911334 1.20514 1.48752 1.76241 2.03199 2.29754 2.55983 2.81910];

N = N(3:end);
dvals = dvals(3:end);

p = polyfit(N,dvals,1);

figure('DefaultAxesFontSize',20);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on;
plot(N, dvals, '.', 'MarkerSize', 30);
plot(N,polyval(p,N), 'LineWidth', 2);
xlabel('$N$', 'Interpreter','latex');
ylabel('$d_0$', 'Interpreter','latex');

%% kink-antikink

load kak05_8;
x = [-n/2:n/2-1];

figure('DefaultAxesFontSize',20);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(lambda, '.','MarkerSize',30);
axis( [-1e-10,1e-10,-2,2] );
xlabel('Re $\lambda$', 'Interpreter','latex');
ylabel('Im $\lambda$', 'Interpreter','latex');
ax1=gca;
% inset
ax2 = axes('Position',[.7 .7 .2 .2])
box on;
plot(lambda, '.','MarkerSize',30);
axis( [-1, 1 ,0.57,0.573] );

% % inset
% ax3 = axes('Position',[.7 .7 .7 .7])
% box on;
% plot(lambda, '.','MarkerSize',30);
% axis( [-1, 1 ,0.57,0.573] );

figure('DefaultAxesFontSize',20);
hold on;
plot(x,V(:,1),'ob','MarkerSize',12, 'LineWidth',1.5)
plot(x,V(:,3),'.r','MarkerSize',30, 'LineWidth',1.5)
plot(x,V(:,1),'-b','LineWidth',1.5);
plot(x,V(:,3),'--r','LineWidth',1.5)
xlabel('$n$', 'Interpreter','latex');
ylabel('$v_n$', 'Interpreter','latex');
legend({ '$\lambda = 0.5713i$',  '$\lambda = 0.5718i$'}, 'Interpreter','latex');


figure('DefaultAxesFontSize',20);
hold on;
plot(x,V(:,45),'ob','MarkerSize',10, 'LineWidth',1.5)
plot(x,V(:,47),'.r','MarkerSize',35, 'LineWidth',1.5)
plot(x,V(:,45),'-b','LineWidth',1.5);
plot(x,V(:,47),'--r','LineWidth',1.5)
xlabel('$n$', 'Interpreter','latex');
ylabel('$v_n$', 'Interpreter','latex');
legend({ '$\lambda = 0.9908i$',  '$\lambda = 0.9999i$'}, 'Interpreter','latex', 'location','northwest');


%% eigenvalue accuracy plots, 2p

load goldstoned025;

cutoff = 4;
d = 0.25;

l1p = imag( sqrt( wp(1,1:cutoff) ) );
l1a = imag( sqrt( wa(1,1:cutoff) ) );
err1 = abs( (l1p-l1a)./l1a );

l2p = imag( sqrt( wp(2,1:cutoff) ) );
l2a = imag( sqrt( wa(2,1:cutoff) ) );
err2 = abs( (l2p-l2a)./l2a );

x = NN(1:cutoff);
r = (1 + 2*d + sqrt(1 + 4*d) ) / (2*d);

figure('DefaultAxesFontSize',20);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on;
plot( x, log(err1), '.b', x, log(err2), '.r', 'MarkerSize', 30);

p1 = polyfit(x,log(err1),1);
p2 = polyfit(x,log(err2),1);
plot( x, polyval(p1,x), '-b', x, polyval(p2,x), '--r', 'LineWidth',2);
xlabel('$N$', 'Interpreter','latex');
ylabel('log relative error', 'Interpreter','latex');

%% eigenvalue accuracy plots, 3p

load goldstone3p;

cutoff = 4;

l1p = imag( sqrt( wp(1,1:cutoff) ) );
l1a = imag( sqrt( wa(1,1:cutoff) ) );
err1 = abs( (l1p-l1a)./l1a );

l2p = imag( sqrt( wp(2,1:cutoff) ) );
l2a = imag( sqrt( wa(2,1:cutoff) ) );
err2 = abs( (l2p-l2a)./l2a );

l3p = imag( sqrt( wp(3,1:cutoff) ) );
l3a = imag( sqrt( wa(3,1:cutoff) ) );
err3 = abs( (l3p-l3a)./l3a );

x = NN(1:cutoff);
r = (1 + 2*d + sqrt(1 + 4*d) ) / (2*d);
r0 = (1 + 2*d + w0 + sqrt( (1+w0)*(w0 + 4*d) ) ) / (2*d);

figure('DefaultAxesFontSize',20);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on;
plot( x, log(err1), '.b', x, log(err2), '.r', x, log(err3), '.k', 'MarkerSize', 30 );

p1 = polyfit(x,log(err1),1);
p2 = polyfit(x,log(err2),1);
p3 = polyfit(x,log(err3),1);
plot( x, polyval(p1,x), '-b', x, polyval(p2,x), '--r', x, polyval(p3,x), '-.k', 'LineWidth',2);
xlabel('$N$', 'Interpreter','latex');
ylabel('log relative error', 'Interpreter','latex');

%% eigenvalue accuracy plots

load goldstoneN4;

N = 4;
start = 2;
cutoff = 11;

l1p = imag( sqrt( wp(1,start:cutoff) ) );
l1a = imag( sqrt( wa(1,start:cutoff) ) );
err1 = abs( (l1p-l1a)./l1a );

l2p = imag( sqrt( wp(2,start:cutoff) ) );
l2a = imag( sqrt( wa(2,start:cutoff) ) );
err2 = abs( (l2p-l2a)./l2a );

dplot = dvals(start:cutoff);

figure('DefaultAxesFontSize',20);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on;
plot( dplot, log10(err1), '.', 'MarkerSize', 35);
xlabel('$d$', 'Interpreter','latex');
ylabel('$\log_{10}$ relative error', 'Interpreter','latex');

%% eigenvector and kink decays

load intersited05;

figure('DefaultAxesFontSize',20);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on;

v = pi - uk;
vv = v(n/2+2 : end-4 );
nv = 1:length(vv);
p = polyfit( nv, log(vv), 1);
plot( nv, log( vv ), '.b', 'MarkerSize',30, 'LineWidth', 2)

v0 = Vk(:,1);     % Goldstone mode 
w0 = lambdak(1)^2;
vv0 = v0(n/2+2 : end-4 );
nv0 = 1:length(vv0);
p0 = polyfit( nv0, log(vv0), 1);
plot( nv0, log( vv0 ), 'or', 'MarkerSize',10, 'LineWidth', 2)

eindex = 41;      % edge mode
v1 = abs( Vk(:,eindex) );    
w1 = lambdak(eindex)^2;
vv1 = v1(n/2+1 : end-5 );
nv1 = 1:length(vv1);
p1 = polyfit( nv1, log(vv1), 1);
plot( nv1, log( vv1 ), 'xk', 'MarkerSize',10, 'LineWidth', 2)

% regression lines
plot( nv, polyval(p, nv), '-b', 'MarkerSize',30, 'LineWidth', 2)
plot( nv0, polyval(p0, nv0), '--r', 'MarkerSize',10, 'LineWidth', 2)
plot( nv1, polyval(p1, nv1), '-.k', 'MarkerSize',10, 'LineWidth', 2)

xlabel('$N$', 'Interpreter','latex');
ylabel('log of solution', 'Interpreter','latex');
legend({ '$K(n) - \pi$' 'Goldstone mode' 'edge mode' }, 'location', 'southwest', 'Interpreter','latex');

r = (1 + 2*d + sqrt( 1 + 4*d ))/(2*d);
r0 = (1 + w0 + 2*d + sqrt( (1 + w0)*(1 + w0 + 4*d ) ) ) /(2*d);
r1 = (1 + w1 + 2*d + sqrt( (1 + w1)*(1 + w1 + 4*d ) ) ) /(2*d);

err = ( p(1) + log(r) )/ log(r);
err0 = (p0(1) + log(r0))/log(r0);
err1 = (p1(1) + log(r1))/log(r1);

disp([err err0 err1]);
