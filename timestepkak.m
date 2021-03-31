clear all;

load kink05;
load kak05_8;

n = length(uk);
x = [-n/2:n/2-1];

% difference operators
D1 = fddiffeasy(n, 1, 1, 'none');
D2 = fddiffeasy(n, 2, 1, 'none');
D2(1,2) = 2;
D2(n,n-1) = 2;

u0 = u1;
v0 = 0*uk;

% kak
centers1 = [10 11];
centers2 = [18 19];
u0(10) = u0(10) + 0.1;
u0(11) = u0(11) + 0.1;
u0(18) = u0(18) - 0.1;
u0(19) = u0(19) - 0.1;


ustart = [u0' ; v0'];

f = @(u) par.eps*D2*u + sin(u);

t = linspace(0,150,5000);

F = @(t, u) [ u(n+1:end) ; f( u(1:n) ) ];

u = rk4(F, ustart, t);

uvals = u(1:n,:);

%%

figure('DefaultAxesFontSize',16);
p = waterfall(t,x,uvals);
p.EdgeColor = 'k';
xlabel('t');
ylabel('n');
zlabel('u_n');

%%

udiff = uvals - repmat(uk', 1, length(t));

y1L = u(10,:) - u1(10);
y2L = u(11,:) - u1(11);
y1R = u(18,:) - u1(18);
y2R = u(19,:) - u1(19);

figure('DefaultAxesFontSize',20);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

subplot(2,1,1);
plot( t, y1L, '-b', t, y2L, '--r', 'LineWidth', 2 );
axis([ 0, 50, -0.1, 0.1] );
xlabel('$t$','interpreter','latex');
ylabel('deviation','interpreter','latex');
title('left kink','interpreter','latex');
subplot(2,1,2);
plot( t, y1R, '-b', t, y2R, '--r', 'LineWidth', 2 );
axis([ 0, 50, -0.1, 0.1] );
xlabel('$t$','interpreter','latex');
ylabel('deviation','interpreter','latex');
title('right kink','interpreter','latex');

% plot( t, vecnorm(udiff) );

[pks,locs] = findpeaks( y1R , 'MinPeakHeight', 0.01);
tlocs = t(locs);
td = tlocs(3:10) - tlocs(2:10-1);
mtd = mean(td);




%%

% Runge-Kutta 4 ODE solver
% t is time grid
function u = rk4(f, u0, t)
    u = u0;
    h = t(2) - t(1);
    for index = 1:(length(t) - 1)
       k1 = h*f( t(index), u(:,end) );
       k2 = h*f( t(index)+h/2, u(:,end)+0.5*k1 );
       k3 = h*f( t(index)+h/2, u(:,end)+0.5*k2 ); 
       k4 = h*f( t(index)+h, u(:,end)+k3 );
       u = [ u  u(:,end)+(k1 + 2*k2 + 2*k3 + k4)/6 ];
    end
end