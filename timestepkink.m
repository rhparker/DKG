clear all;

load kink05_200;
% load data/kink05;
% load kink070_200;

n = length(uk);
x = [-n/2:n/2-1];

% difference operators
D1 = fddiffeasy(n, 1, 1, 'none');
D2 = fddiffeasy(n, 2, 1, 'none');
D2(1,2) = 2;
D2(n,n-1) = 2;

u0 = uk;
v0 = 0*uk;

% kink
centers = [n/2-1 n/2];
% centers = [14 16];
u0 = uk + 0.1*Vk(:,1)';
% u0 = uk + 0.1*Vk(:,83)';
% u0(centers(1)) = u0(centers(1)) + 0.1;
% u0(centers(2)) = u0(centers(2)) + 0.1;

ustart = [u0' ; v0'];

f = @(u) par.eps*D2*u + sin(u);

t = linspace(0,2000,1000);
h = t(2) - t(1);

F = @(t, u) [ u(n+1:end) ; f( u(1:n) ) ];

% u = rk4(F, ustart, t);

usol = ode45( F, [0 2000], ustart);

y1 = usol.y(centers(1),:)  - uk(centers(1));
y2 = usol.y(centers(2),:)  - uk(centers(2));
tplot = usol.x;

figure('DefaultAxesFontSize',20);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

plot( tplot, y1, '-b', tplot, y2, '--r', 'LineWidth', 2 );
% axis([ 0, 50, -0.1, 0.1] );
xlabel('$t$','interpreter','latex');
ylabel('deviation from primary kink','interpreter','latex');

% uvals = u(1:n,:);


%% waterfall plot

% figure('DefaultAxesFontSize',20);
% set(gca,'fontname','times');
% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% p = waterfall(t,x,uvals);
% p.EdgeColor = 'k';
% xlabel('$t$','interpreter','latex');
% ylabel('$n$','interpreter','latex');
% zlabel('$u_n$','interpreter','latex');

%%

udiff = uvals - repmat(uk', 1, length(t));
udiffnorms = vecnorm(udiff);

y1 = u(centers(1),:)  - uk(centers(1));
y2 = u(centers(2),:)  - uk(centers(2));

figure('DefaultAxesFontSize',20);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

plot( t, y1, '-b', t, y2, '--r', 'LineWidth', 2 );
% axis([ 0, 50, -0.1, 0.1] );
xlabel('$t$','interpreter','latex');
ylabel('deviation from primary kink','interpreter','latex');

% norm plot
figure('DefaultAxesFontSize',20);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

plot( t, udiffnorms, '-b', 'LineWidth', 2 );
xlabel('$t$','interpreter','latex');
ylabel('$\| u(n, t) - u_1(n) \|$','interpreter','latex');

%%

[pks,locs] = findpeaks( y1 , 'MinPeakHeight', 0.01);
tlocs = t(locs);
td = tlocs(2:end) - tlocs(1:end-1);
mtd = mean(td);

start = 1;
tstart = tlocs(start);

x = tlocs(start:end);
y = pks(start:end);

lx = log(x);
ly = log(y);

% p = polyfit(lx,ly,1);

% figure('DefaultAxesFontSize',20);
% set(gca,'fontname','times');
% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% hold on;

% tindex = 20;
% t0 = tlocs(tindex);
% 
% lx2 = log( x/t0 - 1 );
% lx2 = lx2(tindex:end);
% ly = ly(tindex:end);

% plot( lx2, ly, '.', 'MarkerSize',30) 
% % plot( lx, polyval(p,lx), 'LineWidth',2 );
% xlabel('$\log t$','interpreter','latex')
% ylabel('$\log$ amplitude','interpreter','latex')
% 
% p2 = polyfit(lx2,ly,1);



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

