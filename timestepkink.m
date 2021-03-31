clear all;

load kink05;

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
centers = [29 30];
% u0 = uk + 0.1*Vk(:,1)';
% u0 = uk + 0.1*Vk(:,83)';
u0(29) = u0(29) + 0.1;
u0(30) = u0(30) + 0.1;

ustart = [u0' ; v0'];

f = @(u) par.eps*D2*u + sin(u);

t = linspace(0,100,10000);
F = @(t, u) [ u(n+1:end) ; f( u(1:n) ) ];

u = rk4(F, ustart, t);
uvals = u(1:n,:);


%% waterfall plot

figure('DefaultAxesFontSize',20);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
p = waterfall(t,x,uvals);
p.EdgeColor = 'k';
xlabel('$t$','interpreter','latex');
ylabel('$n$','interpreter','latex');
zlabel('$u_n$','interpreter','latex');

%%

udiff = uvals - repmat(uk', 1, length(t));

y1 = u(29,:)  - uk(29);
y2 = u(30,:)  - uk(30);

figure('DefaultAxesFontSize',20);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

plot( t, y1, '-b', t, y2, '--r', 'LineWidth', 2 );
% axis([ 0, 50, -0.1, 0.1] );
xlabel('$t$','interpreter','latex');
ylabel('deviation from primary kink','interpreter','latex');

% plot( t, vecnorm(udiff) );

[pks,locs] = findpeaks( y1 , 'MinPeakHeight', 0.01);
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