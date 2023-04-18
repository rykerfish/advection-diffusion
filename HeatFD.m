set(0,'defaulttextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex'); 
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultLineLineWidth',3);
set(0,'defaultAxesFontSize',35)


clf
clear all
%%% Setup grid in x-y
N = 200; % number of grid point along one direction
x = linspace(0,1,N+1); % type 1 grid
dx = x(2)-x(1);
% convert to type 2 grid
x = x(1:end-1) + dx/2.0;

y = x;
dy = dx;


[X,Y] = meshgrid(x,y); % make 2D grid

%%% make finite difference laplacian

e = (1/dx^2)*ones(N,1);
D2 = spdiags([e -2*e e], -1:1, N, N);
% neumann boundary conditions
D2(1,1) = -(1/dx^2);
D2(end,end) = -(1/dx^2);
% periodic boundary conditions
% D2(1,end) = (1/dx^2)*1;
% D2(end,1) = (1/dx^2)*1;


I_n = speye(N);

D_full = kron(I_n, D2) + kron(D2, I_n);

%%% 
epsilon = 0.1; % interface width
CFL = 0.001;
dt = CFL*(dx);
tFinal = 50.0;
Nsteps = round(tFinal/dt);
t_plot = 100;

% Lp = speye(N*N) + (dt/2.0)*epsilon*D_full;
% Lm = speye(N*N) - (dt/2.0)*epsilon*D_full;

Lm = speye(N*N) - (dt)*epsilon*D_full;

dLm = decomposition(Lm);


U = 8*exp(-200*((X-0.5).^2+(Y-0.5).^2));
U = U(:);
Uprev = U;

F = @(u) 0;
img = 0;
for n = 1:Nsteps
    RHS = U;
    U = dLm\RHS;
    if mod(n,t_plot) == 0
       h = surf(X,Y,reshape(U,N,N));
       set(h,'edgecolor','none')
       zlim([0 10])
       title(['t = ' num2str(n*dt)])
       colorbar
       camlight
       drawnow
%        img = img + 1;
%        print('-dpng',['./pngs/heat_frame_' num2str(img) '.png'],'-r100')
    end
end












