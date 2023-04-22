set(0,'defaulttextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex'); 
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultLineLineWidth',3);
set(0,'defaultAxesFontSize',35)

clf
clearvars


%%% Setup grid in x-y
N = 100; % number of grid point along one direction
x = linspace(-3,3,N+1); % type 1 grid
dx = x(2)-x(1);
% convert to type 2 grid
x = x(1:end-1) + dx/2.0;

y = x;
dy = dx;
h = dx;


[X,Y] = meshgrid(x,y); % make 2D grid

% initial condition
F = @(x, y) 3*exp(-(x.^2+y.^2));
u_0 = F(X, Y);

A = readmatrix("serial_c_code/out.txt");

norm(A - u_0) % should be small
% for me, only 10e-4. I think that's fine- it's due to the difference in
% precision between matlab and c