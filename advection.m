set(0,'defaulttextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex'); 
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultLineLineWidth',3);
set(0,'defaultAxesFontSize',35)

clf
clearvars


%%% Setup grid in x-y
N = 10; % number of grid point along one direction
x = linspace(-5,5,N+1); % type 1 grid
dx = x(2)-x(1);
% convert to type 2 grid
x = x(1:end-1) + dx/2.0;

y = x;
dy = dx;
h = dx;


[X,Y] = meshgrid(x,y); % make 2D grid

% initial condition
F = @(x, y) exp(-(x.^2+y.^2));
u_0 = F(X, Y);

% time stepping parameters
CFL = 0.001;
dt = CFL*(dx);
T_final = 10;
N_steps = T_final/dt;


u = u_0;
for n = 1:N_steps

    u_lap = zeros(size(u));
    % 2D finite laplacian
    for x_i = 1:N
        for y_i = 1:N

            u_center = u(x_i, y_i);
            if x_i == 1 % left edge
                u_west = u_center;
            else
                u_west = u(x_i-1, y_i);
            end
            if x_i == N % right edge
                u_east = u_center;
            else
                u_east=u(x_i+1, y_i);
            end
            if y_i == 1 % bottom edge
                u_south = u_center;
            else
                u_south = u(x_i,y_i-1);
            end
            if y_i == N % top edge
                u_north = u_center;
            else
                u_north=u(x_i, y_i+1);
            end

            % finite difference laplacian w 5 point stencil
            u_lap(x_i, y_i) = (u_west + u_east + u_south + u_north - 4*u_center)*(1/h^2);

        end
    end

    % still needs convection term added below
    u = u + dt*D*u_lap;

end
