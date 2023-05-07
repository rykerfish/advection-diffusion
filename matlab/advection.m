set(0,'defaulttextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex'); 
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultLineLineWidth',3);
set(0,'defaultAxesFontSize',35)

clf
clearvars


%%% Setup grid in x-y
N = 200; % number of grid point along one direction
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
% u_0 = ic(X, Y);



% time stepping parameters
% CFL = 0.001;
CFL = 1/(10*max(max(u_0)));
dt = CFL*(dx);
T_final = 50;
N_steps = T_final/dt;

% physical parameters
D = 0.01; % diffusion
vel = @(x) 0.1;



t_plot = 1000;

u = u_0;
for n = 1:N_steps

    u_lap = zeros(size(u));
    u_adv = zeros(size(u));
    for x_i = 1:N
        for y_i = 1:N

            % 2D finite laplacian with periodic BC
            u_center = u(x_i, y_i);
            if x_i == 1 % left edge
                u_west = u(N, y_i);
            else
                u_west = u(x_i-1, y_i);
            end
            if x_i == N % right edge
                u_east = u(1, y_i);
            else
                u_east=u(x_i+1, y_i);
            end
            if y_i == 1 % bottom edge
                u_south = u(x_i, N);
            else
                u_south = u(x_i,y_i-1);
            end
            if y_i == N % top edge
                u_north = u(x_i, 1);
            else
                u_north=u(x_i, y_i+1);
            end

            % finite difference laplacian w 5 point stencil
            u_lap(x_i, y_i) = (u_west + u_east + u_south + u_north - 4*u_center)*(1/h^2);


            % periodic BC for advection
            if x_i == 1
                u_left = u(N, y_i);
            else
                u_left = u(x_i-1, y_i);
            end
            if x_i == N
                u_right = u(1, y_i);
            else
                u_right = u(x_i+1, y_i);
            end

             % upwinding
            if(vel(x_i) > 0)
                u_adv(x_i, y_i) = vel(x_i)*(u(x_i, y_i) - u_left)/dx;
            else
                u_adv(x_i, y_i) = vel(x_i)*(u_right - u(x_i, y_i))/dx;
            end
                % centered difference
%                 u_adv(x_i, y_i) = vel(x_i)*(u(x_i+1, y_i) - u(x_i-1,y_i))/(2*dx); % centered difference
        end 
    end

    % forward euler
    u = u + dt*D*u_lap + u_adv*dt;

    if mod(n,t_plot) == 0
       surface = surf(X,Y,reshape(u,N,N));
       set(surface,'edgecolor','none')
       zlim([0 10])
       title(['t = ' num2str(n*dt)])
       colorbar
       camlight
       drawnow
    end

end

% makes a cylinder as the initial condition
function u_0 = ic(X, Y)

    r = 1;
    u_0 = zeros(size(X));
    
    N = length(X);
    for i = 1:N
        for j = 1:N
            if X(i, j)^2 + Y(i,j)^2 <= r
                u_0(i, j) = 1;
            else
                u_0(i, j) = 0;
            end
        end
    end


end