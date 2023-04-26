Input = 2D Steady State Conduction

clc; clear all;

Lx = 3; % length of radiator in meters
Ly = 2; % height of radiator in meters
nx = 25; % number of grid points
ny = 25;

dx = Lx/nx; % grid spacing
dy = Ly/ny;

% Define the boundary conditions for the radiator in Celsius
T_left = 80; 
T_right = 90; 
T_top = 60; 
T_bottom = 95; 

% Initialize the temperature matrix
T = zeros(ny,nx);

%% Set the boundary conditions in the temperature matrix
T(:,1) = T_left;
T(:,end) = T_right;
T(1,:) = T_top;
T(end,:) = T_bottom;

% Define the convergence criteria
error = 1;
tol = 1e-6;

% Solve for the temperature distribution using finite difference method
while error > tol
    
    T_old = T;

    for i = 2:nx-1
        for j = 2:ny-1
            T(i,j) = (1/4)*(T_old(i+1,j) + T_old(i-1,j) + T_old(i,j+1) + T_old(i,j-1));
        end
    end

    error = max(max(abs(T - T_old)));
end

%% Generate the x and y coordinates for the contour plot
x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);
[X,Y] = meshgrid(x,y);

% Generate the contour plot of the temperature distribution
contourf(X,Y,T,20,'LineColor','none');
colorbar;
xlabel('Length (m)');
ylabel('Height (m)');
title('Temperature Distribution in Car Radiator');
 

%% To Find the temperature at the specified point using interpolation
x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);
[X,Y] = meshgrid(x,y);
T_interp = interp2(X,Y,T,x_coord,y_coord);

% Display the temperature at the specified point
disp(['Temperature at (', num2str(x_coord), ',', num2str(y_coord), '): ', num2str(T_interp), '°C']);

Output =
Enter the x coordinate: 1.75
Enter the y coordinate: 0.75
Temperature at (1.75,0.75): 78.1433°C
