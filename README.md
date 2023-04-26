 

COMPUTATIONAL ENGINEERING LAB
CASE STUDY – 2



20BME103	Digvijay Gohil
20BME110	Keval Patel





(B.Tech. Mechanical Engineering)
School of Technology
Pandit Deendayal Energy University
Gandhinagar, Gujarat.
Metal radiator in car engine that undergoes Heat Transfer

The metal radiator is in contact with the coolant flowing through the radiator and that heat is transferred from the metal to the coolant through convection.
Conduction is also a mode of heat transfer that takes place in a radiator. When hot water or steam flows through the radiator, it heats up the metal surface of the radiator through conduction. The heated metal surface then transfers heat to the air surrounding the radiator through convection.

However, the contribution of conduction to the overall heat transfer in a radiator is typically much smaller compared to convection. This is because the thermal conductivity of air is relatively low, and the metal surface of the radiator is usually designed to have a large surface area to enhance convective heat transfer.

The code models the temperature distribution in a metal radiator in a car engine. The metal radiator is assumed to be a rectangular object with a length of Lx in the x-direction and Ly in the y-direction. The code discretizes the radiator into a grid of nx by ny points, with a grid spacing of dx in the x-direction and dy in the y-direction.

The temperature at the boundaries of the radiator are specified as T_left, T_right, T_top, and T_bottom. The code initializes the temperature distribution in the radiator as a matrix of zeros with dimensions ny by nx.

The code then solves for the temperature distribution using an iterative method. The method updates the temperature at each grid point based on the temperatures at neighbouring grid points, until the solution converges to within a tolerance.


Assumptions:

1.	The radiator is a rectangular object with a uniform cross-section
2.	The radiator is made of a homogeneous material
3.	The temperature at the boundaries of the radiator are constant and do not change over time
4.	The temperature distribution in the radiator is in steady-state, i.e. it does not change over time
5.	The radiator is at a constant temperature throughout, with no internal heat generation or sources/sinks


 
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
Parameter	Description	Value
Lx	Length of radiator	3 m
Ly	Height of radiator	2 m
nx	Number of grid points in x-direction	25
ny	Number of grid points in y-direction	25
dx	Grid spacing in x-direction	Lx/nx
dy	Grid spacing in y-direction	Ly/ny
T_left	Temperature at the left boundary	80 °C
T_right	Temperature at the right boundary	90 °C
T_top	Temperature at the top boundary	60 °C
T_bottom	Temperature at the bottom boundary	95 °C
error	Convergence criteria for the finite difference method	1e-6
tol	Tolerance for the convergence of the temperature matrix	1e-6
 

