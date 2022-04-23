%% PURPOSE: Find the geometry area, centroid and inertia moment

%% CLEANING
clear all; close all; clc

%% COORDINATE MATRIX
%   coordinate_matrix = [ nº | X coordinate | Y coordinate ] 

nnos = 31;                              % Number of nodes
nel = nnos - 1;                         % Number of elements

coordinate_matrix = zeros(nnos,3);      % Coordinate matrix pre-allocation
coordinate_matrix(1,1) = 1;             % Node reference
coordinate_matrix(1,2) = 0;             % X coordinate
coordinate_matrix(1,3) = 0;             % Y coordinate

theta = (pi/2)/(nnos-1);

for i = 1:nnos
    coordinate_matrix(i+1,1) = i+1;                                      %   Node reference
    coordinate_matrix(i+1,2) = cos(theta*(i-1))+coordinate_matrix(1,2);  %   X coordinate
    coordinate_matrix(i+1,3) = sin(theta*(i-1))+coordinate_matrix(1,3);  %   Y coordinate
end

figure(1)                                                                %   First window
xlim([0 1])
ylim([0 1])
scatter(coordinate_matrix(:,2),coordinate_matrix(:,3),100,[1 0 0], 'k*') %   Nodes scatter
title('Nodes')
xlabel('X coordinate')
ylabel('Y coordinate')

figure(2)                                                                %  Second window
fill(coordinate_matrix(:,2),coordinate_matrix(:,3),[0 0 0]) 
title('Geometry area')
xlabel('X coordinate')
ylabel('Y coordinate')

%% INCIDENCE MATRIX

incidence_matrix = zeros(nnos-1,3); 

for i = 1:nnos-1 
    incidence_matrix(i,1) = 1;          %	Central node
    incidence_matrix(i,2) = i+1;        %   First node
    incidence_matrix(i,3) = i+2;        %   Second node
end

%% Area

A = 0;                      %   Area somatory

x1 = coordinate_matrix(incidence_matrix(1,1),2);    %   Central node coordinates
y1 = coordinate_matrix(incidence_matrix(1,1),3);

f = zeros(nnos-1,1);

for i = 1:nnos-1
    x2 = coordinate_matrix(incidence_matrix(i,2),2); % X coordinates of the nodes on the edge
    x3 = coordinate_matrix(incidence_matrix(i,3),2);
    
    y2 = coordinate_matrix(incidence_matrix(i,2),3); % Y coordinates of the nodes on the edge
    y3 = coordinate_matrix(incidence_matrix(i,3),3);
    
    a = 0.5*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));        % Element area computation
    
    f(i) = a;
    
    A = A + a;                                        % Area summation of all elements
end

error_Area = ((pi)/4-A)/(pi/4)*100;

%% Centroid 

k = 0;

j = 0;

for i = 1:nnos-1
    x2 = coord(inci(i,2),2); % X coordinates of the nodes on the edge
    x3 = coord(inci(i,3),2);
    
    y2 = coord(inci(i,2),3); % Y coordinates of the nodes on the edge
    y3 = coord(inci(i,3),3);
    
    k = k + ((x1+x2+x3)/3)*f(i);
    
    j = j + ((y1+y2+y3)/3)*f(i); 
end

X = k/A; %centróide X

Y = j/A; %centróide Y

Exact_centroid = (4/(pi*3))+coord(1,2);

error_X = (Exact_centroid-X)/Exact_centroid*100; % Error 
error_Y = (Exact_centroid-Y)/Exact_centroid*100; % Error 

%% Inertia moment

Ix = 0; %   Inertia moment X axis

Iy = 0; %   Inertia moment Y axis

for i = 1:nnos-1
    x2 = coord(inci(i,2),2); % X coordinates of the nodes on the edge
    x3 = coord(inci(i,3),2);
    
    y2 = coord(inci(i,2),3); % Y coordinates of the nodes on the edge
    y3 = coord(inci(i,3),3);
    
    Ix = Ix + 1/12*(y1^2+y2^2+y3^2+y1*y2+y1*y3+y2*y3)*(x2*y3-x3*y2); 
    Iy = Iy + 1/12*(x1^2+x1*x2+x1*x3+x2^2+x2*x3+x3^2)*(x2*y3-x3*y2); 
end

Iexato = pi/16;             % Exact value inertia moment

erro_Ix = (Iexato-Ix)/Iexato*100;
erro_Iy = (Iexato-Iy)/Iexato*100;