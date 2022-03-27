%% Program: Torsion
%
%   @DESCRIPTION: Define the different vibration modes of a free-free or free-clamped bar
%   @AUTHOR: Renan Miranda Portela
%   
%% Clen memory and close all windows open
clear all
close all
clc

%% COORDINATE MATRIX
%   coord = [number | X-position |  Y-position] 

node = 49;             % number of nodes
L = 1;                 % length

coord = zeros(node-1,3);

for i = 1:node
  coord(i,1) = i;                % node number
  coord(i,2) = L/(node-1)*(i-1); % X-coordinate
  coord(i,3) = 0; % Y-coordinate
end


%% INCIDENCE MATRIX
% inci = [number | material | geometry | node-1 | node-2 ]

inci = zeros(node-1,5); 
for i = 1:node-1
  inci(i,1) = i;           % node number
  inci(i,2) = 1;           % material
  inci(i,3) = 1;           % geometry
  inci(i,4) = i;           % node-1
  inci(i,5) = i +1;        % node-2
end

%% MATERIALS TABLE
%  Tmat = [ Material 1 | Material 2 | Material 3 ...]
%
%  column 1: Young's modulus
%  column 2: Poisson ratio
  
tabmat = [210E9 0.33];

%% GEOMETRY TABLE
%  Tgeo = [ Geometria 1 | Geometria 2 | Geometria 3 ...]
%
%  column 1: area - square meter
%  column 2: density - kg/m³

tabgeo = [0.1 7860];

%% BOUNDARY CONDITIONS
%   bc = [node | degree of freedom (DF) | value]
%
%   DF 1 --> x
%   DF 2 --> y
%   DF 3 --> z
%   DF 4 --> ox
%   DF 5 --> oy
%   DF --> oz

% bc=[];          % Free-Free bar
bc=[1 1 0];       % Clamped-Free bar

%% EXTERNAL LOAD
%   F = [node | DF | value]
%
%   DF 1 --> Fx
%   DF 2 --> Fy
%   DF 3 --> Fz
%   DF 4 --> Mx
%   DF 5 --> My
%   DF 6 --> Mz

Load=[];

%   VECTOR SIZES
nbc = size(bc,1);  % number of boundary conditions
nF = size(Load,1); % number of external loads
neq = 0;           % number of equations
ngdl = 1;          % number of degrees of freedom


%%  Assembly of global mass and stiffness matrices
id=ones(1,node);

for i=1:nbc
    id(bc(i,2),bc(i,2))=bc(i,3);
end

for i= 1:node
    for j = 1:ngdl
        if id(j,i)== 1
            neq = neq +1;
            id(j,i) = neq;
        end
    end
end

kg = zeros(neq,neq); % pre-allocation of global stiffness matrix
mg = zeros(neq,neq); % pre-allocation of global mass matrix

for i = 1:node-1     % matrices assembly
    no1 = inci(i,4);
    no2 = inci(i,5);
    x1 = coord(no1,2);
    x2 = coord(no2,2);
    l = x2 - x1;
    mat = inci(i,2);
    geo = inci(i,3);
    E = tabmat(1,mat);
    a = tabgeo(1);
    rho = tabgeo(2);
    ke = 1/l*[1 -1; -1 1];       % local stiffness matrix
    me = 1*l/6*[2 1;1 2];        % local mass matrix
    loc = [id(1,no1),id(1,no2)]; 
    for j = 1:2
        if loc(j) ~= 0;
            for k = 1:2
                if loc(k) ~=0
                    kg(loc(j),loc(k))=kg(loc(j),loc(k))+ke(j,k);
                    mg(loc(j),loc(k))=mg(loc(j),loc(k))+me(j,k);
                end
            end
        end
    end
end

[theta,D]=eig(kg,mg); % theta = vibration mode

for i = 1:size(theta,1)
    theta(:,i) = theta(:,i)/max(abs(theta(:,i)));    
end

lambda = diag(D); % eigenvalues [K] - lambda*[M] = 0 

poisson = tabmat(2);

G = E/(2*(1+poisson));
 
b = 1;%sqrt(G/rho);

omega = sqrt(lambda)*b; % vibration frequency

coord_2 = coord; %coord_2 = node coordinate after vibration

y_1 = zeros(node,1);

if isempty(bc)
    
    for i = 1 : node-1 % first mode of vibration
        coord_2(i+1,3) = coord_2(i+1,3) + theta(i,2);
        y_1(i+1) = sin((pi*coord_2(i+1,2))/2);
    end
    
    figure(1)
    plot(coord_2(:,2),coord_2(:,3),'*-')
    
    title('Modos de vibração')
    xlabel('Comprimento da barra sob torção (m)') % x-axis label
    ylabel('Deslocamento vertical do nó')         % x-axis label

    hold on
    
    coord_2 = coord; 

    if node >= 3

        y_2 = zeros(node,1);

        for i = 1 : node-1 % second mode of vibration
            coord_2(i+1,3) = coord_2(i+1,3) - theta(i,3);
            y_2(i+1) = sin((3*pi*coord_2(i+1,2))/2);
        end

        erro_2 = norm(abs(coord_2(:,3))-abs(y_2));

        plot(coord_2(:,2),coord_2(:,3),'k')

        coord_2 = coord; 
    end

    if node >= 4

        y_3 = zeros(node,1);

        for i = 1 : node-1 % third mode of vibration
            coord_2(i+1,3) = coord_2(i+1,3) + theta(i,4);
            y_3(i+1) = sin((5*pi*coord_2(i+1,2))/2);
        end

        erro_3 = norm(abs(coord_2(:,3))-abs(y_3));

        plot(coord_2(:,2),coord_2(:,3),'r--')

        legend('Primeiro','Segundo','Terceiro')

        coord_2 = coord; 

    end
    
else

    for i = 1 : node-1 % first mode of vibration
        coord_2(i+1,3) = coord_2(i+1,3) - theta(i,1);
        y_1(i+1) = sin((pi*coord_2(i+1,2))/2);
    end

    erro = norm(abs(coord_2(:,3))-abs(y_1));

    figure(1)
    plot(coord_2(:,2),coord_2(:,3),'*-')

    title('VIbration modes')
    xlabel('Bar length under torsion (m)') % x-axis label
    ylabel('Vertical node displacement (m)') % y-axis label

    hold on

    coord_2 = coord; 

    if node >= 3

        y_2 = zeros(node,1);

        for i = 1 : node-1 % second vibration mode
            coord_2(i+1,3) = coord_2(i+1,3) + theta(i,2);
            y_2(i+1) = sin((3*pi*coord_2(i+1,2))/2);
        end

        erro_2 = norm(abs(coord_2(:,3))-abs(y_2));

        plot(coord_2(:,2),coord_2(:,3),'k')

        coord_2 = coord; 
    end

    if node >= 4

        y_3 = zeros(node,1);

        for i = 1 : node-1 % third vibration mode
            coord_2(i+1,3) = coord_2(i+1,3) - theta(i,3);
            y_3(i+1) = sin((5*pi*coord_2(i+1,2))/2);
        end

        erro_3 = norm(abs(coord_2(:,3))-abs(y_3));

        plot(coord_2(:,2),coord_2(:,3),'r--')

        legend('Primeiro','Segundo','Terceiro')

        coord_2 = coord;

    end
end

fprintf('\n\n******* Eigenvalues *******\n')
fprintf('       %f\n',lambda)
fprintf('\n\n******* Vibration frequencies *******\n')
fprintf('       %f\n',omega)
fprintf('\n\n******* First vibration mode *******\n')
fprintf('       %f\n',theta(:,1))
