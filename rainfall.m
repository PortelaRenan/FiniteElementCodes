%% Program: RAINFALL 
%
%   @DESCRIPTION: Define the amount of rainfall 
%   @AUTHOR: Renan Miranda Portela
%
%% Clean memory and close windows
clear all
close all
clc

%% COORDINATE MATRIX
%   coord = [number | X-position | Y-position] 

coord = [1 0 33.3;
         2 13.2 62.3;
         3 39.3 84.5;
         4 22.2 30.1;
         5 49.9 57.6;
         6 78.8 78.2;
         7 39.3 10.0;
         8 59.7 34.3;
         9 73.9 36.2;
         10 69.8 5.1;
         11 28.0 50.0;
         12 33.3 55.0;
         13 45.0 49.0;
         14 35.0 30.0];
     
scatter(coord(:,2),coord(:,3),10,[1 0 0]) % presenting the nodes
     
%% INCIDENCE MATRIX
% inci = [number of the element | node-1 | node-2 | node-3 | node-4 ]

inci = [1 1 4 5 2;
        2 2 5 6 3;
        3 4 7 8 5;
        4 5 8 9 6;
        5 7 10 9 8;
        6 11 14 13 12];

%% RAINFALL MATRIX
u = [4.62 3.81 4.76 5.45 4.90 10.35 4.96 4.26 18.36 15.69];
    
%% COLOR MATRIX
cores = [1 0 0;
         0 1 0;
         0 0 1;
         1 1 0;
         0 0 0;
         0 1 1];

%% AREA PLOTS     
 for i = 1:size(inci,1)
     x = [coord(inci(i,2),2),coord(inci(i,3),2),coord(inci(i,4),2),coord(inci(i,5),2)];
     y = [coord(inci(i,2),3),coord(inci(i,3),3),coord(inci(i,4),3),coord(inci(i,5),3)];
     patch(x,y,cores(i,:))
 end
 
title('Elements areas')
xlabel('distance in [km]')
ylabel('distance in [km]')
hold on 
scatter(70,40,30,'*')

%% TOTAL AREA
A = 0;            % STARTING SUM
a_i = zeros(5,1); % ELEMENT AREA VALUE 

for i = 1:5                  % COMPUTING THE AREA BY THE JACOBIAN 
    x1 = coord(inci(i,2),2); % X-COORDINATES 
    x2 = coord(inci(i,3),2);
    x3 = coord(inci(i,4),2);
    x4 = coord(inci(i,5),2);
    
    y1 = coord(inci(i,2),3); % Y-COORDINATES
    y2 = coord(inci(i,3),3);
    y3 = coord(inci(i,4),3);
    y4 = coord(inci(i,5),3);
    
    a = 0.5*(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3);
    a_i(i) = a;
    A = A + a;
end

A2 = 0;

for i = 1:5                  % COMPUTING AREA BY GAUSSIAN QUADRATURE 
    x1 = coord(inci(i,2),2); % X-COORDINATES
    x2 = coord(inci(i,3),2);
    x3 = coord(inci(i,4),2);
    x4 = coord(inci(i,5),2);
    x = [x1,x2,x3,x4];
    
    y1 = coord(inci(i,2),3); % Y-COORDINATES
    y2 = coord(inci(i,3),3);
    y3 = coord(inci(i,4),3);
    y4 = coord(inci(i,5),3);
    y = [y1,y2,y3,y4];
    
    for j = 1:2
        if j == 1
            e = -0.57735;
        else 
            e = 0.57735;
        end
        for m = 1:2
            if m == 1
                n = -0.57735;
            else 
                n = 0.57735;
            end
            
            a = 0.25*[-(1-n),(1-n),(1+n),-(1+n);
                      -(1-e),-(1+e),(1+e),(1-e)];
            b(:,1) = x';
            b(:,2) = y';
            
            J = a*b;
%             [ J ] = jacobiano_4P( e,n,x,y );
            
            A2 = A2 + det(J);
        end
    end
end

%% TOTAL RAINFALL
 
Q = 0;           % TOTAL RAINFALL VARIABLE
g = zeros(5,1);

for i = 1:5                  % BY THE JACOBIAN
    x1 = coord(inci(i,2),2); % X-COORDINATES
    x2 = coord(inci(i,3),2);
    x3 = coord(inci(i,4),2);
    x4 = coord(inci(i,5),2);
        
    y1 = coord(inci(i,2),3); % Y-COORDINATES
    y2 = coord(inci(i,3),3);
    y3 = coord(inci(i,4),3);
    y4 = coord(inci(i,5),3);
    
    u1 = u(inci(i,2));
    u2 = u(inci(i,3));
    u3 = u(inci(i,4));
    u4 = u(inci(i,5));
    
    q = (4*(u1/4 - u2/4 - u3/4 + u4/4)*((x1*y3)/8 - (x3*y1)/8 - (x1*y4)/8 - (x2*y3)/8 + (x3*y2)/8 + (x4*y1)/8 + (x2*y4)/8 - (x4*y2)/8))/3 + (4*(u1/4 + u2/4 - u3/4 - u4/4)*((x1*y2)/8 - (x2*y1)/8 - (x1*y3)/8 + (x3*y1)/8 + (x2*y4)/8 - (x4*y2)/8 - (x3*y4)/8 + (x4*y3)/8))/3 + 4*(u1/4 + u2/4 + u3/4 + u4/4)*((x1*y2)/8 - (x2*y1)/8 - (x1*y4)/8 + (x2*y3)/8 - (x3*y2)/8 + (x4*y1)/8 + (x3*y4)/8 - (x4*y3)/8);
    g(i) = q;
    Q = Q + q;
end

fprintf('\n\n******* TOTAL RAINFALL *******\n')
fprintf('        %f\n',Q)

%% AVERAGE RAINFALL

MQ = Q/A;

fprintf('\n\n******* AVERAGE RAINFALL *******\n')
fprintf('       %f\n',MQ)

%% INTERPOLATION

x=[49.9,59.7,73.9,78.8]; % POINTS TO INTERPOLATE AROUND X
y=[57.6,34.3,36.2,78.2]; % POINTS TO INTERPOLATE AROUND Y

X=70; % REFERENCE POINT
Y=40;

U_Ponto_Interpolado = u(5)*n(1)+u(8)*n(2)+u(9)*n(3)+u(6)*n(4);

fprintf('\n\n******* AMOUNT OF RAINFALL AT (70,40) *******\n')
fprintf('       %f\n',U_Ponto_Interpolado)

%% RAINFALL IN REGION A

A = zeros(2,4);

for i = 1:4
    A(1,i) = coord(inci(6,i+1),2); % X-COORDINATES
end

for j = 1:4
    A(2,j) = coord(inci(6,j+1),3); % Y-COORDINATES
end
 
B = zeros(2,4);

for i = 1:4
    B(1,i) = coord(inci(1,i+1),2); 
end

for j = 1:4
    B(2,j) = coord(inci(1,j+1),3); 
end

C = zeros(2,4);

for i = 1:4
    C(1,i) = coord(inci(3,i+1),2); 
end

for j = 1:4
    C(2,j) = coord(inci(3,j+1),3); 
end

x = B(1,:); % ELEMENT 1 REGION
y = B(2,:);

X = A(1,1);
Y = A(2,1);

[ n ] = NR_search( x,y,X,Y ); % A1 POINT

U_A1 = u(1)*n(1)+u(4)*n(2)+u(5)*n(3)+u(2)*n(4);

X = A(1,4);
Y = A(2,4);

[ n ] = NR_search( x,y,X,Y ); % A4 POINT

U_A4 = u(1)*n(1)+u(4)*n(2)+u(5)*n(3)+u(2)*n(4);

x = C(1,:); %região do elemento 3
y = C(2,:);

X = A(1,2);
Y = A(2,2);

[ n ] = NR_search( x,y,X,Y ); % A4 POINT

U_A2 = u(4)*n(1)+u(7)*n(2)+u(8)*n(3)+u(5)*n(4); % A2 POINT

X = A(1,3);
Y = A(2,3);

[ n ] = NR_search( x,y,X,Y ); % A4 POINT

U_A3 = u(4)*n(1)+u(7)*n(2)+u(8)*n(3)+u(5)*n(4); % A3 POINT

x1 = A(1,1); x2 = A(1,2); x3 = A(1,3); x4 = A(1,4);
y1 = A(2,1); y2 = A(2,2); y3 = A(2,3); y4 = A(2,4);

Q_A = (4*(U_A1/4 - U_A2/4 - U_A3/4 + U_A4/4)*((x1*y3)/8 - (x3*y1)/8 - (x1*y4)/8 - (x2*y3)/8 + (x3*y2)/8 + (x4*y1)/8 + (x2*y4)/8 - (x4*y2)/8))/3 + (4*(u1/4 + u2/4 - u3/4 - u4/4)*((x1*y2)/8 - (x2*y1)/8 - (x1*y3)/8 + (x3*y1)/8 + (x2*y4)/8 - (x4*y2)/8 - (x3*y4)/8 + (x4*y3)/8))/3 + 4*(u1/4 + u2/4 + u3/4 + u4/4)*((x1*y2)/8 - (x2*y1)/8 - (x1*y4)/8 + (x2*y3)/8 - (x3*y2)/8 + (x4*y1)/8 + (x3*y4)/8 - (x4*y3)/8);

Area_A = 0.5*(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3);

MQ_A = Q_A/Area_A;

fprintf('\n\n******* TOTAL AMOUNT OF RAINFALL IN ELEMENT A *******\n')
fprintf('       %f\n',Q_A)

fprintf('\n\n******* AVERAGE AMOUNT OF RAINFALL IN ELEMENT A *******\n')
fprintf('       %f\n',MQ_A)
