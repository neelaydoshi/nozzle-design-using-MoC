%%%%%%%%%%%%%%
% Author : Neelay Doshi
% Project: Using Method of Characteristics for SERN Design
%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
clc
clear 
format short g


%%%%%%%%%%%%%%
% Timer start
tic

%%%%%%%%%%%%%%
% Constants
%%%%%%%%%%%%%%
% only change the variables below 
%%%%%%%%%%%%%%
n   = 25; % number of characteristic lines 
M0  = 1.5; % inlet mach number (>1)
Me  = 3; % exit mach number (desired)
gamma= 1.4; % specific heat ratio


%%%%%%%%%%%%%%
% do NOT alter the code below (unless you know what you are doing)
%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
% turning angle for first ramp
nu_0= PrandtlMeyer(M0, gamma); 
nu_e= PrandtlMeyer(Me, gamma);
theta_turn= (nu_e - nu_0)/2;
theta_n= theta_turn/(n-1);

%%%%%%%%%%%%%%
theta_arr   = 0 : theta_n : theta_turn;
nu_arr      = theta_arr + nu_0;
K_arr       = nu_arr + theta_arr;

%%%%%%%%%%%%%%
K_mat= ones(n);
THETA= zeros(n);
for i= 1:n
    K_mat(i, i:end)= -1;
    THETA(i, i:end)= theta_arr(1:end-(i-1) );
end
THETA= THETA + THETA';


%%%%%%%%%%%%%%
K   = K_arr'.*K_mat; % each row of K_mat multiplied with same element of K_arr
NU  = abs(K-THETA);


%%%%%%%%%%%%%%
MACH= zeros(n);
for i= 1:n
    for j= i:n
        MACH(i, j)=  InversePrandtlMeyer(M0, Me, NU(i, j), gamma);                       
    end
end
MACH= MACH + MACH' - eye(n).*MACH;

MU  = asind(1./MACH);
MU  = K_mat'.*MU + 2*eye(n).*MU;

SLOPE= MU + THETA; %in degrees


%%%%%%%%%%%%%%
% Computing Internal Grid Points 
%%%%%%%%%%%%%%
% 
% Equation of two lines:
% y = m1*x1 + c1
% y = m2*x2 + c2
%
% Computing intersection points using matrix multiplication:
% [-m1, 1] [x] = [c1]
% [-m2, 1] [y] = [c2]
%    [A] x [B] = [C]
%          [B] = inv([A]) x [C]
%
% Extracting "x" and "y" value from matrix "B":
% x = B(1)
% y = B(2)
%

m_ALL= tand(SLOPE);
m_start= -abs( m_ALL(:,1) )' ; 
c_start= ones(n, 1);

m_symmetry= 0;
c_symmetry= 0;

m1= m_start(1);
c1= c_start(1);
m2= m_symmetry;
c2= c_symmetry;

A= [-m1, 1; -m2, 1];
C= [c1; c2];
B= A\C;

X= zeros(n);
X(1, 1)= B(1);
Y= zeros(n);
Y(1, 1)= B(2);

%%%%%%%%%%%%%%
for i= 2:n
    m1= m_start(i);
    c1= c_start(i);
    m2= m_ALL(1, i-1) ;
    c2= Y(1, i-1) - m_ALL(1, i-1)*X(1, i-1);
    
    A= [-m1, 1; -m2, 1];
    C= [c1; c2];

    B= A\C;
    X(1, i)= B(1);
    Y(1, i)= B(2);
end

X= X + X' - eye(n).*X;
Y= Y + Y' - eye(n).*Y;

%%%%%%%%%%%%%%
m_symmetry= 0;
c_symmetry= 0;
for i= 2:n
    for j= i:n
        if i == j % diagonal 
            m1= m_ALL(i, j-1);
            c1= Y(i, j-1) - m1*X(i, j-1);
            m2= m_symmetry;
            c2= c_symmetry;

            A= [-m1, 1; -m2, 1];
            C= [c1; c2];

            B= A\C;
            X(i, j)= B(1);
            Y(i, j)= B(2);
            
        else
            m1= m_ALL(i, j-1);
            c1= Y(i, j-1) - m1*X(i, j-1);
            m2= m_ALL(j, i-1);
            c2= Y(j, i-1) - m2*X(j, i-1);
            
            A= [-m1, 1; -m2, 1];
            C= [c1; c2];

            B= A\C;
            X(i, j)= B(1);
            Y(i, j)= B(2);
            
            X(j, i)= X(i, j);
            Y(j, i)= Y(i, j);
        end

    end
end


%%%%%%%%%%%
% ROOF
m_ROOF= zeros(1, n);
c_ROOF= zeros(1, n);

X_ROOF= zeros(1, n);
Y_ROOF= zeros(1, n);

%%%%%%%%%%%
% First Intersection Point
m_ROOF(1)= tand(theta_turn);
c_ROOF(1)= 1;

m1= m_ROOF(1);
c1= c_ROOF(1);
m2= m_ALL(1, end) ;
c2= Y(1, end) - m2*X(1, end);

A= [-m1, 1; -m2, 1];
C= [c1; c2];

B= A\C;
X_ROOF(1)= B(1);
Y_ROOF(1)= B(2);


%%%%%%%%%%%%%%
for i= 2:n
    m_ROOF(i)= tand( THETA(i-1, end) );
    c_ROOF(i)= Y_ROOF(i-1) - m_ROOF(i)*X_ROOF(i-1);
    
    m1= m_ROOF(i);
    c1= c_ROOF(i);
    m2= m_ALL(i, end) ;
    c2= Y(i, end) - m2*X(i, end);
    
    A= [-m1, 1; -m2, 1];
    C= [c1; c2];

    B= A\C;
    X_ROOF(i)= B(1);
    Y_ROOF(i)= B(2);
end


%%%%%%%%%%%
X_ALL= zeros(n, n+2);
Y_ALL= zeros(n, n+2);

X_ALL(:, 1)= 0;
Y_ALL(:, 1)= 1;

X_ALL(:, 2:n+1)= X;
Y_ALL(:, 2:n+1)= Y;

X_ALL(:, end)= X_ROOF';
Y_ALL(:, end)= Y_ROOF';


%%%%%%%%%%%%%%
% Printing 
fprintf('############## \n');
toc % Timer off
fprintf('Number of characteristic lines: %d \n', n);
fprintf('Nozzle exit coordinates (x, y): (%.4f, %.4f)\n', X_ROOF(end), Y_ROOF(end));
fprintf('############## \n');


%%%%%%%%%%%%%
% Pressure Profile and Force Calculation

Ps_inlet    = 1.25e5; % static pressure at inlet
Pt_inlet    = Ps_inlet*(1 + (gamma-1)/2*M0^2)^(gamma/(gamma-1)); % total pressure at inlet

theta_roof  = atand(m_ROOF);
mach_roof   = [MACH(1, end), MACH(1:end-1, end)'];

term_1      = 1 + (gamma-1)/2.*mach_roof.^2;
term_2      = gamma/(gamma-1);
Ps_roof     = Pt_inlet ./ ( term_1.^term_2 );

x_roof  = [0, X_ROOF];
y_roof  = [1, Y_ROOF];
dx_roof = ( x_roof(2:end) - x_roof(1:end-1) ).^2;
dy_roof = ( y_roof(2:end) - y_roof(1:end-1) ).^2;
len_roof= sqrt( dx_roof + dy_roof );

Fx_roof = -Ps_roof.*(len_roof.*1).*sind(theta_roof);
Fy_roof = Ps_roof.*(len_roof.*1).*cosd(theta_roof);

scaling = 10;
Fx = sum(Fx_roof)/scaling;
Fy = sum(Fy_roof)/scaling;
fprintf('############## \n');
fprintf('Force in x-direction on nozzle roof = %.3f kN \n', Fx/1e3);
fprintf('Force in y-direction on nozzle roof = %.3f kN \n', Fy/1e3);
fprintf('############## \n');



%%%%%%%%%%%%%
% Plotting Characteristic Lines
figure()
for i= 1:n
    plot(X_ALL(i, :), Y_ALL(i, :), 'b');
    hold on
end
axis equal

%%%%%%%%%%%%%%
% Plotting ROOF and BASE
%
% Base
plot([-1, X_ROOF(end)], [0, 0], '-k', LineWidth= 2)
hold on
%
% Inlet Roof
plot([-1, 0], [1, 1], '-k', LineWidth= 2); 
%
% Nozzle Contour Roof 
plot([0, X_ROOF(1)], [1, Y_ROOF(1)], '-k', LineWidth= 2);
plot(X_ROOF, Y_ROOF, 'k', LineWidth= 2);


%% %%%%%%%%%%%%%%
% Generating (x, y, z) coordinates for plotting in ICEM-CFD
% manually add the first coordinate (0, 1, 0)
Z= zeros(n, 1);
nozzle_coord= [X_ALL(:, end), Y_ALL(:, end), Z];



