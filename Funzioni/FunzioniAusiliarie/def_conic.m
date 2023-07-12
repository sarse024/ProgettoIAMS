function [f_conic, e, a, b, pos_c, th_rotation] = def_conic(f1,f2,a)
% def_conic.m - Function to calculate the parametre of a conic from 2 foci and semimajor
% axis (a)
%
% PROTOTYPE:
% [f_conic, e, a, b, pos_c, th_rotation] = def_conic(f1,f2,a)
%
% DESCRIPTION:
% The function generate the conic equation and calculate the most importatn parametre. 
% Pos_c and th_rotation are useful to rotate and translate the conic in the correct position. 
% Angles in degree
%
% INPUT:
% f1    [2x1]   First focus in (x,y)                                            [-]
% f2    [2x1]   Second focus in (x,y)                                           [-]                       
% a     [1x1]   Semimajor axis                                                  [-]
%
% OUTPUT:
% f_conic       [-]     Function of the conic section in quadratic formula      [-]
% e             [1x1]   Eccentricity                                            [-]
% a             [1x1]   Semi-major axis                                         [-]
% b             [1x1]   Semi-minor axis                                         [deg]
% pos_c         [1x2]   Position of the center of conic                         [-]
% th_rotation   [1x1]   Rotation of the conic in the plane (clockwise angle)    [deg]
    
    % formula to find the coefficients of the conic genearl equation with
    % given foci and a
    A = 16*a^2-4*(f1(1) - f2(1))^2;
    B = -8*(f1(1) - f2(1))*(f1(2) - f2(2));
    C = 16*a^2-4*(f1(2) - f2(2))^2;
    D = 4*(f1(1) - f2(1))*(f1(1)^2 - f2(1)^2 + f1(2)^2 - f2(2)^2) - 16*(a^2)*(f1(1) + f2(1));
    E = 4*(f1(2) - f2(2))*(f1(1)^2 - f2(1)^2 + f1(2)^2 - f2(2)^2) - 16*(a^2)*(f1(2) + f2(2));
    F = 4*(f1(1)^2 + f1(2)^2)*(f2(1)^2 + f2(2)^2) - (f1(1)^2 + f2(1)^2 + f1(2)^2 + f2(2)^2 - 4*a^2)^2;
    
    % define the conic equation
    f_conic = @(x,y) A*x.^2 + B*x.*y + C*y.^2 + D*x + E*y + F;
    
    matrix = [A B/2 D/2; B/2 C E/2; D/2 E/2 F];
    Q = [A B D; B C E; D E F];

    if(det(Q) ~= 0) %non denegere conic
        if(det(matrix)>0)
            eta = -1;
        else
            eta = 1;
        end
        e = sqrt(2*sqrt((A-C)^2 + B^2)/(eta*(A+C) + sqrt((A-C)^2+B^2)));
    else
        e = -1; % test to indicate a degenere conic
    end

    if(B^2-4*A*C>0) % hyperbola case
        a = -a;
    end
    
    % calculate the position of the center
    K = det([A,B/2;B/2,C]);
    xc = -1/K * det([D/2,B/2; E/2,C]);
    yc = -1/K * det([A,D/2; B/2,E/2]);
    pos_c = [xc,yc];
    
    % return b with geometrical form
    b = abs(a)*sqrt(e.^2 -1);

    % Avoid problem of rotation in vertical hyperbole
    correction = 0; % angle in deg
    vect_pos = f1 - f2;
    test = rad2deg(acos(vect_pos(1)/norm(vect_pos)));
    if(abs(test) > 45 && abs(test) < 135)
        correction = 90; % angle in deg
    end

    th_rotation = 1/2*rad2deg(acot((A-C)/B)) + correction; % clockwise angle

end