function [kepT, dv, th] = findOrbit(r1, r2, kep1, kep2)
% findOrbit.m - Function to calculate the most efficienty elliptical orbit
% that connect two point in two define orbite
%
% PROTOTYPE:
% [kepT, dv, th] = findOrbit(r1, r2, kep1, kep2)
%
% DESCRIPTION:
% The function calculates different elliptical orbits that connect two point on two different orbit and calculates the dvs in the two maneuvering points using the vector difference. 
% The function returns the orbit with the minimum total dv (dv1 + dv2) of
% the two pulses. If the algorithm don't find any ellipse with the correct precision without crashing the Earth, return all zero. Angles in degree
%
% INPUT:
% r1    [3x1]   First position vector                               [km]
% r2    [3x1]   Second position vector                              [km]
% kep1  [1x6]   Orbital Parametre of fist orbit                     [a,e,i,OM,om,th]
% kep2  [1x6]   Orbital Parametre of second orbit                   [a,e,i,OM,om,th]
%
% OUTPUT:
% kepT  [1x6]   Orbital Parametre of transfert orbit                [a,e,i,OM,om,th]
% dv    [1x1]   Minimum delta v                                     [km/s]
% th    [1x4]   Contains: [th1, new_th1, new_th2, th2]              [deg]
% - th1     is the true anomlay of r1 in kep1
% - new_th1 is the true anomlay of r1 in kepT
% - new_th2 is the true anomlay of r2 in kepT
% - th2     is the true anomlay of r2 in kep2

    % Costant
    ALLERT = 6378 + 300; %Km to avoid collision with Earth
    DIM = 1e3; % to improve the precision of formula
    toll = 1e-2; %precision of 1 m
    PREC = 5e3;
    COST = 10;
    
    % Orbital value of first orbit
    a1 = kep1(1);
    e1 = kep1(2);
    i1 = kep1(3);
    OM1 = kep1(4);
    om1 = kep1(5);
    p1 = a1*(1-e1^2);
    
    % Orbital value of second orbit
    a2 = kep2(1);
    e2 = kep2(2);
    i2 = kep2(3);
    OM2 = kep2(4);
    om2 = kep2(5);
    p2 = a2*(1-e2^2);

    % Gravitational parameter [km^3/s^2]
    mu = 398600;
    
    % Define normal vector of the plane that contain the new orbit
    N = cross(r1,r2);
    if(N(3)<0)
        N = cross(r2,r1);
    end
    
    % Inclination of the plane
    iT = rad2deg(acos(N(end)/norm(N)));
    
    % Node line
    k = [0 0 1];
    n = cross(k,N);

    %RAAN (Ascensione retta del nodo ascendente)
    if n(2)>=0
        OMT=rad2deg(acos(n(1)/norm(n)));
    else
        OMT=rad2deg(2*pi-acos(n(1)/norm(n)));
    end
    
    % Find th1
    th1 = rad2deg(acos(1/e1*(p1/norm(r1)-1)));
    [new_r1,v1] = kep2car(a1,e1,i1,OM1,om1,th1,mu);
    if(norm(new_r1 - r1) > toll) % to avoid problem of acos that return from 0 to 180 
        th1 = 360 - th1;
        [~,v1] = kep2car(a1,e1,i1,OM1,om1,th1,mu);
    end
    
    % Find th2
    th2 = rad2deg(acos(1/e2*(p2/norm(r2)-1))); 
    [new_r2,v2] = kep2car(a2,e2,i2,OM2,om2, th2,mu);
    if(norm(new_r2 - r2) > toll) % to avoid problem of acos that return from 0 to 180 
        th2 = 360 - th2;
        [~,v2] = kep2car(a2,e2,i2,OM2,om2,th2,mu);
    end
    
    % Calcolate r1 and r2 in the orbital plane. 
    R = rotationMatrix(iT,OMT,0); % define the R matrix
    p1 = R*r1/DIM; % downgrade of DIM time. The most important thing is to calculate the angle
    p2 = R*r2/DIM;
    
    % Define the geometric locus of the positions of the second focus (form
    % of conic section usually hyperbola). p1 and p2 are the new foci of
    % the new conic
    a = 1/2*abs(norm(p2) - norm(p1));
    [~, ~, a, b, pos_center, th_rotation] = def_conic(p1,p2,a);

    % Extrapolate the points that belong to the hyperbola

    x1 = linspace(abs(a), COST*abs(a), PREC);  %From 0 to a, Upper
    ytop =  b*sqrt((x1).^2/a^2-1);  %Upper Hyperbola Part
    ybot = -b*sqrt((x1).^2/a^2-1);  %Lower Hyperbola Part
    
    % define the rotation matrix for clockwise angle (th_rotation)
    R = @(th) [cos(th), -sin(th); sin(th), cos(th)];
    
    % create the hyperbola by joining the x and y points
    x_point = [x1, x1, -x1, -x1];
    y_point = [ytop, ybot, ytop, ybot];
    
    % rotate hyperbola of th_rotation
    focal = R(deg2rad(th_rotation))*[x_point; y_point];

    % translate the hyperbola to the position of the center
    focal(1,:) = focal(1,:) + pos_center(1);
    focal(2,:) = focal(2,:) + pos_center(2);
    
    % define vector to find minimum  
    orbit_shape = [];
    dv = [];
    th = [];
    
    % for each point in focal
    for k = 1:length(focal)

        fc = focal(:, k);
        % find a 
        a_1 = 1/2* (norm(p1) - sqrt((p1(1) - fc(1))^2 + (p1(2) - fc(2))^2));
        a_2 = 1/2* (norm(p2) - sqrt((p2(1) - fc(1))^2 + (p2(2) - fc(2))^2));
        if(abs(a_1 - a_2) <= toll)
            a = a_1;
        else
            a = 1/2* (norm(p1) + sqrt((p1(1) - fc(1))^2 + (p1(2) - fc(2))^2));
        end

        % use the geometrical formula
        c = norm(fc/2);
        e = c/a;

        if(e>=0 && e<1) %considero solo le ellissi 
            
            a = a*DIM; % upscale a to the real value

            % check if the pericenter is above the risk of collision
            if(a*(1-e) >= ALLERT)
                
                % find omT in the transfert plane
                omT = rad2deg(atan2(fc(2),fc(1))) - 180; % non è detto che vali sempre perchè si potrebbe girare il fuoco!!!
                
                if(omT < 0)
                    omT = omT + 360;
                end
                
                % find new_th1
                new_th1 = rad2deg(atan2(p1(2),p1(1)))-omT;
                if(new_th1 < 0)
                    new_th1 = new_th1 + 360;
                end
                
                % find new_th2
                new_th2 = rad2deg(atan2(p2(2),p2(1)))-omT;
                if(new_th2 < 0)
                    new_th2 = new_th2 + 360;
                end
                
                [new_r1,new_v1] = kep2car(a,e,iT,OMT,omT,new_th1,mu);       
                [new_r2,new_v2] = kep2car(a,e,iT,OMT,omT,new_th2,mu);
                
                
                % plotting same example in direttissima
                %{
                if(mod(k,100) == 0)
                    kepT = [a,e,iT,OMT,omT,0];
                    [X,Y,Z] = plotOrbit(kepT,mu,360,0.1);
                    plot3(X,Y,Z);
                end
                %}


                % check the precisione with vectorial difference of r1 and new_r1 (same with r2 and new_r2)
                COND = norm(new_r2 - r2) < toll && norm(new_r1 - r1) < toll;
                if(COND)
                    % save all the possible elliptical transfert orbit in a vector and then valutate the dv min
                    orbit_shape = [orbit_shape; a e omT];
                    dv1 = norm(new_v1 - v1);
                    dv2 = norm(v2 - new_v2);
                    dv = [dv; dv1 + dv2];
                    th = [th; new_th1 new_th2];
                end
            end
        end
    end

    if(isempty(dv)) % don't exist any orbit that is an ellipse with the correct precision of 1 m without crashing the Earth
        kepT = [0, 0, 0, 0, 0, 0];
        dv = 0;
        th = [th1, 0, 0, th2];
    else
        % search minimum 
        [dv, k_vmin] = min(dv);
        
        % return the data of the minimum orbit
        kepT = [orbit_shape(k_vmin,1), orbit_shape(k_vmin,2),iT, OMT, orbit_shape(k_vmin,3), th(k_vmin,1)];
        th = [th1, th(k_vmin,1), th(k_vmin,2), th2];
    end