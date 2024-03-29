function [delta_v, om2, th, delta_t] = changeOrbitalPlane(a, e, i1, OM1, om1, i2, OM2, th0)
% changeOrbitalPlane.m - Orbital plane change maneuver at the closest point. Angles in input in degree
%
% PROTOTYPE:
% [delta_v, om2, th, delta_t] = changeOrbitalPlane(a, e, i1, OM1, om1, i2, OM2, th0)
%
% DESCRIPTION:
% Function to calculate the speed required to make a change
% plane at the closest intersection point. Then it find other parameters modified during the maneuver (om and th).
% The function also calculates the times flight from the position of the
% satellite to the maneuver.
%
% INPUT:
% a         [1x1]   Semi-major axis                             [km]
% e         [1x1]   Eccentricity                                [-]
% i1        [1x1]   Inclination of first orbit                  [deg]
% OM1       [1x1]   RAAN of first orbit                         [deg]
% om1       [1x1]   Pericentre anomaly of first orbit           [deg]
% i2        [1x1]   Inclination of final orbit                  [deg]
% OM2       [1x1]   RAAN of final orbit                         [deg]
% th0       [1x1]   True anomaly                                [deg]
%
% OUTPUT:
% delta_v   [1x1]   Speed ​​difference module                     [km/s]
% om2       [1x1]   Pericentre anomaly of final orbit           [deg]
% th        [1x1]   True anomaly final and point of maneuver    [deg]
% delta_t   [1x1]   Time of flight from th0 to th               [s]

% Gravitational parameter [km^3/s^2]
global mu;

% Tollerance 
toll = 1e-5;

% Conversion deg to rad
i1 = deg2rad(i1);
OM1 = deg2rad(OM1);
om1 = deg2rad(om1);
i2 = deg2rad(i2);
OM2 = deg2rad(OM2);
th0 = deg2rad(th0);

dOM = OM2-OM1;

if(dOM == 0)
    % find th of descending node
    th = pi-om1;
    
    % check if satellite posizion in already pass the descending node
    if(th0 > th)
        % find th of ascending node
        th = 2*pi - om1;
    end

    % find transverse velocity in th
    v_t = sqrt(mu/(a*(1-e^2)))*(1+e*cos(th));

    % Speed ​​difference module of maneuver  
    delta_v = abs(2 * v_t * sin((i2-i1)/2));
    
    % Final pericentre anomaly
    om2 = om1;
else
    % To avoid error in spherical triangle definition
    if(i1 == 0)
        i1 = i1+toll; %introduce little error (some metres)
    end
    if(i2 == 0)
        i2 = i2+toll;
    end
    
    % Find dihedral angle
    alpha = acos(cos(i1)*cos(i2) + sin(i1)*sin(i2)*cos(dOM));
    
    % Check if it's better a biellitical strategies
    if(rad2deg(alpha) >= 38.4)
      fprintf("\nRecommend plan change transfer with bielliptical manoeuvre\n")
    end

    %find sides u1 and u2 of spherical triangle
    cosU1 = (cos(alpha)*cos(i1) - cos(i2))/(sin(alpha)*sin(i1));
    sinU1 = sin(i2)* sin(dOM)/sin(alpha);

    cosU2 = (cos(i1) - cos(alpha)*cos(i2))/(sin(alpha)*sin(i2));
    sinU2 = sin(i1)* sin(dOM)/sin(alpha);

    u1 = atan2(sinU1, cosU1); %SPEZZARE IN DUE CASI
    u2 = atan2(sinU2, cosU2);

    % Find point of maneuver (th)
    if(u1 >= om1)
        th = u1 - om1;
    else
        th = 2*pi + u1 - om1;
    end
    
    % Final pericentre anomaly (om2)
    if(u2 >= th)
        om2 = u2 - th;
    else
        om2 = 2*pi + (u2 - th);
    end
    
    delta_t_ascendent_node = timeOfFlight(a,e,rad2deg(th0), rad2deg(th),mu);
    delta_t_discendent_node = timeOfFlight(a,e,rad2deg(th0), rad2deg(mod(pi+th,2*pi)) ,mu);

    if(delta_t_discendent_node <= delta_t_ascendent_node)
        th = mod(pi+th,2*pi);
    end

    % Find transverse velocity in th
    v_t = sqrt(mu/(a*(1-e^2)))*(1+e*cos(th));

    % Speed ​​difference module of maneuver  
    delta_v = abs(2*v_t*sin(alpha/2)); 
end

% Conversion output rad to deg
th = rad2deg(th);
om2 = rad2deg(om2);

% Time of flight from th0 to th
delta_t = timeOfFlight(a,e,rad2deg(th0),th,mu);


end