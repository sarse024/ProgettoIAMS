function  [delta_v, th_f, delta_t] = changeOrbitShape(aI, eI, omI, aF, eF, omF, th0, option, rb)
% changeOrbitShape.m - Function to calculate the speed and time required to change the shape of the orbit with tangent transfers. Angles in input in degree
%
% PROTOTYPE:
% Hohmann transfer orbit:
% [delta_v, delta_v_f, th_f, delta_t] = changeOrbitShape(aI, eI, omI, aF, eF, omF, th0, option)
%
% Biellitical transfer orbit:
% [delta_v, th_f, delta_t] = changeOrbitShape(aI, eI, omI, aF, eF, omF, th0, option, rb)
%
% DESCRIPTION:
% Function to calculate the speed and time required to change the orbit shape with tangent transfers. 
% Based on the inputs is possible whether to carry out a Homann transfer (2 pulses) or a biellitic transfer (3 pulses). 
% Through the input 'option' it's possible to choose the starting point of the maneuver (apocenter or pericenter).
% The function also calculates the times from the start position to the final point.
% Output are vector. For total value you need to sum in your code.
%
% INPUT:
% Hohmann transfer orbit and Biellitcal transfer
% aI        [1x1]       Semi-major axis of first orbit                      [km]
% eI        [1x1]       Eccentricity of first orbit                         [-]
% omI       [1x1]       Pericentre anomaly of first orbit                   [deg]
% aF        [1x1]       Semi-major axis of target orbit                     [km]
% eF        [1x1]       Eccentricity of target orbit                        [-]
% omF       [1x1]       Pericentre anomaly of target orbit                  [deg]
% th0       [1x1]       True anomaly on fist orbit                          [deg]
% option    [string]    Choose of the point of initial maneuver             [-]
%                       'apo' -> will start at Apocenter
%                       'per' -> will start at Pericenter
% Only Biellitical transfer orbit:
% rb       [1x1]        Apocenter radius of transfert orbit in biellitical  [Km]
%
% OUTPUT:
% Hohmann transfer orbit:
% delta_v   [2x1]       Speed ​​difference module in each point of maneuver   [km/s]
% th        [1x1]       True anomaly final and point of maneuver            [deg]
% delta_t   [2x1]       Time of flight from th0 to th_f for each transfer   [s]
%
% Biellitical transfer orbit:
% delta_v   [3x1]       Speed ​​difference module in each point of maneuver   [km/s]
% th        [1x1]       True anomaly final and point of maneuver            [deg]
% delta_t   [3x1]       Time of flight from th0 to th_f for each transfer   [s]

% Gravitational parameter [km^3/s^2]
global mu;

% Conversion deg to rad
omI = deg2rad(omI);
omF = deg2rad(omF);
th0 = deg2rad(th0);

% check if orbit are tangent
if(~(mod(omI,pi) == mod(omF,pi)))
   error("Orbit not tangent\n");
end

if (~strcmp(option,'per') && ~strcmp(option, 'apo'))
    error(['Maneuvering point not identified, enter:' ...
           '\napo : to start from apocenter' ...
           '\nper : to start from pericenter\n'])
end

if(nargin == 8)
    %%% Hohmann transfer orbit %%%
    fprintf('\nHomann trasfer strategy\n')


    % initialize output vector
    delta_v = zeros(2,1);
    delta_t = zeros(2,1);

    % Find apocenter and pericenter radius of the transfer orbit
    if(omF == omI)
        if(strcmp(option,'per'))% Case pericenter-apocenter
            r1 = aI*(1 - eI);
            r2 = aF*(1 + eF);
            th_f = pi;
        else
            r1 = aI*(1 + eI); % Case apocenter-pericenter
            r2 = aF*(1 - eF);
            th_f = 0;
        end
    else
        if(strcmp(option,'per'))% Case pericenter-pericenter
            r1 = aI*(1 - eI);
            r2 = aF*(1 - eF);
            th_f = 0;
        else
            r1 = aI*(1 + eI); % Case apocenter-apocenter
            r2 = aF*(1 + eF);
            th_f = pi;
        end
    end

    % Semi-major axis of transfer orbit
    a_t = (r1 + r2)/2;
    
    % Speed ​​difference module at first point of maneuver
    delta_v(1) = sqrt(2*mu *((1/r1) - 1/(2*a_t))) - sqrt(2*mu *((1/r1) - 1/(2*aI)));
    % Speed ​​difference module at second point of maneuver
    delta_v(2) = sqrt(2*mu*((1/r2) - 1/(2*aF))) - sqrt(2*mu *((1/r2) - 1/(2*a_t)));
    
    % Time of flight from th0 to first point of maneuver
    if(strcmp(option,'per'))
        delta_t(1) = timeOfFlight(aI,eI,th0, 0); % Case pericenter
    else
        delta_t(1) = timeOfFlight(aI,eI,th0, pi); % Case apocenter
    end

    delta_t(2) = pi*sqrt(a_t^3/mu);

else 
    %%% Biellitical transfer %%%
    fprintf('\nBiellitical trasfer strategy\n')

    %initialize output vector
    delta_v = zeros(3,1);
    delta_t = zeros(3,1);

    % Find apocenter and pericenter radius of the two transfer orbit
    if(omF == omI)
        if(strcmp(option,'per'))%caso pericenter-pericentro
            r1 = aI*(1 - eI);
            r3 = aF*(1 - eF);
            th_f = 0;
        else
            r1 = aI*(1 + eI); %caso apocenter-apocenter
            r3 = aF*(1 + eF);
            th_f = pi;
        end
    else
        if(strcmp(option,'per'))%caso pericenter-apocenter
            r1 = aI*(1 - eI);
            r3 = aF*(1 + eF);
            th_f = pi;
        else
            r1 = aI*(1 + eI); %caso apocenter-pericenter
            r3 = aF*(1 - eF);
            th_f = 0;
        end
    end
    % Semi-major axis of first transfer orbit
    a_t1 = (r1 + rb)/2;

    % Semi-major axis of second transfer orbit
    a_t2 = (rb + r3)/2;

    % Speed ​​difference module at first point of maneuver
    delta_v(1) = sqrt(2*mu *((1/r1) - 1/(2*a_t1))) - sqrt(2*mu *((1/r1) - 1/(2*aI)));

    % Speed ​​difference module at second point of maneuver
    delta_v(2) = sqrt(2*mu*((1/rb) - 1/(2*a_t2))) - sqrt(2*mu *((1/rb) - 1/(2*a_t1)));

    % Speed ​​difference module at third point of maneuver
    delta_v(3) = sqrt(2*mu*((1/r3) - 1/(2*aF))) - sqrt(2*mu *((1/r3) - 1/(2*a_t2)));

    % Time of flight from th0 to first point of maneuver
    if(strcmp(option,'per'))
        delta_t(1) = timeOfFlight(aI,eI,th0, 0); % Case pericenter
    else
        delta_t(1) = timeOfFlight(aI,eI,th0, pi); % Case apocenter
    end

    % Time of flight of first trasfer orbit
    delta_t(2) = pi*sqrt(a_t1^3/mu);

    % Time of flight of first trasfer orbit
    delta_t(3) = pi*sqrt(a_t2^3/mu);
end

% Conversion rad to deg
th_f = deg2rad(th_f);