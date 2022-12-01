function [dv, omf, thf, dt] = changePeriapsisArg(a, e, omi, dom, th0)
% NB: forse mancano delle condizioni per evitare casini nel codice

% changePeriapsisArg.m - Maneuver to change the orientation in the orbital 
% plane at the closest point (change of pericentre anomaly)
%
% PROTOTYPE:
%   [dv, omf, thf, dt] = changePeriapsisArg(ai, ei, omi, dom, th0)
%
% DESCRIPTION:
% Function to calculate the speed required to change pericentre anomaly, given the 
% initial pericentre anomaly and the variation required.
% It also calculates the maneuver point and the time of flight from the initial
% position of the satellite to the maneuver.
% Angles in degrees.
% 
% INPUT: 
% a         [1x1] Semi-major axis                    [km]
% e         [1x1] Eccentricity                       [-]
% omi       [1x1] Initial pericentre anomaly         [°]
% dom       [1x1] Variation of pericentre anomaly    [°]
% th0       [1x1] Starting true anomaly              [°]
%
% OUTPUT:
% dv        [1x1] Velocity variation          [km/s]
% omf       [1x1] Final pericentre anomaly    [°]
% thf       [2x1] Final true anomaly vector   [°]
% dt        [1x1] Time from th0 to thf        [s]

% Gravitational Parameter [km^3/s^2]
mu = 398600.433;
% conversion deg to rad:
omi = deg2rad(omi);
dom = deg2rad(dom);
th0 = deg2rad(th0);
% omf in degrees:
omf = rad2deg(omi + dom);

% Since the possible intersections between the initial and final orbit are 
% in two different points, we group them in a single vector thf: 
% the first component is the true anomaly of the manouver point for the initial 
% orbit, the second component is the true anomaly for the final orbit. 
if th0 <= dom/2 || th0 > dom/2 + pi
    thf = [rad2deg(dom/2); rad2deg(2*pi - dom/2)];
else
    thf = [rad2deg(pi + dom/2); rad2deg(pi - dom/2)];
end

% Speed difference module of maneuver:
dv = 2*sqrt(mu/(a*(1-e^2)))*e*sin(dom/2);
% Time of Flight to reach the maneuver point:
dt = timeOfFlight(a, e, rad2deg(th0), thf(1), mu);
