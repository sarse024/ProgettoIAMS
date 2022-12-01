function [X,Y,Z] = plotOrbit(kepEl, mu, deltaTh, stepTh)
% plotOrbit.m - Plot the arc length deltaTh of the orbit described by
% kepEl.
%
% PROTOTYPE:
% plotOrbit(kepEl, mu, deltaTh, stepTh)
%
% DESCRIPTION:
% Plot the arc length of the orbit described by a set of orbital
% elements for a specific arc length.
%
% INPUT:
% kepEl [1x6] orbital elements [km,deg]
% mu [1x1] gravitational parameter [km^3/s^2]
% deltaTh [1x1] arc length [deg]
% stepTh [1x1] arc length step [deg]
%
% OUTPUT:
% X [1xn] X position [km]
% Y [1xn] Y position [km]
% Z [1xn] Z position [km]
    %expliciting the elements of the kepEl vector
    a=kepEl(1);
    e=kepEl(2);
    i=kepEl(3);
    OM=kepEl(4);
    om=kepEl(5);
    th0=kepEl(6);
    %changing the angles in radians
    i=deg2rad(i);
    OM=deg2rad(OM);
    om=deg2rad(om);
    th0=deg2rad(th0);
    deltaTh=deg2rad(deltaTh);
    stepTh=deg2rad(stepTh);
    %finding the final theta value and creating the theta vector which
    %contains every theta value
    thf=th0+deltaTh;
    theta=th0:stepTh:thf;
    r=[];
    %going back to the angles in degrees
    i=rad2deg(i);
    OM=rad2deg(OM);
    om=rad2deg(om);
    theta=rad2deg(theta);
    for k=1:length(theta)
        [rk] = kep2car(a, e, i, OM, om, theta(k), mu);
        r=[r;rk'];      %matrix with all the radii
    end
    %the final values are the three columns of the r matrix
    X=r(:,1);
    Y=r(:,2);
    Z=r(:,3);