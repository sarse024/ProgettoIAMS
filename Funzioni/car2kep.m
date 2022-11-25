function [a,e,i,OM,om,th]=car2kep(r,v,mu)
% car2kep.m - Conversion from Cartesian coordinates to Keplerian elements
%
% PROTOTYPE:
% [a, e, i, OM, om, th] = car2kep(r, v, mu)
%
% DESCRIPTION:
% Conversion from Cartesian coordinates to Keplerian elements. Angles in
% degrees.
%
% INPUT:
% r [3x1] Position vector [km]
% v [3x1] Velocity vector [km/s]
% mu [1x1] Gravitational parameter [km^3/s^2]
%
% OUTPUT:
% a [1x1] Semi-major axis [km]
% e [1x1] Eccentricity [-]
% i [1x1] Inclination [deg]
% OM [1x1] RAAN [deg]
% om [1x1] Pericentre anomaly [deg]
% th [1x1] True anomaly [deg]
    r_norm=norm(r); %norm position vector
    v_norm=norm(v); %norm velocity vector
    %angular momentum and his vector
    h_vect=cross(r,v);
    h_norm=norm(h_vect);
    %inclination
    i=acos(h_vect(end)/h_norm);
    %eccentricity vector and eccentricity
    e_vect=(1/mu)*((v_norm^2-mu/r_norm)*r-dot(r,v)*v);
    e=norm(e_vect);
    %specific mechanic energy and major semiaxis
    eps=0.5*v_norm^2-mu/r_norm;
    a=-mu/(2*eps);
    %nodes line and his module
    k=[0 0 1];
    n_vect=cross(k,h_vect);
    n_norm=norm(n_vect);
    %RAAN (Ascensione retta del nodo ascendente)
    if n_vect(2)>=0
        OM=acos(n_vect(1)/n_norm);
    else
        OM=2*pi-acos(n_vect(1)/n_norm);
    end
    %pericentre anomaly
    if e_vect(end)>=0
        om=acos(dot(n_vect,e_vect)/(n_norm*e));
    else
        om=2*pi-acos(dot(n_vect,e_vect)/(n_norm*e));
    end
    %radius velocity
    vr=dot(r,v)/r_norm;
    %true anomaly
    if vr>=0
        th=acos(dot(e_vect,r)/(e*r_norm));
    else
        th=2*pi-acos(dot(e_vect,r)/(e*r_norm));
    end
    %change angles in degrees
    i=rad2deg(i);
    OM=rad2deg(OM);
    om=rad2deg(om);
    th=rad2deg(th);
