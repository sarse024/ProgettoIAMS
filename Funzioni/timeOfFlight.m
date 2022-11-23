function [deltaT]=timeOfFlight(a,e,thi,thf,mu)
% timeOfFlight: Time elapsed to go from one point to another. Same orbit
%
% PROTOTYPE:
% [deltaT] = timeOfFlight(a,e,thi,thf,mu)
%
% DESCRIPTION:
% This function calculates the time elapsed to go from one point 
% to another on the same orbit
% 
% INPUT:
% a [1x1] Semi-major axis [km]
% e [1x1] Eccentricity [-]
% th1 [1x1] True anomaly first point [rad]
% th2 [1x1] True anomaly second point[rad]
%
% OUTPUT:
% deltaT [1x1] Time elapsed [s]
    %changing angles from degrees to radians
    thi=deg2rad(thi);
    thf=deg2rad(thf);
    %Orbit period
    T=2*pi*sqrt(a^3/mu);
    %eccentric anomaly
    E1=2*atan(sqrt((1-e)/(1+e))*tan(thi/2));
    if E1<0
        E1=E1+2*pi;
    end
    E2=2*atan(sqrt((1-e)/(1+e))*tan(thf/2));
    if E2<0
        E2=E2+2*pi;
    end
    %medium anomaly
    M1=E1-e*sin(E1);
    M2=E2-e*sin(E2);
    %elapsed times
    ti=sqrt(a^3/mu)*M1;
    tf=sqrt(a^3/mu)*M2;
    %time of flight
    if thf >= thi
        deltaT=tf-ti;
    else
        deltaT=tf-ti+T;
    end