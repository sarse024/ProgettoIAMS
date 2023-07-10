%% STANDARD STRATEGY
clear;
clc;
close all;

% Costant
global mu;

mu = 398600;

option = 'per'; % useful for changing start point in changeShape
precision = 0.1; % step precision in plotOrbit [deg]
step_animation = 50; % variable for animation of trajectory

% Initialize main figure
figure
Terra3d
hold on
grid on
title('Strategia base')

%%%%%%%% GRUPPO B7 %%%%%%%%
dati_elaborati = [5088.9118 -3196.5659 -8222.7989 1.9090 5.6220 -1.0700 14020.0000 0.3576 1.3220 0.9764 1.8130 0.4336];

%altri gruppi per prova da cancellare
%dati_elaborati = [-7394.8822 849.6501 4181.2069 -3.0970 -5.2110 -3.5690 15660.0000 0.2659 0.8321 1.2440 2.9320 3.0970];
%dati_elaborati = [-6096.8804 -1361.2133 4888.2313 -0.2075 -7.2520 -1.4010 12860.0000 0.2842 1.4520 0.3340 2.7340 2.9360];

% --- Initial Orbit ---
rI = dati_elaborati(1:3)';
vI = dati_elaborati(4:6)';

% Orbital parametre of initial orbit
[aI, eI, iI, OMI, omI, thI] = car2kep(rI, vI, mu);
kepI = [aI,eI,iI,OMI,omI,thI];

% plot initial orbit
[X,Y,Z] = plotOrbit(kepI,mu,360,precision);
orbitI = plot3(X,Y,Z, 'LineWidth', 2);
% plot starting point
start = plot3(rI(1), rI(2), rI(3), 'xb', 'LineWidth', 4);

% --- Final Orbit ---
aF = dati_elaborati(7);
eF = dati_elaborati(8);
iF = rad2deg(dati_elaborati(9));
OMF = rad2deg(dati_elaborati(10));
omF = rad2deg(dati_elaborati(11));
thF = rad2deg(dati_elaborati(12));

% Orbital parametre of final orbit
[rF, vF] = kep2car(aF,eF,iF,OMF,omF,thF,mu);
kepF = [aF,eF,iF,OMF,omF,thF];

% plot final orbit
[X,Y,Z] = plotOrbit(kepF,mu,360,precision);
orbitF = plot3(X,Y,Z, 'LineWidth', 2);
target = plot3(rF(1), rF(2), rF(3), 'xr', 'LineWidth', 4);

% vector for store point of the trajectory
X_traj = [];
Y_traj = [];
Z_traj = [];



fprintf('IMPLEMENTAZIONE STRATEGIA BASE:\n\n')

fprintf(['Parametri orbitali punto di partenza:\n' ...
        '\ta\t\t\te\t\ti\t\tOM\t\tom\t\tth\n ' ...
        '%4.2f\t %.4f\t %2.2f\t %2.2f\t %2.2f\t %2.2f\n\n'], kepI(1), kepI(2), kepI(3),kepI(4),kepI(5), kepI(6))
fprintf(['Parametri orbitali punto di arrivo:\n' ...
        '\ta\t\t\te\t\ti\t\tOM\t\tom\t\tth\n ' ...
        '%4.2f\t %.4f\t %2.2f\t %2.2f\t %2.2f\t %2.2f\n'], kepF(1), kepF(2), kepF(3),kepF(4),kepF(5), kepF(6))

% --- 1° Manoeuvre: change plane --- 
[dv1, om1, th1, dt1] = changeOrbitalPlane(aI, eI, iI, OMI, omI, iF, OMF, 200); %forse problema dt

% update trajectory from thI to point of maneuvre th1
[X,Y,Z] = plotOrbit(kepI,mu,abs(thI-th1),precision);
X_traj = X;
Y_traj = Y;
Z_traj = Z;

kep1 = [aI, eI, iF, OMF, om1, th1]; % Orbital parametre vector of first orbit transfert

% Plot first orbit transfert
[X,Y,Z] = plotOrbit(kep1,mu,360,precision);
orbit1 = plot3(X,Y,Z,'--r', 'LineWidth', 1);

% Plot first transfert manouvre point
[r1, v1] = kep2car(aI, eI, iF, OMF, om1, th1, mu);
pt_changePlane = plot3(r1(1), r1(2), r1(3), 'ok', 'LineWidth', 2);

% Output change plane
fprintf('\n---- CAMBIO DI PIANO ----\n')
stampInfoManovra(kep1, thI, dt1,dv1)