% Circular Bielliptic Transfer strategy (with optimal radius)
clear;
clc;
close all;
% defining parameters
global mu;
mu = 398600;
precision = 0.1;
step_animation = 10;

%%%%%%%% GRUPPO B7 %%%%%%%%
dati_elaborati = [5088.9118 -3196.5659 -8222.7989 1.9090 5.6220 -1.0700 14020.0000 0.3576 1.3220 0.9764 1.8130 0.4336];

% Initial Orbit data
rI = dati_elaborati(1:3)';
vI = dati_elaborati(4:6)';

% Other parameters of the initial orbit:
[aI, eI, iI, OMI, omI, thI] = car2kep(rI, vI, mu);
kepI = [aI,eI,iI,OMI,omI,thI];

% Final Orbit Data:
aF = dati_elaborati(7);
eF = dati_elaborati(8);
iF = rad2deg(dati_elaborati(9));
OMF = rad2deg(dati_elaborati(10));
omF = rad2deg(dati_elaborati(11));
thF = rad2deg(dati_elaborati(12));

% Other parameters of the final orbit:
[rF, vF] = kep2car(aF,eF,iF,OMF,omF,thF,mu);
kepF = [aF,eF,iF,OMF,omF,thF];

% defining figure:
fig = figure;
Terra3d
hold on
grid on
title('Circular Bielliptic Strategy', fontsize=25)

% vector for store point of the trajectory
X_traj = [];
Y_traj = [];
Z_traj = [];

% Initial Orbit
[X,Y,Z] = plotOrbit(kepI,mu,360,precision);
orbitI = plot3(X,Y,Z, '--','LineWidth', 1);
start = plot3(rI(1), rI(2), rI(3), 'xb', 'LineWidth', 4);

% Final Orbit
[X,Y,Z] = plotOrbit(kepF,mu,360,precision);
orbitF = plot3(X,Y,Z, '--', 'LineWidth', 1);
target = plot3(rF(1), rF(2), rF(3), 'xr', 'LineWidth', 4);

fprintf('IMPLEMENTAZIONE STRATEGIA BIELLITTICA OTTIMALE:\n\n')
fprintf(['Parametri orbitali punto di partenza:\n' ...
        '\ta\t\t\te\t\ti\t\tOM\t\tom\t\tth\n ' ...
        '%4.2f\t %.4f\t %2.2f\t %2.2f\t %2.2f\t %2.2f\n\n'], kepI(1), kepI(2), kepI(3),kepI(4),kepI(5), kepI(6))
fprintf(['Parametri orbitali punto di arrivo:\n' ...
        '\ta\t\t\te\t\ti\t\tOM\t\tom\t\tth\n ' ...
        '%4.2f\t %.4f\t %2.2f\t %2.2f\t %2.2f\t %2.2f\n'], kepF(1), kepF(2), kepF(3),kepF(4),kepF(5), kepF(6))

% Computing the angle of orbital plane change in relation to the initial
% orbit: this is done to give the impulse to start the bielliptic strategy
% 180° before this angle. This allows the manuever of orbital plane change
% to be done at the apocentre of the elliptic orbit.
[~, ~, th_cp, ~] = changeOrbitalPlane(aI, eI, iI, OMI, omI, iF, OMF, thI);
% Starting point of the bielliptic:
th_in_be = th_cp - 180;

% First maneuver: Circularization of the initial orbit (at the apocentre)
raI = aI*(1+eI);
a1 = raI;
e1 = 0;
[dv1, th1, dt1] = changeOrbitShape(aI, eI, omI, a1, e1, omI, thI, 'apo');
thman1 = 180;
dth1 = thman1 - thI;
t_tot1 = sum(dt1);

% update trajectory from thI to point of maneuvre th1
[X,Y,Z] = plotOrbit(kepI,mu, dth1,precision);
orbitI = plot3(X,Y,Z, 'b', 'LineWidth', 2);
X_traj = X;
Y_traj = Y;
Z_traj = Z;

% plotting point of circularization
[r1, v1] = kep2car(aI, eI, iI, OMI, omI, thman1, mu);
pt_circ1 = plot3(r1(1), r1(2), r1(3), 'oc', 'MarkerSize', 5, 'MarkerFaceColor', 'c');
kep1 = [a1, e1, iI, OMI, omI, thman1];



% Output and plot:
%[X,Y,Z] = plotOrbit(kep1,mu,360,precision);
%orbit1 = plot3(X,Y,Z,'--r', 'LineWidth', 1);
fprintf('\n---- CIRCOLARIZZAZIONE 1 ----\n')
stampInfoManovra(kep1,thman1, sum(abs(thman1)), sum(abs(dv1)), t_tot1)

% Defining parameters for the first half of the bielliptic
rb = 1.903405340534054e+04; % optimal bielliptic radius computed
rp_t1 = a1; % rp_t1 = r1
ra_t1 = rb;
a_t1 = (rp_t1+ra_t1)/2;
e_t1 = (ra_t1-rp_t1)/(ra_t1+rp_t1);
om_t1 = omI + th_in_be;

% Second maneuver: First half of the bielliptic transfer
[dv2, th2, dt2] = changeOrbitShape(a1, e1, om_t1, a_t1, e_t1, om_t1, th_in_be, 'per');
t_tot2 = sum(dt1) + sum(dt2);

% update trajectory from th1 to point of maneuvre th2
[X,Y,Z] = plotOrbit(kep1,mu, abs(th_in_be+180) ,precision);
orbit1 = plot3(X,Y,Z, 'c', 'LineWidth', 2);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];


% Plot and output of the second maneuver:
kep_t1 = [a_t1, e_t1, iI, OMI, om_t1, 0];
%[X,Y,Z] = plotOrbit(kep_t1,mu,360,precision);
%orbit2 = plot3(X,Y,Z,'--m', 'LineWidth', 1);
fprintf('\n---- PRIMO IMPULSO BIELLITTICA ----\n')
stampInfoManovra(kep_t1,th_in_be, sum(abs(dt2)), sum(abs(dv2)), t_tot2)

% Plot starting point of the bielliptic:
[r_in_be, v_in_be] = kep2car(a_t1, e_t1, iI, OMI, om_t1, th1, mu);
in_be = plot3(r_in_be(1), r_in_be(2), r_in_be(3), 'og','MarkerSize', 5, 'MarkerFaceColor', 'g');

% Third maneuver: change of plane
[dv3, om3, th3, dt3] = changeOrbitalPlane(a_t1, e_t1, iI, OMI, om_t1, iF, OMF, (th2-1e-7)); 
t_tot3 = sum(dt1) + sum(dt2) + dt3;
kep2 = [a_t1, e_t1, iF, OMF, om3, th3]; 

% update trajectory from th2 to point of maneuvre th3
[X,Y,Z] = plotOrbit(kep_t1,mu, abs(th3) ,precision);
orbit2 = plot3(X,Y,Z, 'g', 'LineWidth', 2);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];


% plotting point of plane change
[r2, v2] = kep2car(a_t1, e_t1, iF, OMF, om3, th3, mu);
pt_changePlane = plot3(r2(1), r2(2), r2(3), 'or', 'MarkerSize', 5, 'MarkerFaceColor', 'r');

% Output change of plane
fprintf('\n---- CAMBIO DI PIANO ----\n')
stampInfoManovra(kep2, th2, dt3,dv3, t_tot3)

% Fourth maneuver: second half of the bielliptic
ra3 = rb;
rp3 = aF*(1+eF);
a_t2 = (ra3 + rp3)/2;
e_t2 = (ra3-rp3)/(ra3+rp3);

[dv4, th4, dt4] = changeOrbitShape(a_t1, e_t1, om3, a_t2, e_t2, om3, (th3-1e-7), 'apo');
t_tot4 = sum(dt1) + sum(dt2) + dt3 + sum(dt4);


% Plotting second half of the bielliptic with new plane and output:
kep_t2 = [a_t2, e_t2, iF, OMF, om3, th3];

%[X,Y,Z] = plotOrbit(kep_t2,mu,360,precision);
%orbit3 = plot3(X,Y,Z,'--b', 'LineWidth', 1);

fprintf('\n---- SECONDO IMPULSO BIELLITTICA ----\n')
stampInfoManovra(kep_t2,th3, sum(abs(dt4)), sum(abs(dv4)), t_tot4)

% Fifth maneuver: ricircularization of the orbit to intersect the final
% orbit in its apocentre
a5 = aF*(1+eF); % raggio dell'orbita circolare = raggio dell'apocentro finale
e5 = 0;
[dv5, th5, dt5] = changeOrbitShape(a_t2, e_t2, om3, a5, e5, om3, th4, 'per');
t_tot5 = sum(dt1) + sum(dt2) + dt3 + sum(dt4) + sum(dt5);
kep5 = [a5, e5, iF, OMF, om3, th5];
% Plotting point of circularization:
%[r3, v3] = kep2car(a_t2, e_t2, iF, OMF, om3, th4, mu);
%pt_circ21 = plot3(r3(1), r3(2), r3(3), 'o', 'MarkerSize', 5, 'MarkerFaceColor', '#EDB120');

% update trajectory from th2 to point of maneuvre th3
[X,Y,Z] = plotOrbit(kep_t2,mu, abs(360 - abs(omF-om3)) ,precision);
orbit3 = plot3(X,Y,Z, 'r', 'LineWidth', 2);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

% Plot final circular orbit
%[X,Y,Z] = plotOrbit(kep5,mu,360,precision);
%orbit4 = plot3(X,Y,Z,'--k', 'LineWidth', 1);

% Output second circularization
fprintf('\n---- SECONDA CIRCOLARIZZAZIONE ----\n')
stampInfoManovra(kep5,th4, sum(dt5), sum(abs(dv5)), t_tot5)

% Sixth maneuver: decircularization into the final orbit
[dv6, th6, dt6] = changeOrbitShape(a5, e5, omF, aF, eF, omF, th5, 'apo');
t_tot6 = sum(dt1) + sum(dt2) + dt3 + sum(dt4) + sum(dt5) + sum(dt6);
kep6 = [aF, eF, iF, OMF, omF, 180];

% Plot point of decircularization
[r5, v5] = kep2car(aF, eF, iF, OMF, omF, th5, mu);
pt_circ22 = plot3(r5(1), r5(2), r5(3), 'o', 'Color', '#EDB120','MarkerSize', 5, 'MarkerFaceColor', '#EDB120');

% update trajectory from th5 to point of maneuvre thF
[X,Y,Z] = plotOrbit(kep6, mu, abs(180+thF) ,precision);
orbitF = plot3(X,Y,Z, 'Color', "#EDB120", 'LineWidth', 2);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

% Output decircularization
fprintf('\n---- DECIRCOLARIZZAZIONE ----\n')
stampInfoManovra(kep6,th5, sum(dt6), sum(abs(dv6)), t_tot6)

% Ultimately, we reach the target point
dt7 = timeOfFlight(aF, eF, th6, thF,mu);

if(thF < th6)
    dTh = 360 + thF - th6;
else
    dTh = thF - th6; 
end

% Total cost and transfer time:
dt_tot = sum(dt1) + sum(dt2) + dt3 + sum(dt4) + sum(dt5) + sum(dt6) + dt7;
dv_tot = sum(abs(dv1)) + sum(abs(dv2)) + abs(dv3) + sum(abs(dv4)) + sum(abs(dv5)) + sum(abs(dv6));

% initialize satellite and camera
h = plot3(nan,nan,nan,"om", 'MarkerSize',10, 'MarkerFaceColor', 'm');
set(gca, "CameraPosition", 10^4*[30 8 10]);

fprintf('\nAttesa fino al punto finale: %4.2f s\n', dt_tot)
% legend
legend([orbitI start pt_circ1 orbit1 in_be orbit2 pt_changePlane orbit3 pt_circ22 orbitF target h], 'Initial orbit', 'Initial point', 'First circularization point', 'Circularized initial orbit', 'Beginning of bi-elliptic transfer point', 'First part of bi-elliptic transfer', 'Plane change and shape change point', 'Second part of bi-elliptic transfer', 'Decircularization point', 'Final orbit', 'Final point', 'Satellite')
% Final Output with results:
fprintf('\n---- Riassunto strategia biellittica Circolarizzata ----\n')
fprintf('Variazione di velocità totale: %2.3f km/s\n', dv_tot)
fprintf('Tempo totale di manovra:\n')
fprintf('Secondi: %5.2f s\n', dt_tot)
fprintf('Minuti: %5.2f min\n', dt_tot/60)
fprintf('Ore: %5.2f h\n', dt_tot/60/60)
fprintf('Giorni: %5.2f d\n', dt_tot/60/60/24)

%% ANIMATION

xlabel('X axis [Km]'), ylabel('Y axis [Km]'), zlabel('Z axis [Km]')

%gifFile = 'Biellittica.gif';
%exportgraphics(fig1, gifFile);

for i = 1:step_animation:length(X_traj)
    set(h,'XData',X_traj(i),'YData',Y_traj(i),'ZData',Z_traj(i));
    drawnow
    %exportgraphics(fig1, gifFile, Append=true);
end
set(h,'XData',X_traj(end),'YData',Y_traj(end),'ZData',Z_traj(end));
drawnow
