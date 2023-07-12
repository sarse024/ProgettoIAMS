%% 2 IMPULSE STRATEGY "DIRETTISSIMA"
clear;
clc;
close all;

% costant
global mu;

mu = 398600;
precision = 0.1;
step_animation = 10;
t_tot = 0;

%%%%%%%% GRUPPO B7 %%%%%%%%
dati_elaborati = [5088.9118 -3196.5659 -8222.7989 1.9090 5.6220 -1.0700 14020.0000 0.3576 1.3220 0.9764 1.8130 0.4336];

% --- Initial Orbit ---
rI = dati_elaborati(1:3)';
vI = dati_elaborati(4:6)';

% Orbital parametre of initial orbit
[aI, eI, iI, OMI, omI, thI] = car2kep(rI, vI, mu);
kepI = [aI,eI,iI,OMI,omI,thI];

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

% Initialize main figure
fig1 = figure;
Terra3d;
hold on
grid on
title('"Direttissima"  Manoeuvre', fontsize=25)

% plot intial orbit
[X,Y,Z] = plotOrbit(kepI,mu,360,precision);
orbitI = plot3(X,Y,Z, '--', 'LineWidth', 1.5);
start = plot3(rI(1), rI(2), rI(3), 'xb', 'LineWidth', 4);

% plot final orbit
[X,Y,Z] = plotOrbit(kepF,mu,360,precision);
orbitF = plot3(X,Y,Z, '--', 'LineWidth', 1.5);
target = plot3(rF(1), rF(2), rF(3), 'xr', 'LineWidth', 4);

% find the best transfert orbit from rI to rF
[kepT, dv, th] = findOrbit(rI,rF,kepI,kepF);

% plot the full transfert orbit in dashed line
[X,Y,Z] = plotOrbit(kepT,mu,360, precision);
plot3(X,Y,Z, '--r','LineWidth',1);

% plot trajectory in continuos line
[X_traj,Y_traj,Z_traj] = plotOrbit(kepT, mu, mod(th(3)-th(2),360),precision);
orbitT= plot3(X_traj,Y_traj,Z_traj,'r','LineWidth',2);

% Calculate the two speeds by vector difference
[r1,v1] = kep2car(kepT(1),kepT(2),kepT(3),kepT(4),kepT(5),th(2),mu);
[r2,v2] = kep2car(kepT(1),kepT(2),kepT(3),kepT(4),kepT(5),th(3),mu);

% time of flight from r1 to r2 on kepT
t = timeOfFlight(kepT(1),kepT(2),th(2),th(3),mu);

%%%% OUTPUT %%%%
fprintf('IMPLEMENTAZIONE STRATEGIA 2 IMPULSI DIRETTA:\n')

fprintf('\n--- Prima manovra ---\n')
stampInfoManovra(kepT,th(1),0,norm(v1-vI),t_tot)
fprintf('Errore sul primo punto di manovra: %2.2f m\n', norm(r1-rI)*1000)

fprintf('\n--- Seconda manovra ---\n')
t_tot = t;
stampInfoManovra(kepF,th(3),t,norm(vF-v2),t_tot)
fprintf('Errore sul secondo punto di manovra: %2.2f m\n', norm(r2-rF)*1000)

% Summary output
fprintf('\n---- Riassunto strategia 2 impulsi diretta ----\n')
fprintf('Variazione di velocità totale: %2.3f m/s\n', norm(dv)*1000)
fprintf('Tempo totale di manovra:\n')
fprintf('Secondi: %5.2f s\n', t)
fprintf('Minuti: %5.2f m\n', t/60)
fprintf('Ore: %5.2f h\n', t/60/60)
fprintf('Giorni: %5.2f d\n', t/60/60/24)

%%% ANIMATION %%%
% satellite
h = plot3(nan,nan,nan,"om", 'MarkerSize',10, 'MarkerFaceColor', 'm');

legend([orbitI, start, orbitF, target, orbitT, h], 'Initial Orbit', 'Starting point', 'Final Orbit', 'Target point', 'Transfer Trajectory', 'Satellite', 'Location','southeast')
xlabel('X axis [Km]'), ylabel('Y axis [Km]'), zlabel('Z axis [Km]')
set(gca, "CameraPosition", 10^4*[1 1.4 0.1]);
%%
gifFile = 'Direttissima.gif';
exportgraphics(fig1, gifFile);

% draw satellite on the trajectory
for i = 1:step_animation:length(X_traj)
    set(h,'XData',X_traj(i),'YData',Y_traj(i),'ZData',Z_traj(i));
    drawnow
    exportgraphics(fig1, gifFile, Append=true);
end
set(h,'XData',X_traj(end),'YData',Y_traj(end),'ZData',Z_traj(end));
drawnow


