%% Prima strategia alternativa: cambio forma e poi cambio piano
clear;
clc;
close all;

global mu;

mu = 398600;

option = 'per';
precision = 0.1;
step_animation = 50;
t_tot = 0;

%initialize main figure
fig1 = figure;
Terra3d
hold on
grid on
title('Alternative standard strategy', fontsize=25)

%%%%%%%% GRUPPO B7 %%%%%%%%
dati_elaborati = [5088.9118 -3196.5659 -8222.7989 1.9090 5.6220 -1.0700 14020.0000 0.3576 1.3220 0.9764 1.8130 0.4336];

% --- Initial Orbit ---
rI = dati_elaborati(1:3)';
vI = dati_elaborati(4:6)';

% Orbital parameters of initial orbit
[aI, eI, iI, OMI, omI, thI] = car2kep(rI, vI, mu);
kepI = [aI,eI,iI,OMI,omI,thI];

% plot initial orbit
[X,Y,Z] = plotOrbit(kepI,mu,360,precision);
orbitI = plot3(X,Y,Z, '--', 'LineWidth', 1);
% plot starting point
start = plot3(rI(1), rI(2), rI(3), 'xb', 'LineWidth', 4);

% --- Final Orbit ---
aF = dati_elaborati(7);
eF = dati_elaborati(8);
iF = rad2deg(dati_elaborati(9));
OMF = rad2deg(dati_elaborati(10));
omF = rad2deg(dati_elaborati(11));
thF = rad2deg(dati_elaborati(12));

% Orbital parameters of final orbit
[rF, vF] = kep2car(aF,eF,iF,OMF,omF,thF,mu);
kepF = [aF,eF,iF,OMF,omF,thF];

% plot final orbit
[X,Y,Z] = plotOrbit(kepF,mu,360,precision);
orbitF = plot3(X,Y,Z, '--', 'LineWidth', 1);
target = plot3(rF(1), rF(2), rF(3), 'xr', 'LineWidth', 4);

% vectors to store point of the trajectory
X_traj = [];
Y_traj = [];
Z_traj = [];

fprintf('IMPLEMENTAZIONE STRATEGIA ALTERNATIVA 1:\n\n')
fprintf(['Parametri orbitali punto di partenza:\n' ...
        '\ta\t\t\te\t\ti\t\tOM\t\tom\t\tth\n ' ...
        '%4.2f\t %.4f\t %2.2f\t %2.2f\t %2.2f\t %2.2f\n\n'], kepI(1), kepI(2), kepI(3),kepI(4),kepI(5), kepI(6))
fprintf(['Parametri orbitali punto di arrivo:\n' ...
        '\ta\t\t\te\t\ti\t\tOM\t\tom\t\tth\n ' ...
        '%4.2f\t %.4f\t %2.2f\t %2.2f\t %2.2f\t %2.2f\n'], kepF(1), kepF(2), kepF(3),kepF(4),kepF(5), kepF(6))

% --- 1° Manoeuvre: change shape --- 

% to handle strange situations of transfer orbit
omT = omI;
if(strcmp(option, 'per'))
    th_man = 0;
    th_manT = 0;
    rt1 = aI*(1 - eI);
    rt2= aF*(1 + eF); %major
    if(rt2 < rt1)
        omT = mod(omF + 180,360);
        th_manT = 180;
    end
elseif (strcmp(option, 'apo'))
    th_man = 180;
    th_manT = 180;
    rt1 = aI*(1 + eI);
    rt2= aF*(1 - eF); %minor
    if(rt2 > rt1)
        omT = mod(omF + 180,360);
        th_manT = 0;
    end
else
    error('Error. Not specified manoeuvre point for change orbital shape')
end

%update trajectory from thI to first manoeuvre point th_man
if(th_man < thI)
    dTh = 360 + th_man - thI;
else
    dTh = th_man - thI;
end
[X,Y,Z] = plotOrbit(kepI,mu,dTh,precision);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

%plot from th0 to th1
orbitI = plot3(X,Y,Z, '-b', 'LineWidth',2);

% plot first manoeuvre point
[r1, v1] = kep2car(aI, eI, iI, OMI, omI, th_man, mu);
pt_changeShape1 = plot3(r1(1), r1(2), r1(3), 'ok', 'MarkerSize', 8, 'MarkerFaceColor', 'black');

%change shape -- pericenter-apocenter
[dv1, th1, dt1] = changeOrbitShape(aI, eI, omI, aF, eF, omI, thI, option);

% transfer orbit of change orbit shape manoeuvre
aT = (rt1 + rt2)/2;
eT = abs((rt2-rt1))/(rt1+rt2);

kepT = [aT, eT, iI, OMI, omI, th_manT]; %after first impulse

%plot first transfer orbit 
[X,Y,Z] = plotOrbit(kepT, mu, 180 ,precision);
orbit1 = plot3(X,Y,Z, 'k', 'LineWidth', 2);
% update trajectory from thI to th_manT
[X,Y,Z] = plotOrbit(kepT, mu, 180 ,precision);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

kep1 = [aF, eF, iI, OMI, omI, th1];  % first transfer orbit --> final form

% Plot second manoeuvre point
[r2, v2] = kep2car(aT, eT, iI, OMI, omI, 180, mu);
pt_changeShape2 = plot3(r2(1), r2(2), r2(3), 'oc', 'MarkerFaceColor','c');

% output change orbital shape
fprintf('\n---- CHANGE ORBITAL SHAPE ----\n')
fprintf(['\n.... First manoeuvre in ',option,' : ....\n'])
t_tot = dt1(1);
stampInfoManovra(kepT, th_man, dt1(1), dv1(1), t_tot);
fprintf('\n.... Second manoeuvre ....\n')
t_tot = dt1(2) + t_tot;
stampInfoManovra(kep1, th1, dt1(2), dv1(2), t_tot);

% --- 2° Manoeuvre: change Orbital Plane --- 
[dv2, om1, th2, dt2] = changeOrbitalPlane(aF, eF, iI, OMI, omI, iF, OMF, th1);

%update trajectory from th1 to second manoeuvre point th2
[X,Y,Z] = plotOrbit(kep1,mu,abs(th1-th2),precision);
orbit2 = plot3(X,Y,Z,'-c', 'LineWidth', 2);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

kep2 = [aF, eF, iF, OMF, om1, th2]; % Orbital parametre vector of second transfer orbit

% Plot second transfer orbit 
%[X,Y,Z] = plotOrbit(kep2,mu,360,precision);
%orbit2 = plot3(X,Y,Z,'--r', 'LineWidth', 1);

% Plot third manoeuvre point
[r2, v2] = kep2car(aF, eF, iF, OMF, om1, th2, mu);
pt_changePlane = plot3(r2(1), r2(2), r2(3), 'or', 'LineWidth', 3);

% output change orbital plane
fprintf('\n---- CHANGE ORBITAL PLANE ----\n')
t_tot = t_tot + dt2;
stampInfoManovra(kep2, th2, dt2, dv2, t_tot)

% --- 3° Manoeuvre: change periapsisArg --- 
[dv3, omPR, vec_th, dt3] = changePeriapsisArg(aF,eF,om1, (omF-om1), th2); 

% Check if new om is correct (equal to omF)
if(omPR ~= omF) 
    error('Errore nei calcoli')
end

th31 = vec_th(1); % th before the maneuver
th3 = vec_th(2); % th after the change of velocity

% plot third manoeuvre point
[r31, v31] = kep2car(aF, eF, iF, OMF, om1, th31, mu);
pt_changePeriapsis = plot3(r31(1), r31(2), r31(3), 'o', 'LineWidth', 2, 'Color', "#EDB120");

%update trajectory from th2 to point of maneuvre th31
if(th31 < th2)
    dTh = 360 + th31 - th2;
else
    dTh = th31-th2;
end

[X,Y,Z] = plotOrbit(kep2,mu,dTh,precision);
orbit3 = plot3(X,Y,Z,'-r', 'LineWidth', 2);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

% New orbital parameters after change of velocity
kep3 = [aF, eF, iF, OMF, omF, th3]; 

% plot third transfer orbit 
%[X,Y,Z] = plotOrbit(kep3,mu,360,precision);
%orbit3 = plot3(X,Y,Z,'--g', 'LineWidth', 1);

% output change periapsis
fprintf('\n---- CHANGE PERIAPSIS ----\n')
t_tot = t_tot + dt3;
stampInfoManovra(kep3,th31, dt3, dv3, t_tot)

% --- Go to target --- 

% time from th3 to thF on the final orbit
dt4 = timeOfFlight(aF, eF, th3, thF, mu);

% update trajectory
if(thF < th3)
    dTh = 360 + thF - th3;
else
    dTh = thF - th3; %serve abs nel caso in cui th3 sia uguale a 0
end
[X,Y,Z] = plotOrbit(kep3, mu, dTh,precision);
orbitF = plot3(X,Y,Z, 'Color', "#EDB120",'LineWidth',2);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

fprintf('\n Time to go to target: %4.2f s\n', dt4)

dt_tot = sum(dt1) + dt2 + dt3 + dt4;
dv_tot = sum(abs(dv1)) + abs(dv2) + abs(dv3);

% initialize satellite
h = plot3(nan,nan,nan,"om", 'MarkerSize',10, 'MarkerFaceColor', 'm');

%legend
legend([orbitI start orbit1 pt_changeShape1 orbit2 pt_changeShape2 orbit3 pt_changePlane  pt_changePeriapsis  orbitF target h], 'Initial Orbit', 'Initial point', 'Change shape orbit 1', 'Change shape manoeuvre point 1', 'Change shape orbit 2',  'Change shape manoeuvre point 2', 'Change plane orbit', 'Change plane manoeuvre point', 'Change periapsis manoeuvre point', 'Final Orbit', 'Target Point', 'Satellite') %sistemare plot orbita

%summary output
fprintf('\n---- Summary alternative base strategy  ----\n')
fprintf('Variazione di velocità totale: %2.3f km/s\n', dv_tot)
fprintf('Tempo totale di manovra:\n')
fprintf('Secondi: %5.2f s\n', dt_tot)
fprintf('Minuti: %5.2f m\n', dt_tot/60)
fprintf('Ore: %5.2f h\n', dt_tot/60/60)
fprintf('Giorni: %5.2f g\n', dt_tot/60/60/24)

%% Animation 
%gifFile = 'AlternativeStrategy.gif';
%exportgraphics(fig1, gifFile);

% animation
for i = 1:step_animation:length(X_traj)
    set(h,'XData',X_traj(i),'YData',Y_traj(i),'ZData',Z_traj(i));
    drawnow
    %exportgraphics(fig1, gifFile, Append=true);
end
set(h,'XData',X_traj(end),'YData',Y_traj(end),'ZData',Z_traj(end));
drawnow
