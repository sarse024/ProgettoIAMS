%% strategia 2 impulsi
%%%%%%%% INIZIALIZZAZIONE DEL PROBLEMA %%%%%%%%%%%
clear all;
close all;
clc;

global mu;

mu = 398600;
option = 'asc'; %asc o disc to choise the line node in which intersecate the orbit
precision = 0.1;
toll = 1e-4;
step_animation = 50;

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

% vector for store point of the trajectory
X_traj = [];
Y_traj = [];
Z_traj = [];

% Initialize main figure
fig1 = figure('WindowState','maximized');
Terra3d;
hold on
grid on
title('Optimized Secant Manoeuvre', fontsize=20)

% plot intial orbit
[X,Y,Z] = plotOrbit(kepI,mu,360,precision);
plot3(X,Y,Z, '--', 'LineWidth', 1);
start = plot3(rI(1), rI(2), rI(3), 'xb', 'LineWidth', 4);

% plot final orbit
[X,Y,Z] = plotOrbit(kepF,mu,360,precision);
plot3(X,Y,Z, '--', "MarkerFaceColor", "#EDB120", 'LineWidth', 1);
target = plot3(rF(1), rF(2), rF(3), 'xr', 'LineWidth', 4);

fprintf('IMPLEMENTAZIONE STRATEGIA 2 IMPULSI OTTIMIZZATA punto di partenza:\n\n')

%% find the th for the asc or disc node

% define plane 
hI = cross(rI, vI);
hF = cross(rF, vF);

if(strcmp(option,'asc'))
    vector = cross(hI,hF); % whit the minus will obtain the other intersection
else
    vector = -cross(hI,hF); % for 'disc'
end

% rotation matrix 
R_F = rotationMatrix(iF,OMF,omF);

vec = R_F*vector; %position vector in the perifocae system
th_intersection = atan2(vec(2),vec(1)); %find th

% adjust angle for atan
if(th_intersection < 0)
    th_intersection = 2*pi+th_intersection;
end

r_inters = aF*(1-eF^2)/(1+eF*cos(th_intersection)); %find the radius of the intesection

% define the r2 on vector line
t = r_inters/norm(vector);
r2 = vector*t;

% calculate the position vector of the point of maneuvre
[r_point,v_point]  = kep2car(aF,eF,iF,OMF,omF, rad2deg(th_intersection),mu);

% plot the point of maneuvre
second_point_man = plot3(r_point(1),r_point(2),r_point(3),'xm', 'LineWidth',3);

% check if there is some error
if(abs(r_point-r2) > toll)
    error('Errore di precisione')
end

%% Start to calculate the best orbit
th_initial = 0:10:359; % vector with angle that we want to start

% initialize the risult matrix (usefull for search minimum)
risult = [th_initial', zeros(length(th_initial),1), zeros(length(th_initial),1)]; %[th_initial, delta_v, delta_t]

% initialize waitbar
wb = waitbar(0,'Calcolo Orbite');

for k = 1:length(th_initial)
    % find r1
    r1 = kep2car(aI,eI,iI,OMI,omI,th_initial(k),mu);

    % find the best transfert orbit from r1 to r2
    [kepT, dv, th] = findOrbit(r1,r2,kepI,kepF);

    % for plotting same example 
    %{
    %if (mod(k,6) == 0)
        [X,Y,Z] = plotOrbit(kepT,mu,360,0.1);
        plot3(X,Y,Z);
        plot3(r1(1),r1(2),r1(3), 'or','LineWidth',4);
        text(r1(1),r1(2),r1(3), strcat(num2str(th_initial(k),3), '°'), fontsize=20)
    %end
    %}

    % calculate the different time from thI to thF
    t1 = timeOfFlight(aI,eI,thI,th(1),mu);
    t2 = timeOfFlight(kepT(1),kepT(2),th(2),th(3),mu);
    t3 = timeOfFlight(aF,eF,th(4),thF,mu);

    % save result
    risult(k,2) = dv;
    risult(k,3) = t1 + t2 + t3;
    waitbar(k/length(th_initial), wb, sprintf('Calcolo Orbite: %2.0f %%', k/length(th_initial)*100));
end

delete(wb)

% clean the risult matrix from null risult
for k = length(risult):-1:1
    if(risult(k,2) == 0)
        risult(k,:) = []; %cancello la riga
    end
end

% plot the best transfert orbit
figure(fig1)
[vmin, kvmin] = min(risult(:,2));
[r1,v1] = kep2car(aI,eI,iI,OMI,omI,risult(kvmin,1),mu);
[kepT, dv, th] = findOrbit(r1,r2,kepI,kepF);

[X,Y,Z] = plotOrbit(kepT,mu,360, precision);
plot3(X,Y,Z, '--r','LineWidth',1);

[X,Y,Z] = plotOrbit(kepT,mu, mod(th(3)-th(2),360), precision);
orbitT = plot3(X,Y,Z,'r','LineWidth', 2);

%% create graphs for comparison
figure
tiledlayout(2,1)
nexttile
plot(risult(:,1),risult(:,2),'-o');
title('Total Velocity variation with True Anomaly at Departure  ', fontsize=20);
xlabel('True Anomaly at Departure [deg]', fontsize=15), ylabel('Velocity [km/s]', fontsize=15)
nexttile
plot(risult(:,1),risult(:,3)/60/60,'-o');
xlabel('True Anomaly at Departure [deg]', fontsize=15);
ylabel('Time [h]', fontsize=15)
title('Total Time variation with True Anomaly at Departure', fontsize=20);

%% graphic test and animation
figure(fig1)

% to make the code more readable
th1 = th(1);
new_th1 = th(2);
new_th2 = th(3);
th2 = th(4);

[rprov1, vprov1] = kep2car(kepT(1),kepT(2),kepT(3),kepT(4), kepT(5), new_th1,mu);
first_point_man = plot3(rprov1(1),rprov1(2),rprov1(3),'xc', 'LineWidth',5);

[rprov2, vprov2] = kep2car(kepT(1),kepT(2),kepT(3),kepT(4), kepT(5), new_th2,mu);

 % calculate the different time from thI to thF of the best transfert orbit
t1 = timeOfFlight(aI,eI,thI,th1,mu);
t2 = timeOfFlight(kepT(1),kepT(2),new_th1,new_th2,mu);
t3 = timeOfFlight(aF,eF,th2,thF,mu);
t = t1 + t2 + t3;

%%% OUTPUT %%%

fprintf('%%%% MANOVRA CON V_tot MINIMA %%%%\n\n')
fprintf('--- Prima manovra ---\n')
stampInfoManovra(kepT, th1, t1, norm(vprov1-v1))
fprintf('Errore sul primo punto di manovra: %2.2f m\n\n', norm(rprov1-r1)*1000)

kep2 = [aF,eF,iF,OMF,omF,th(4)];
fprintf('--- Seconda manovra ---\n')
stampInfoManovra(kep2, new_th2, t2, norm(v_point-vprov2))
fprintf('Errore sul secondo punto di manovra: %2.2f m\n\n', norm(r_point-rprov2)*1000)

fprintf('Tempo di attesa fino al target della missione dopo 2 impulso: %4.2f\n\n', t3)

fprintf('--- Riassunto della manovra ---\n')
fprintf('Costo Totale della manovra: %2.4f\n', vmin)
fprintf('Tempo totale di manovra:\n')
fprintf('Secondi: %5.2f s\n', t)
fprintf('Minuti: %5.2f m\n', t/60)
fprintf('Ore: %5.2f h\n', t/60/60)
fprintf('Giorni: %5.2f d\n', t/60/60/24)

%%% ANIMATION %%%

if(th1 < thI)
    dTh = 360 + th1 - thI;
else
    dTh = th1 - thI;
end

% drawn on the initial orbit to reach the first maneuvering point
[X_traj,Y_traj,Z_traj] = plotOrbit(kepI,mu,dTh,precision);
orbitI = plot3(X_traj,Y_traj,Z_traj, 'b', 'LineWidth', 2);

if(new_th2 < new_th1)
    dTh = 360 + new_th2 - new_th1;
else
    dTh = new_th2 - new_th1; 
end

% drawn on the transfer orbit to reach the second maneuver point
[X,Y,Z] = plotOrbit(kepT,mu,dTh,precision);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

kepF = [aF,eF,iF,OMF,omF,th2];

if(thF < th2)
    dTh = 360 + thF - th2;
else
    dTh = thF - th2;
end

% stretched on the final orbit to reach the final point
[X,Y,Z] = plotOrbit(kepF,mu,dTh,precision);
orbitF = plot3(X,Y,Z, "MarkerFaceColor", "#EDB120", 'LineWidth', 2);

X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

% initialize satellite
h = plot3(nan,nan,nan,"om", 'LineWidth',4);

% legend
legend([orbitI, start, orbitF, target, orbitT, first_point_man, second_point_man, h], 'Initial Orbit', 'Starting point', 'Final Orbit', 'Target point', 'Transfer Trajectory', 'First point of manoeuvre', 'Second point of manoeuvre', 'Satellite', 'Location','best')
xlabel('X axis [Km]'), ylabel('Y axis [Km]'), zlabel('Z axis [Km]')

% to save gif animation
%pause(10)
%gifFile = '2impulsi.gif';
%exportgraphics(fig1, gifFile);

% animation
for i = 1:step_animation:length(X_traj)
    set(h,'XData',X_traj(i),'YData',Y_traj(i),'ZData',Z_traj(i));
    %exportgraphics(fig1, gifFile, Append=true);
    drawnow;
end
set(h,'XData',X_traj(end),'YData',Y_traj(end),'ZData',Z_traj(end));
drawnow



