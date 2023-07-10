%% strategia 2 impulsi
%%%%%%%% INIZIALIZZAZIONE DEL PROBLEMA %%%%%%%%%%%
clear all;
close all;
clc;

global mu;

mu = 398600;
option = 'disc'; %asc o disc to choise the line node in which intersecate the orbit
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
fig1 = figure;
Terra3d;
hold on
grid on
title('Strategia 2 impulsi Secante Ottimizzando punto di partenza')

% plot intial orbit
[X,Y,Z] = plotOrbit(kepI,mu,360,precision);
orbitI = plot3(X,Y,Z, 'LineWidth', 2);
start = plot3(rI(1), rI(2), rI(3), 'xb', 'LineWidth', 4);

% plot final orbit
[X,Y,Z] = plotOrbit(kepF,mu,360,precision);
orbitF = plot3(X,Y,Z, 'LineWidth', 2);
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
point_man = plot3(r_point(1),r_point(2),r_point(3),'xm');

% check if there is some error
if(abs(r_point-r2) > toll)
    error('Errore di precisione')
end

%% Start to calculate the best orbit
th_initial = 1:3:360; % vector with angle that we want to start
th_final = 1:3:360;

% initialize the risult matrix (usefull for search minimum)
risult = [th_initial', zeros(length(th_initial),1), zeros(length(th_initial),1), zeros(length(th_initial),1)]; %[th_initial, th_final, delta_v, delta_t]

% initialize waitbar
wb = waitbar(0,'Calcolo Orbite');

for k = 1:length(th_initial)
    % find r1
    r1 = kep2car(aI,eI,iI,OMI,omI,th_initial(k),mu);
    temp = [th_final', zeros(length(th_final),1), zeros(length(th_final),1)]; %[th_final, delta_v, delta_t]
    for j = 1:length(th_final)
        % find r1
        r2 = kep2car(aF,eF,iF,OMF,omF,th_final(j),mu);
    
        % find the best transfert orbit from r1 to r2
        [kepT, dv, th] = findOrbit(r1,r2,kepI,kepF);
    
        % for plotting same example 
        %{
        if (mod(k,100) == 0)
            [X,Y,Z] = plotOrbit(kepT,mu,360,0.1);
            plot3(X,Y,Z);
            plot3(r1(1),r1(2),r1(3), 'or','LineWidth',4);
            text(r1(1),r1(2),r1(3), num2str(th_initial(k),3))
        end
        %}
    
        % calculate the different time from thI to thF
        t1 = timeOfFlight(aI,eI,thI,th(1),mu);
        t2 = timeOfFlight(kepT(1),kepT(2),th(2),th(3),mu);
        t3 = timeOfFlight(aF,eF,th(4),thF,mu);
    
        % save result
        temp(k,2) = dv;
        temp(k,3) = t1 + t2 + t3;
    end

    % clean the risult matrix from null risult
    for z = length(temp):-1:1
        if(temp(z,2) == 0)
            temp(z,:) = []; %cancello la riga
        end
    end

    % incollo le righe
    % clean the risult matrix from null risult
    if(isempty(temp))
        temp = [0,0,0];
    else
        [vmin, kmin] = min(temp(:,2));
        risult(k,2:end) = temp(kmin,:);
    end

    waitbar(k/length(th_initial), wb, sprintf('Calcolo Orbite: %2.0f %%', k/length(th_initial)*100));
end

delete(wb)

% clean the risult matrix from null risult
for k = length(risult):-1:1
    if(risult(k,3) == 0)
        risult(k,:) = []; %cancello la riga
    end
end

% plot the best transfert orbit
figure(fig1)
[vmin, kvmin] = min(risult(:,3));
[r1,v1] = kep2car(aI,eI,iI,OMI,omI,risult(kvmin,1),mu);
[kepT, dv, th] = findOrbit(r1,r2,kepI,kepF);
[X,Y,Z] = plotOrbit(kepT,mu,360,0.1);
orbitaT = plot3(X,Y,Z);

fprintf('\nRisultato migliore è di:')
fprintf('dv: %2.2f km/s\n', vmin)
fprintf('th1: %2.2f gradi\n', th(1))
fprintf('th2: %2.2f gradi\n', th(4))