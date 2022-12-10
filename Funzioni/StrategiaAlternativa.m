%% Prima strategia alternativa: cambio forma e poi cambio piano
clear;
clc;
close all;

global mu;

mu = 398600;

option = 'per'; %apo or per
precision = 0.1;
step_animation = 10;

%%%%%%%% GRUPPO B7 %%%%%%%%
dati_elaborati = [5088.9118 -3196.5659 -8222.7989 1.9090 5.6220 -1.0700 14020.0000 0.3576 1.3220 0.9764 1.8130 0.4336];

%altri gruppi per prova
%dati_elaborati = [-7394.8822 849.6501 4181.2069 -3.0970 -5.2110 -3.5690 15660.0000 0.2659 0.8321 1.2440 2.9320 3.0970];

%dati orbita iniziale
rI = dati_elaborati(1:3)';
vI = dati_elaborati(4:6)';

%calcolo altri parametri orbita iniziale
[aI, eI, iI, OMI, omI, thI] = car2kep(rI, vI, mu);
kepI = [aI,eI,iI,OMI,omI,thI];

%dati orbita finale
aF = dati_elaborati(7);
eF = dati_elaborati(8);
iF = rad2deg(dati_elaborati(9));
OMF = rad2deg(dati_elaborati(10));
omF = rad2deg(dati_elaborati(11));
thF = rad2deg(dati_elaborati(12));

%calcolo altri parametri orbita finale
[rF, vF] = kep2car(aF,eF,iF,OMF,omF,thF,mu);
kepF = [aF,eF,iF,OMF,omF,thF];

%vettori per tacciare orbita
X_traj = [];
Y_traj = [];
Z_traj = [];

%definizione fig
figure
Terra3d
hold on
grid on
title('Strategia Alternativa')

%orbita iniziale
[X,Y,Z] = plotOrbit(kepI,mu,360,precision);
orbitI = plot3(X,Y,Z, 'LineWidth', 2);
start = plot3(rI(1), rI(2), rI(3), 'xb', 'LineWidth', 4);

%orbita finale
[X,Y,Z] = plotOrbit(kepF,mu,360,precision);
orbitF = plot3(X,Y,Z, 'LineWidth', 2);
target = plot3(rF(1), rF(2), rF(3), 'xr', 'LineWidth', 4);

fprintf('IMPLEMENTAZIONE STRATEGIA ALTERNATIVA 1:\n\n')
fprintf(['Parametri orbitali punto di partenza:\n' ...
        '\ta\t\t\te\t\ti\t\tOM\t\tom\t\tth\n ' ...
        '%4.2f\t %.4f\t %2.2f\t %2.2f\t %2.2f\t %2.2f\n\n'], kepI(1), kepI(2), kepI(3),kepI(4),kepI(5), kepI(6))
fprintf(['Parametri orbitali punto di arrivo:\n' ...
        '\ta\t\t\te\t\ti\t\tOM\t\tom\t\tth\n ' ...
        '%4.2f\t %.4f\t %2.2f\t %2.2f\t %2.2f\t %2.2f\n'], kepF(1), kepF(2), kepF(3),kepF(4),kepF(5), kepF(6))

% cambio di forma: Manovra 1
% dal pericentro dell'orbita 1 all'apocentro dell'orbita 2:

% gestione casi strani orbtia di trasferimento
omT = omI;
if(strcmp(option, 'per'))
    th_man = 0;
    th_manT = 0;
    rt1 = aI*(1 - eI);
    rt2= aF*(1 + eF); %maggiore
    if(rt2 < rt1)
        omT = mod(omF + 180,360);
        th_manT = 180;
    end
elseif (strcmp(option, 'apo'))
    th_man = 180;
    th_manT = 180;
    rt1 = aI*(1 + eI);
    rt2= aF*(1 - eF); %minore
    if(rt2 > rt1)
        omT = mod(omF + 180,360);
        th_manT = 0;
    end
else
    error('Error. Punto di manovra per cambio forma non specificato')
end

%aggiorno traettoria da thI a punto di manovra th_man
if(th_man < thI)
    dTh = 360 + th_man - thI;
else
    dTh = th_man - thI;
end
[X,Y,Z] = plotOrbit(kepI,mu,dTh,precision);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

%stampo punto di manovra cambio forma
[r1, v1] = kep2car(aI, eI, iI, OMI, omI, th_man, mu);
pt_changeShape = plot3(r1(1), r1(2), r1(3), 'oc', 'LineWidth', 1);

%trasferimento pericentro-apocentro
[dv1, th1, dt1] = changeOrbitShape(aI, eI, omI, aF, eF, omI, thI, option);

% orbita di trasferimento cambio forma
aT = (rt1 + rt2)/2;
eT = abs((rt2-rt1))/(rt1+rt2);

kepT = [aT, eT, iI, OMI, omT, th_manT];%dopo primo impulso
%plot orbita di trasferimento e aggiorno le coordinate
[X,Y,Z] = plotOrbit(kepT, mu, 360,precision);
orbit1 = plot3(X,Y,Z, '--k', 'LineWidth', 2);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

kep1 = [aF, eF, iI, OMI, omI, th1];  % 1 orbita di trasferimento --> forma finale

% output cambio di forma
fprintf('\n---- CAMBIO DI FORMA ----\n')
fprintf(['\n.... Prima manovra in ',option,' : ....\n'])
stampInfoManovra(kepT, th_man, dt1(1), dv1(1));
fprintf('\n.... Seconda manovra ....\n')
stampInfoManovra(kep1, th1, dt1(2), dv1(2));

% cambio di piano : Manovra 2

[dv2, om1, th2, dt2] = changeOrbitalPlane(aF, eF, iI, OMI, omI, iF, OMF, th1); %forse problema dt

%aggiorno traettoria da th1 a punto di manovra th2
[X,Y,Z] = plotOrbit(kep1,mu,abs(th1-th2),precision);
X_traj = X;
Y_traj = Y;
Z_traj = Z;

kep2 = [aF, eF, iF, OMF, om1, th2]; % 2 orbita di trasferimento

% Disegno la seconda orbita di trasferimento 
[X,Y,Z] = plotOrbit(kep2,mu,360,precision);
orbit2 = plot3(X,Y,Z,'--r', 'LineWidth', 1);

% disegno il punto punto di manovra 2
[r2, v2] = kep2car(aF, eF, iF, OMF, om1, th2, mu);
pt_changePlane = plot3(r2(1), r2(2), r2(3), 'ok', 'LineWidth', 2);

% output cambio di piano
fprintf('\n---- CAMBIO DI PIANO ----\n')
stampInfoManovra(kep2, th2, dt2, dv2)

% cambio di anomalia del pericentro: Manovra 3
[dv3, omPR, vec_th, dt3] = changePeriapsisArg(aF,eF,om1, (omF-om1), th2); 

if(omPR ~= omF) 
    error('Errore nei calcoli')
end

th31 = vec_th(1); % anomalia di manovra rispetto all'orbita iniziale
th3 = vec_th(2); % anomalia di manovra sull'orbita finale

% disegno punto di manovra 3
[r31, v31] = kep2car(aF, eF, iF, OMF, om1, th31, mu);
pt_changePeriapsis = plot3(r31(1), r31(2), r31(3), '^g', 'LineWidth', 2);

%aggiorno traettoria da th2 a punto di manovra th21
if(th31 < th2)
    dTh = 360 + th31 - th2;
else
    dTh = th31-th2;
end

[X,Y,Z] = plotOrbit(kep2,mu,dTh,precision);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

kep3 = [aF, eF, iF, OMF, omF, th3]; %Nuova orbita dopo cambio di velocità

% Disegno seconda orbita di trasferimento 
[X,Y,Z] = plotOrbit(kep3,mu,360,precision);
orbit3 = plot3(X,Y,Z,'--g', 'LineWidth', 1);

% output cambio di periasse
fprintf('\n---- CAMBIO DI PERIASSE ----\n')
stampInfoManovra(kep3,th31, dt3, dv3)

% ultimo pezzo: arrivo fino al target
dt4 = timeOfFlight(aF, eF, th3, thF, mu);

if(thF < th3)
    dTh = 360 + thF - th3;
else
    dTh = thF - th3; %serve abs nel caso in cui th3 sia uguale a 0
end
[X,Y,Z] = plotOrbit(kep3, mu, dTh,precision);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

fprintf('\nAttesa fino al punto finale: %4.2f s\n', dt4)

dt_tot = sum(dt1) + dt2 + dt3 + dt4;
dv_tot = sum(abs(dv1)) + abs(dv2) + abs(dv3);

%satellite
h = plot3(nan,nan,nan,"om", 'LineWidth',4);

%legend
legend([orbitI orbitF start orbit1 pt_changeShape orbit2 pt_changePlane target pt_changePeriapsis orbit3 h], 'Orbita Iniziale', 'Orbita Finale','Punto Iniziale', 'Orbita di cambio forma', 'Punto di manovra cambio forma', 'Orbita di cambio piano', 'Punto di manovra cambio piano', 'Punto finale', 'Punto di manovra cambio periasse', 'Manovra di trasferimento cambio forma', 'Satellite') %sistemare plot orbita

%output riassuntivo
fprintf('\n---- Riassunto strategia alternativa ----\n')
fprintf('Variazione di velocità totale: %2.3f km/s\n', dv_tot)
fprintf('Tempo totale di manovra:\n')
fprintf('Secondi: %5.2f s\n', dt_tot)
fprintf('Minuti: %5.2f h\n', dt_tot/60)
fprintf('Ore: %5.2f d\n', dt_tot/60/60)
fprintf('Giorni: %5.2f g\n', dt_tot/60/60/24)

%disegno satellite che segue la traettoria
for i = 1:step_animation:length(X_traj)
    set(h,'XData',X_traj(i),'YData',Y_traj(i),'ZData',Z_traj(i));
    drawnow
end
set(h,'XData',X_traj(end),'YData',Y_traj(end),'ZData',Z_traj(end));
drawnow

%% plot meno incasinato:
%definizione fig
figure
Terra3d
hold on
grid on
title('Strategia Alternativa')

%orbita iniziale
[X,Y,Z] = plotOrbit(kepI,mu,360,precision);
orbitI = plot3(X,Y,Z, 'LineWidth', 2);
start = plot3(rI(1), rI(2), rI(3), 'xb', 'LineWidth', 4);

%orbita finale
[X,Y,Z] = plotOrbit(kepF,mu,360,precision);
orbitF = plot3(X,Y,Z, 'LineWidth', 2);
target = plot3(rF(1), rF(2), rF(3), 'xr', 'LineWidth', 4);

% traettoria da thI a punto di manovra th_man
[X,Y,Z] = plotOrbit(kepI,mu,abs(thI-th_man),precision);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

% punto di manovra cambio forma
[r1, v1] = kep2car(aI, eI, iI, OMI, omI, th_man, mu);
pt_changeShape = plot3(r1(1), r1(2), r1(3), 'oc', 'LineWidth', 1);

%plot orbita di trasferimento e aggiorno le coordinate
[X,Y,Z] = plotOrbit(kepT, mu, abs(th_man - th2),precision);
orbit1 = plot3(X,Y,Z, '--k', 'LineWidth', 2);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

% cambio di piano : Manovra 2
% traettoria da th1 a punto di manovra th2
[X,Y,Z] = plotOrbit(kep1,mu,abs(th1-th2),precision);
X_traj = X;
Y_traj = Y;
Z_traj = Z;

% punto punto di manovra 2
[r2, v2] = kep2car(aF, eF, iF, OMF, om1, th2, mu);
pt_changePlane = plot3(r2(1), r2(2), r2(3), 'ok', 'LineWidth', 2);

% Disegno la seconda orbita di trasferimento 
[X,Y,Z] = plotOrbit(kep2,mu,abs(th2-th3),precision);
orbit2 = plot3(X,Y,Z,'--r', 'LineWidth', 1);

% punto di manovra 3
[r31, v31] = kep2car(aF, eF, iF, OMF, om1, th31, mu);
pt_changePeriapsis = plot3(r31(1), r31(2), r31(3), '^g', 'LineWidth', 2);

%aggiorno traettoria da th2 a punto di manovra th31
[X,Y,Z] = plotOrbit(kep2,mu,abs(th2-th3),precision);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

% ultimo pezzo: arrivo fino al target
% aggiorno traiettoria
[X,Y,Z] = plotOrbit(kep3, mu, (thF-th3 + 360),precision);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

%satellite
h = plot3(nan,nan,nan,"om", 'LineWidth',4);

%disegno satellite che segue la traettoria
for i = 1:step_animation:length(X_traj)
    set(h,'XData',X_traj(i),'YData',Y_traj(i),'ZData',Z_traj(i));
    drawnow
end
set(h,'XData',X_traj(end),'YData',Y_traj(end),'ZData',Z_traj(end));
drawnow
