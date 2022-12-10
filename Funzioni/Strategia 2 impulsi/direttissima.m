%% strategia 2 impulsi Direttissima
clear;
clc;
close all;

global mu;

mu = 398600;
ALLERT = 6378 + 300; %Km

% Ancora parecchio da sistemare
% Pulire codice da cose non utili
% Controllare le velocità
% Bug: se abbasso il meshdensity non converge la soluzione
% Aggiungere animazione

option = 'asc'; %asc o disc
precision = 0.1;
toll = 1e-4;

%%%%%%%% GRUPPO B7 %%%%%%%%
dati_elaborati = [5088.9118 -3196.5659 -8222.7989 1.9090 5.6220 -1.0700 14020.0000 0.3576 1.3220 0.9764 1.8130 0.4336];

%altri gruppi per prova
%dati_elaborati = [-7394.8822 849.6501 4181.2069 -3.0970 -5.2110 -3.5690 15660.0000 0.2659 0.8321 1.2440 2.9320 3.0970];

%dati_elaborati = [-6096.8804 -1361.2133 4888.2313 -0.2075 -7.2520 -1.4010 12860.0000 0.2842 1.4520 0.3340 2.7340 2.9360];

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
fig1 = figure;
Terra3d;
hold on
grid on
title('Strategia base')

%orbita iniziale
[X,Y,Z] = plotOrbit(kepI,mu,360,precision);
orbitI = plot3(X,Y,Z, 'LineWidth', 2);
start = plot3(rI(1), rI(2), rI(3), 'xb', 'LineWidth', 4);

%orbita finale
[X,Y,Z] = plotOrbit(kepF,mu,360,precision);
orbitF = plot3(X,Y,Z, 'LineWidth', 2);
target = plot3(rF(1), rF(2), rF(3), 'xr', 'LineWidth', 4);



[kepT, dv, th] = findOrbit(rI,rF,kepI,kepF,option);
[X,Y,Z] = plotOrbit(kepT,mu,360,0.1);
plot3(X,Y,Z);

[~,v1] = kep2car(kepT(1),kepT(2),kepT(3),kepT(4),kepT(5),th(2),mu);
[~,v2] = kep2car(kepT(1),kepT(2),kepT(3),kepT(4),kepT(5),th(3),mu);
t = timeOfFlight(kepT(1),kepT(2),th(2),th(3),mu);

fprintf('IMPLEMENTAZIONE STRATEGIA 2 IMPULSI DIRETTA:\n')
fprintf('Costo della manovra: %2.3f km/s\n', dv)
fprintf('Tempo di manovra: %4.1f s\n', t)

fprintf('\n--- Prima manovra ---\n')
stampInfoManovra(kepT,th(1),0,norm(v1-vI))

fprintf('\n--- Seconda manovra ---\n')
stampInfoManovra(kepF,th(3),t,norm(vF-v2))

