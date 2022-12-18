%% STRATEGIA BIELLITTICA

clear;
clc;
close all;

global mu;

mu = 398600;


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
title('Strategia biellittica')

%orbita iniziale
[X,Y,Z] = plotOrbit(kepI,mu,360,precision);
orbitI = plot3(X,Y,Z, 'LineWidth', 2);
start = plot3(rI(1), rI(2), rI(3), 'xb', 'LineWidth', 4);

%orbita finale
[X,Y,Z] = plotOrbit(kepF,mu,360,precision);
orbitF = plot3(X,Y,Z, 'LineWidth', 2);
target = plot3(rF(1), rF(2), rF(3), 'xr', 'LineWidth', 4);

fprintf('IMPLEMENTAZIONE STRATEGIA BIELLITTICA:\n\n')

fprintf(['Parametri orbitali punto di partenza:\n' ...
        '\ta\t\t\te\t\ti\t\tOM\t\tom\t\tth\n ' ...
        '%4.2f\t %.4f\t %2.2f\t %2.2f\t %2.2f\t %2.2f\n\n'], kepI(1), kepI(2), kepI(3),kepI(4),kepI(5), kepI(6))
fprintf(['Parametri orbitali punto di arrivo:\n' ...
        '\ta\t\t\te\t\ti\t\tOM\t\tom\t\tth\n ' ...
        '%4.2f\t %.4f\t %2.2f\t %2.2f\t %2.2f\t %2.2f\n'], kepF(1), kepF(2), kepF(3),kepF(4),kepF(5), kepF(6))

% calcolo l'angolo di manovra del cambio piano rispetto alla prima orbita
% --> in questo modo darò l'impulso 180° prima per arrivare al p.to di
% cambio piano nell'apocentro

[~, ~, th_cp, ~] = changeOrbitalPlane(aI, eI, iI, OMI, omI, iF, OMF, thI);
% disegno th_cp per capire se è il nodo ascendente/discendente:
[r_cp, v_cp] = kep2car(aI, eI, iI, OMI, omI, th_cp, mu);
cp = plot3(r_cp(1), r_cp(2), r_cp(3), 'xb', 'LineWidth', 4); % DOVREBBE essere il nodo ascendente
% punto di inizio della biellittica:
th_in_be = th_cp - 180;
% disegno il punto a cui inizia la biellittica:
[r_in_be, v_in_be] = kep2car(aI, eI, iI, OMI, omI, th_in_be, mu);
in_be = plot3(r_in_be(1), r_in_be(2), r_in_be(3), 'xb', 'LineWidth', 4);
% siccome th_in_be è l'angolo opposto al nodo ascendente, è il nodo
% discendente --> quindi om_T = 180 (ricorda che om è l'angolo dal nodo
% ascendente al pericentro)
% trovato om_T, faccio il cambio di periasse per raddrizzare l'orbita di
% trasferimento

% manovra 1: primo cambio di periasse
om_T = 180;
[dv1, om_T, vec_th1, dt1] = changePeriapsisArg(aI, eI, omI, (om_T-omI), thI);

th11 = vec_th1(1);
th1 = vec_th1(2);

% disegno punto di manovra 2
[r11, v11] = kep2car(aI, eI, iI, OMI, om_T, th11, mu);
pt_changePeriapsis1 = plot3(r11(1), r11(2), r11(3), '^g', 'LineWidth', 2);

%aggiorno traettoria da th1 a punto di manovra th21
if(th11 < thI)
    dTh = 360 + th11 - thI;
else
    dTh = th11-thI;
end

[X,Y,Z] = plotOrbit(kepI,mu,dTh,precision);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

kep1= [aI, eI, iI, OMI, om_T, th1]; %Nuova orbita dopo cambio di periasse


% Disegno nuova orbita dopo cambio di periasse 
[X,Y,Z] = plotOrbit(kep1,mu,360,precision);
orbit1 = plot3(X,Y,Z,'--g', 'LineWidth', 1);

% output primo cambio di periasse
fprintf('\n---- PRIMO CAMBIO DI PERIASSE ----\n')
stampInfoManovra(kep1,th11, dt1, dv1)

% dobbiamo ora raggiungere il pericentro dell'orbita --> da qui inizia la
% biellittica
th2 = 0; % ! va controllato che th2 corrisponda a th_in_be !
[rth2, vth2] = kep2car(aI, eI, iI, OMI, om_T, th2, mu);

if rth2 ~= r_in_be
    disp('ERRORE: BIELLITTICA NON ALLINEATA AD ASSE DEI NODI')
end

errore = abs(rth2-r_in_be);

% qui dobbiamo "scegliere" rb --> per ora lo scelgo a caso, poi va ottimizzato con un for 
% --> DA QUI VA FINITA MA IN TEORIA è OK, IL
% PROBLEMA è CHE DOBBIAMO ALLINEARE LA BIELLITTICA ALL'ASSE DEI NODI (???)
rb = 3*aF;
rp_t1 = aI*(1-eI); % perché rp_t1 = rp_I
ra_t1 = rb;
a_t1 = (rp_t1+ra_t1)/2;
e_t1 = (ra_t1-rp_t1)/(ra_t1+rp_t1);
option = 'per';
% manovra 2: prima ellisse della biellittica
[dv2, th3, dt2] = changeOrbitShape(aI, eI, om_T, a_t1, e_t1, om_T, th2, option);

% aggiorno traettoria da th1 a punto di manovra th2
[X,Y,Z] = plotOrbit(kep1,mu,(abs(th1-th2)),precision);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

% disegno prima parte della biellittica:
kep_t1 = [a_t1, e_t1, iI, OMI, om_T, th3];
[X,Y,Z] = plotOrbit(kep_t1,mu,360,precision);
orbit2 = plot3(X,Y,Z,'--m', 'LineWidth', 1);

% output prima parte di biellittica
fprintf('\n---- PRIMO IMPULSO BIELLITTICA ----\n')
stampInfoManovra(kep_t1,th2, sum(abs(dt2)), sum(abs(dv2)))

% aggiorno traiettoria da th2 a th3
% aggiorno traettoria da th2 a punto di cambio piano th3
[X,Y,Z] = plotOrbit(kep_t1,mu,(abs(th3-th2)),precision);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

% arrivati a th3 --> cambio piano e cambio forma (OVVIAMENTE SICCOME LA
% BIELLITTICA è STORTA, th3 è l'apocentro ma non è l'angolo di manovra del
% cambio piano --> motivo per cui dt3 =/= 0)
% prima solo cambio piano
[dv3, om3, th4, dt3] = changeOrbitalPlane(a_t1, e_t1, iI, OMI, om_T, iF, OMF, th3); 

% aggiorno traettoria
[X,Y,Z] = plotOrbit(kepI,mu,abs(th4-th3),precision);
X_traj = X;
Y_traj = Y;
Z_traj = Z;

kep3 = [a_t1, e_t1, iF, OMF, om3, th4]; %1 orbita di trasferimento

% Disegno prima orbita di trasferimento 
[X,Y,Z] = plotOrbit(kep3,mu,360,precision);
orbit1 = plot3(X,Y,Z,'--r', 'LineWidth', 1);

% disegno il punto punto di manovra 1
[r3, v3] = kep2car(a_t1, e_t1, iF, OMF, om3, th4, mu);
pt_changePlane = plot3(r3(1), r3(2), r3(3), 'ok', 'LineWidth', 2);

% output cambio di piano
fprintf('\n---- CAMBIO DI PIANO ----\n')
stampInfoManovra(kep3, th3, dt3,dv3)

% ... --> bisogna mettere a posto la biellittica prima di andare avanti


