% Biellittica Circolarizzata OTTIMIZZATA

% NB: 1) va fatto il plot fatto bene, 2) ho usato una versione modificata del
% changeOrbitShape per non avere le due componenti di dv e dt da sommare
% (era un casino da gestire col for)

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

fprintf('IMPLEMENTAZIONE STRATEGIA BIELLITTICA 2:\n\n')

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
% punto di inizio della biellittica:
th_in_be = th_cp - 180;
% disegno il punto a cui inizia la biellittica:
[r_in_be, v_in_be] = kep2car(aI, eI, iI, OMI, omI, th_in_be, mu);
in_be = plot3(r_in_be(1), r_in_be(2), r_in_be(3), 'xb', 'LineWidth', 4);

% manovra 1: CIRCOLARIZZAZIONE ORBITA INIZIALE (all'apocentro -> fino al "pericentro" dell'orbita circolare) 

raI = aI*(1+eI);
a1 = raI;
e1 = 0;
[dv1, th1, dt1] = changeOrbitShape(aI, eI, omI, a1, e1, omI, thI, 'apo');
thman1 = 180;
dth1 = thman1 - thI;
kep1 = [a1, e1, iI, OMI, omI, th1];
% output prima parte di biellittica
fprintf('\n---- CIRCOLARIZZAZIONE 1 ----\n')
stampInfoManovra(kep1,thman1, sum(abs(dt1)), sum(abs(dv1)))

% dobbiamo ora raggiungere il punto di inizio della biellittica
rbv = linspace(aF, 5*aF, 10000);
a_t1v = [];
e_t1v = [];
ra_t1v = [];
dv2v = [];
dt2v = []; 
dv3v = [];
dt3v = [];
ra3v = [];
a_t2v = [];
e_t2v = [];
dv4v = [];
dt4v = [];
dv5v = [];
dt5v = [];
dv_totv = [];
dt_totv = [];

%  OTTIMIZZAZIONE
for k = 1 : 10000
rb = rbv(k);
rp_t1 = a1; % perché rp_t1 = r1
ra_t1v(k) = rb;
a_t1v(k) = (rp_t1+ra_t1v(k))/2;
e_t1v(k) = (ra_t1v(k)-rp_t1)/(ra_t1v(k)+rp_t1);
a_t1v = [a_t1v; a_t1v(k)];
e_t1v = [e_t1v; e_t1v(k)];
ra_t1v = [ra_t1v; ra_t1v(k)];


om_t1 = omI + th_in_be;
% manovra 2: prima ellisse della biellittica
[dv2v(k), th2, dt2v(k)] = changeOrbitShapeALT(a1, e1, om_t1, a_t1v(k), e_t1v(k), om_t1, th_in_be, 'per');
dv2v = [dv2v; dv2v(k)];
dt2v = [dt2v; dt2v(k)];

% arrivati a th2 --> cambio piano e cambio forma
% prima solo cambio piano
[dv3v(k), om3, th3, dt3v(k)] = changeOrbitalPlane(a_t1v(k), e_t1v(k), iI, OMI, om_t1, iF, OMF, (th2-1e-7)); 
dv3v = [dv3v; dv3v(k)];
dt3v = [dt3v; dt3v(k)];

% SECONDA PARTE DELLA BIELLITTICA
ra3v(k) = rb;
rp3 = aF*(1+eF);
a_t2v(k) = (ra3v(k) + rp3)/2;
e_t2v(k) = (ra3v(k)-rp3)/(ra3v(k)+rp3);

[dv4v(k), th4, dt4v(k)] = changeOrbitShapeALT(a_t1v(k), e_t1v(k), om3, a_t2v(k), e_t2v(k), om3, (th3-1e-7), 'apo');
dv4v = [dv4v; dv4v(k)];
dt4v = [dt4v; dt4v(k)];

% manovra 4: ricircolarizzazione dell'orbita per intersecare l'orbita
% finale nell'apocentro
a5 = aF*(1+eF); % raggio dell'orbita circolare = raggio dell'apocentro finale
e5 = 0;
[dv5v(k), th5, dt5v(k)] = changeOrbitShapeALT(a_t2v(k), e_t2v(k), om3, a5, e5, om3, th4, 'per');
dv5v = [dv5v; dv5v(k)];
dt5v = [dt5v; dt5v(k)];

% utima manovra: decircolarizzazione
[dv6, th6, dt6] = changeOrbitShape(a5, e5, omF, aF, eF, omF, th5, 'apo');

% ultimo pezzo
dt7 = timeOfFlight(aF, eF, th6, thF,mu);

dt_totv(k) = sum(dt1) + dt2v(k) + dt3v(k) + dt4v(k) + dt5v(k) + sum(dt6) + dt7;
dv_totv(k) = sum(abs(dv1)) + dv2v(k) + dv3v(k) + dv4v(k) + dv5v(k) + sum(abs(dv6));

dv_totv = [dv_totv; dv_totv(k)];
dt_totv = [dt_totv; dt_totv(k)];
end

% stampo il dv ottimale
[dv_tot_opt, iter_v_opt] = min(dv_totv);

% plot dell'andamento di dvtot al variare di rb
figure(2)
plot(rbv, dv_totv(1:10000), 'r', 'LineWidth', 1)
xlabel('r_b [km]')
ylabel('dv [km/s]')

% PARAMETRI DELLA STRATEGIA OTTIMALE:
rb_opt = rbv(iter_v_opt);
a_t1 = a_t1v(iter_v_opt);
e_t1 = e_t1v(iter_v_opt);
dt2 = dt2v(iter_v_opt);
dv2 = dv2v(iter_v_opt);
dt3 = dt3v(iter_v_opt);
dv3 = dv3v(iter_v_opt);
a_t2 = a_t2v(iter_v_opt);
e_t2 = e_t2v(iter_v_opt);
dt4 = dt4v(iter_v_opt);
dv4 = dv4v(iter_v_opt);
dt5 = dt5v(iter_v_opt);
dv5 = dv5v(iter_v_opt);
dt_tot = dt_totv(iter_v_opt);

% PLOT E CALCOLI PARZIALI:


% 2) prima parte di biellittica
% aggiorno traettoria da th1 a punto di manovra th_in_be
[X,Y,Z] = plotOrbit(kep1,mu,(abs(th_in_be-th1)),precision);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

% disegno prima parte della biellittica:
kep_t1 = [a_t1, e_t1, iI, OMI, om_t1, th2];
[X,Y,Z] = plotOrbit(kep_t1,mu,360,precision);
orbit2 = plot3(X,Y,Z,'--m', 'LineWidth', 1);

% output prima parte di biellittica
fprintf('\n---- PRIMO IMPULSO BIELLITTICA ----\n')
stampInfoManovra(kep_t1,th_in_be, sum(abs(dt2)), sum(abs(dv2)))

% aggiorno traiettoria da th_in_be a th2
% aggiorno traettoria da th_in_be a punto di cambio piano th2
[X,Y,Z] = plotOrbit(kep_t1,mu,(abs(th2)),precision);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

% 3) CAMBIO PIANO
% aggiorno traettoria
[X,Y,Z] = plotOrbit(kep_t1,mu,abs(th3-th2),precision);
X_traj = X;
Y_traj = Y;
Z_traj = Z;

kep2 = [a_t1, e_t1, iF, OMF, om3, th3]; % orbita di trasferimento

% non plotto il cambio piano e basta, plotto direttamente quando ci sarà
% anche il cambio forma

% disegno il punto punto di manovra di cambio piano
[r2, v2] = kep2car(a_t1, e_t1, iF, OMF, om3, th3, mu);
pt_changePlane = plot3(r2(1), r2(2), r2(3), 'ok', 'LineWidth', 2);

% output cambio di piano
fprintf('\n---- CAMBIO DI PIANO ----\n')
stampInfoManovra(kep2, th2, dt3,dv3)

% disegno seconda parte della biellittica, nel nuovo piano:
kep_t2 = [a_t2, e_t2, iF, OMF, om3, th4];
[X,Y,Z] = plotOrbit(kep_t2,mu,360,precision);
orbit3 = plot3(X,Y,Z,'--b', 'LineWidth', 1);

% output seconda parte di biellittica
fprintf('\n---- SECONDO IMPULSO BIELLITTICA ----\n')
stampInfoManovra(kep_t2,th3, dt4, dv4)

kep5 = [a5, e5, iF, OMF, om3, th5];
%stampo punto di circolarizzazione
[r3, v3] = kep2car(a_t2, e_t2, iF, OMF, om3, th4, mu);
pt_circ2 = plot3(r3(1), r3(2), r3(3), 'oc', 'LineWidth', 1);
% disegno orbita circolare finale
[X,Y,Z] = plotOrbit(kep5,mu,360,precision);
orbit4 = plot3(X,Y,Z,'--k', 'LineWidth', 1);

% aggiorno traiettoria da th4 a th5:
dth5 = th5 - th4;
[X,Y,Z] = plotOrbit(kep_t2,mu,dth5,precision);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

% output seconda circolarizzazione
fprintf('\n---- SECONDA CIRCOLARIZZAZIONE ----\n')
stampInfoManovra(kep5,th4, dt5, dv5)

kep6 = [aF, eF, iF, OMF, omF, th6];

% stampo punto di decircolarizzazione
[r5, v5] = kep2car(a5, e5, iF, OMF, om3, th5, mu);
pt_5 = plot3(r5(1), r5(2), r5(3), 'oc', 'LineWidth', 1);

% disegno orbita circolare finale
[X,Y,Z] = plotOrbit(kep6,mu,360,precision);
orbit5 = plot3(X,Y,Z,'--y', 'LineWidth', 1);

% aggiorno traiettoria da th5 a th6:
dth5 = abs(th6 - th5);
[X,Y,Z] = plotOrbit(kep6,mu,dth5,precision);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

% output seconda circolarizzazione
fprintf('\n---- DECIRCOLARIZZAZIONE ----\n')
stampInfoManovra(kep6,th5, sum(dt6), sum(abs(dv6)))

% plot ultimo pezzo
if(thF < th6)
    dTh = 360 + thF - th6;
else
    dTh = thF - th6; %serve abs nel caso in cui th3 sia uguale a 0
end
[X,Y,Z] = plotOrbit(kep6, mu, dTh,precision);
X_traj = [X_traj; X];
Y_traj = [Y_traj; Y];
Z_traj = [Z_traj; Z];

fprintf('\nAttesa fino al punto finale: %4.2f s\n', dt7)

















%satellite
h = plot3(nan,nan,nan,"om", 'LineWidth',4);

% legend
% legend([orbitI start orbitF target in_be pt_circ1 orbit1 orbit2 pt_changePlane orbit3 pt_circ2 orbit4 pt_5 orbit5], 'Orbita Iniziale', 'Punto Iniziale', 'Orbita Finale', 'Punto Finale', 'P.to Inizio Biellittica', 'P.to Prima Circolarizzazione', 'Prima Parte Biellittica', 'P.to Cambio Piano e Cambio Forma', 'Seconda Parte Biellittica', 'P.to Seconda circolarizzazione', 'Seconda Orbita Circolare', 'P.to Decircolarizzazione', 'Orbita Finale decircolarizzata')
%output riassuntivo
fprintf('\n---- Riassunto strategia biellittica Circolarizzata ----\n')
fprintf('Raggio ottimale biellittica: %2.3f km\n', rb_opt)
fprintf('Variazione di velocità totale: %2.3f km/s\n', dv_tot_opt)
fprintf('Tempo totale di manovra:\n')
fprintf('Secondi: %5.2f s\n', dt_tot)
fprintf('Minuti: %5.2f h\n', dt_tot/60)
fprintf('Ore: %5.2f d\n', dt_tot/60/60)
fprintf('Giorni: %5.2f g\n', dt_tot/60/60/24)
