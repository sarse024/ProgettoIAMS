%% strategia 2 impulsi
clear all;
close all;
clc;

global mu;

mu = 398600;
ALLERT = 6378 + 300; %Km

% Ancora parecchio da sistemare
% Controllare le velocità
% Aggiungere animazione
% Abbellire grafici
% Sistemare output
% Pulire Codice da cose inutili e non necessarie
% Sistemare def_conic (Capire bene gradiente per definire om)

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

fprintf('IMPLEMENTAZIONE STRATEGIA 2 IMPULSI DIRETTA:\n\n')

hI = cross(rI, vI);
hF = cross(rF, vF);

if(strcmp(option,'asc'))
    vector = cross(hI,hF); %con il meno ottengo l'altra intersezione
else
    vector = -cross(hI,hF);
end

R_F = rotationMatrix(iF,OMF,omF);

vec = R_F*vector; %vettore della posizione nel sistema perifocae
th_intersection = atan2(vec(2),vec(1)); %trovo il theta

if(th_intersection < 0)
    th_intersection = 2*pi+th_intersection;
end

r_inters = aF*(1-eF^2)/(1+eF*cos(th_intersection)); %calcolo il raggio di intersezione
t = r_inters/norm(vector); %il raggio appartiene alla retta di intersezione tra i due piani
r2 = vector*t;

[r_point,v_point]  = kep2car(aF,eF,iF,OMF,omF, rad2deg(th_intersection),mu);
point_man = plot3(r_point(1),r_point(2),r_point(3),'xm');
%controllo
if(abs(r_point-r2) > toll)
    error('Errore di precisione')
end

th_initial = 0:359;
%th_initial = mod(th_initial,360);
%%
risult = [th_initial', zeros(length(th_initial),1), zeros(length(th_initial),1)]; %[th_initial, delta_v, delta_t]

%%
tic
wb = waitbar(0,'Calcolo Orbite');
for k = 1:length(th_initial)
    r1 = kep2car(aI,eI,iI,OMI,omI,th_initial(k),mu);
    [kepT, dv, th] = findOrbit(r1,r2,kepI,kepF,option);
    risult(k,2) = dv;

    %[X,Y,Z] = plotOrbit(kepT,mu,360,0.1);
    %plot3(X,Y,Z);
    
    t1 = timeOfFlight(aI,eI,thI,th(1),mu);
    t2 = timeOfFlight(kepT(1),kepT(2),th(2),th(3),mu);
    t3 = timeOfFlight(aF,eF,th(4),thF,mu);
    risult(k,3) = t1 + t2 + t3;
    waitbar(k/length(th_initial), wb, sprintf('Calcolo Orbite: %2.0f %%', k/length(th_initial)*100));
end
delete(wb)
toc
%pulisco la matrice dei risultati 
%%
risultTest = risult;
%%
risult = risultTest;

for k = length(risult):-1:1
    if(risult(k,2) == 0)
        risult(k,:) = []; %cancello la riga
        k = k;
    end
end
%%
figure(fig1)
[vmin, kvmin] = min(risult(:,2));
[r1,v1] = kep2car(aI,eI,iI,OMI,omI,risult(kvmin,1),mu);
[kepT, dv, th] = findOrbit(r1,r2,kepI,kepF,option);
[X,Y,Z] = plotOrbit(kepT,mu,360,0.1);
plot3(X,Y,Z);

%%
figure
tiledlayout(2,1)
nexttile
plot(risult(:,1),risult(:,2),'-o');
title('Confronto velocità optimizzata')
nexttile
plot(risult(:,1),risult(:,3)/60/60,'-o');
title('Confronto tempo optimizzato')

%% test grafici
figure(fig1)
a = kepT(1);
e = kepT(2);
om = norm(kepT(5));
th1 = th(1);
new_th1 = th(2);
new_th2 = th(3);
th2 = th(4);
[rprov1, vprov1] = kep2car(a, e,iI, OMI, om, new_th1,mu);
punto = plot3(rprov1(1),rprov1(2),rprov1(3),'xm');

fprintf('%%%% MANOVRA CON V_tot MINIMA %%%%\n')
fprintf('Costo Totale della manovra: %2.4f\n', vmin)
fprintf('Tempo Totale di manovra: %4.2f\n', risult(kvmin,3))
fprintf('Prima manovra:%3.2f\n', th1)
fprintf('Punto di manovra:\n')
fprintf('Velocita: %2.4f\n', norm(vprov1-v1))
fprintf('Angolo di velocità: \n') %Da calcolare
fprintf(['Parametri Orbitali:\n' ...
        '\ta\t\t\te\t\ti\t\tOM\t\tom\t\n ' ...
        '%4.2f\t %.4f\t %2.2f\t %2.2f\t %2.2f\n'], kepT(1),kepT(2),kepT(3),kepT(4),kepT(5))
fprintf('Errore sul primo punto di manovra: %2.2f m\n', norm(rprov1-r1))


[rprov2, vprov2] = kep2car(a, e, kepT(3), kepT(4), om, new_th2,mu);
prova2 = plot3(rprov2(1),rprov2(2),rprov2(3),'xm', 'LineWidth',5);
fprintf('Seconda manovra:\n')
fprintf('Punto di manovra:%f\n',new_th2)
fprintf('Velocita: %2.4f\n', norm(v_point-vprov2))
fprintf('Angolo di velocità: \n') %Da calcolare
fprintf('Errore sul secondo punto di manovra: %2.2f m\n', norm(r_point-rprov2))
fprintf('Time of flight: %4.2f\n', timeOfFlight(a,e,new_th1,new_th2,mu))



