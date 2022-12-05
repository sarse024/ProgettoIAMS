%% strategia 2 impulsi
clear;
clc;
close all;

global mu;

mu = 398600;

% Ancora parecchio da sistemare
% Aggiungere condizione per evitare lo schianto a terra con tolleranza;
% Vedere se possibile scrivere una funzione
% Controllare le velocità
% Bug: se abbasso il meshdensity non converge la soluzione
% Aggiungere animazione

option = 'per'; %apo or per
precision = 0.1;
toll = 1e-4;

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

fprintf('IMPLEMENTAZIONE STRATEGIA 2 IMPULSI:\n\n')

%Trovo vettore che identifica l'intersezione
hI = cross(rI, vI);
hF = cross(rF, vF);
vector = -cross(hI,hF); %con il meno ottengo l'altra intersezione

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

%controllo
if(abs(r_point-r2) > toll)
    error('Errore di precisione')
end

punto_inter = plot3([0 r2(1)], [0, r2(2)], [0, r2(3)],'xr');


%% funzione 
tic
DIM = 1e3;
R_I = rotationMatrix(iI,OMI,omI);
p1 = R_I*rI/DIM;
p2 = R_I*r2/DIM;
a = 1/2*abs(norm(p2) - norm(p1));
[focalHyper,a, e] = def_conic(p1,p2,a);

fig = figure;
fp = fimplicit(focalHyper, [-DIM, DIM, -DIM, DIM]/10, 'MeshDensity',5000); %forse più elgante risolvere su a passando per Newton
focal = [fp.XData; fp.YData];
close(fig)
%%
orbit_shape = [];
dv = [];
t = [];
d1 = [];

for k = 1:length(focal)
    fc = focal(:, k);
    a_1 = 1/2* (norm(p1) - sqrt((p1(1) - fc(1))^2 + (p1(2) - fc(2))^2));
    a_2 = 1/2* (norm(p2) - sqrt((p2(1) - fc(1))^2 + (p2(2) - fc(2))^2));
    if(abs(a_1 - a_2) <= toll)
        a = a_1;
    else
        a = 1/2* (norm(p1) + sqrt((p1(1) - fc(1))^2 + (p1(2) - fc(2))^2));
    end

    [fun,a, e] = def_conic([0,0],fc,a);
    if(e>=0 && e<1)
        p = a*(1-e^2);
        om = thI - rad2deg(acos(1/e*(p/norm(p1)-1))) + omI;
        %om = thI - 360 + rad2deg(acos(1/e*(p/norm(p1)-1))) + omI; %trovare
        %questa condizione
        if(om < 0)
            om = 360 + om;
        end
        a = a*DIM;
        new_thI = thI + (omI - om);
        new_th_intersection = rad2deg(acos(1/e*(p/norm(p2)-1)));
        if(new_th_intersection < 0)
            new_th_intersection =  new_th_intersection + 360;
        end
        %calcoliamo le velocità rispetto all'orbita appena trovata
        [r1,v1] = kep2car(a,e,iI,OMI,om,new_thI,mu);
        [r2,v2] = kep2car(a,e,iI,OMI,om,new_th_intersection,mu);
        
        orbit_shape = [orbit_shape; a e om];
        dv1 = norm(v1 - vI);
        dv2 = norm(v_point - v2);
        d1 = [d1; dv1];
        dv = [dv; dv1 + dv2];
        t = [t; timeOfFlight(a,e,new_thI,new_th_intersection,mu)];
        
        %korbit = [orbit_shape(k,1), orbit_shape(k,2),iI, OMI, orbit_shape(k,3), 0];
        %[X,Y,Z] = plotOrbit(korbit,mu,360,0.1);
        %plot3(X,Y,Z)

    end

end
[vmin, k_Vmin] = min(dv);
[tmin, k_Tmin] = min(t);
kep_Vmin = [orbit_shape(k_Vmin,1), orbit_shape(k_Vmin,2),iI, OMI, orbit_shape(k_Vmin,3), 0];

figure(fig1)
hold on
a = orbit_shape(k_Vmin,1);
e = orbit_shape(k_Vmin,2);
om = orbit_shape(k_Vmin,3);
new_thI = thI + (omI - om);
if(new_thI < 0)
    new_thI = 360 + new_thI;
end
new_th_intersection = rad2deg(acos(1/e*(a*(1-e^2)/norm(r_point)-1)));
if(new_th_intersection < 0)
    new_th_intersection = 360+ new_th_intersection;
end
[rprov1, vprov1] = kep2car(a, e,iI, OMI, om, new_thI,mu);
punto = plot3(rprov1(1),rprov1(2),rprov1(3),'xm');

[X,Y,Z] = plotOrbit(kep_Vmin,mu,360,0.1);
orbitVmin = plot3(X,Y,Z);

fprintf('%%%% MANOVRA CON V_tot MINIMA %%%%\n')
fprintf('Costo Totale della manovra: %2.4f\n', vmin)
fprintf('Tempo Totale di manovra: %4.2f\n', t(k_Vmin))
fprintf('Prima manovra:\n')
fprintf('Punto di manovra:\n')
fprintf('Velocita: %2.4f\n', norm(vI-vprov1))
fprintf('Angolo di velocità: \n') %Da calcolare
fprintf(['Parametri Orbitali:\n' ...
        '\ta\t\t\te\t\ti\t\tOM\t\tom\t\n ' ...
        '%4.2f\t %.4f\t %2.2f\t %2.2f\t %2.2f\n'], kep_Vmin(1), kep_Vmin(2), kep_Vmin(3),kep_Vmin(4),kep_Vmin(5))
fprintf('Errore sul primo punto di manovra: %2.2f m\n', norm(rprov1-rI))


[rprov2, vprov2] = kep2car(a, e,iI, OMI, om, new_th_intersection,mu);
prova2 = plot3(rprov2(1),rprov2(2),rprov2(3),'xm', 'LineWidth',5);
fprintf('Seconda manovra:\n')
fprintf('Punto di manovra:\n')
fprintf('Velocita: %2.4f\n', norm(v_point-vprov2))
fprintf('Angolo di velocità: \n') %Da calcolare
fprintf('Errore sul secondo punto di manovra: %2.2f m\n', norm(rprov2-r_point))
fprintf('Time of flight: %4.2f\n', timeOfFlight(a,e,new_thI,new_th_intersection,mu))

toc


legend([punto, orbitVmin, prova2, punto_inter], 'Errore', 'Orbita', 'Errore2', 'Punto di Manovra 2')



