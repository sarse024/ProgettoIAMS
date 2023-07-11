% Code to compute the optimal bielliptic radius in terms of minimum dv
clear;
clc;
close all;
% defining parameters for the study:
global mu;
mu = 398600;

%%%%%%%% B7 DATA %%%%%%%%
dati_elaborati = [5088.9118 -3196.5659 -8222.7989 1.9090 5.6220 -1.0700 14020.0000 0.3576 1.3220 0.9764 1.8130 0.4336];

% Initial Orbit Data
rI = dati_elaborati(1:3)';
vI = dati_elaborati(4:6)';

% Other parameters of the initial orbit:
[aI, eI, iI, OMI, omI, thI] = car2kep(rI, vI, mu);
kepI = [aI,eI,iI,OMI,omI,thI];

% Final Orbit Data:
aF = dati_elaborati(7);
eF = dati_elaborati(8);
iF = rad2deg(dati_elaborati(9));
OMF = rad2deg(dati_elaborati(10));
omF = rad2deg(dati_elaborati(11));
thF = rad2deg(dati_elaborati(12));

% Other parameters of the final orbit:
[rF, vF] = kep2car(aF,eF,iF,OMF,omF,thF,mu);
kepF = [aF,eF,iF,OMF,omF,thF];

fprintf('IMPLEMENTAZIONE STRATEGIA BIELLITTICA 2:\n\n')

fprintf(['Parametri orbitali punto di partenza:\n' ...
        '\ta\t\t\te\t\ti\t\tOM\t\tom\t\tth\n ' ...
        '%4.2f\t %.4f\t %2.2f\t %2.2f\t %2.2f\t %2.2f\n\n'], kepI(1), kepI(2), kepI(3),kepI(4),kepI(5), kepI(6))
fprintf(['Parametri orbitali punto di arrivo:\n' ...
        '\ta\t\t\te\t\ti\t\tOM\t\tom\t\tth\n ' ...
        '%4.2f\t %.4f\t %2.2f\t %2.2f\t %2.2f\t %2.2f\n'], kepF(1), kepF(2), kepF(3),kepF(4),kepF(5), kepF(6))

% Computing the angle of orbital plane change in relation to the initial
% orbit: this is done to give the impulse to start the bielliptic strategy
% 180° before this angle. This allows the manuever of orbital plane change
% to be done at the apocentre of the elliptic orbit.
[~, ~, th_cp, ~] = changeOrbitalPlane(aI, eI, iI, OMI, omI, iF, OMF, thI);
% Starting point of the bielliptic (with plot):
th_in_be = th_cp - 180;
[r_in_be, v_in_be] = kep2car(aI, eI, iI, OMI, omI, th_in_be, mu);

% First maneuver: Circularization of the initial orbit (at the apocentre)
raI = aI*(1+eI);
a1 = raI;
e1 = 0;
[dv1, th1, dt1] = changeOrbitShape(aI, eI, omI, a1, e1, omI, thI, 'apo');
thman1 = 180;
dth1 = thman1 - thI;
kep1 = [a1, e1, iI, OMI, omI, th1];

fprintf('\n---- CIRCOLARIZZAZIONE 1 ----\n')
stampInfoManovra(kep1,thman1, sum(abs(dt1)), sum(abs(dv1)))

% Defining vectors used in the 'for' cicle to optimize the strategy
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

% Optimization
for k = 1 : 10000
rb = rbv(k);
rp_t1 = a1; % since rp_t1 = r1
ra_t1v(k) = rb;
a_t1v(k) = (rp_t1+ra_t1v(k))/2;
e_t1v(k) = (ra_t1v(k)-rp_t1)/(ra_t1v(k)+rp_t1);
a_t1v = [a_t1v; a_t1v(k)];
e_t1v = [e_t1v; e_t1v(k)];
ra_t1v = [ra_t1v; ra_t1v(k)];

om_t1 = omI + th_in_be;

% Second maneuver: First half of the bielliptic transfer
[dv2v(k), th2, dt2v(k)] = changeOrbitShapeALT(a1, e1, om_t1, a_t1v(k), e_t1v(k), om_t1, th_in_be, 'per');
dv2v = [dv2v; dv2v(k)];
dt2v = [dt2v; dt2v(k)];

% At th2: change of orbital plane and second half of the bielliptic
% Third maneuver: change of plane
[dv3v(k), om3, th3, dt3v(k)] = changeOrbitalPlane(a_t1v(k), e_t1v(k), iI, OMI, om_t1, iF, OMF, (th2-1e-7)); 
dv3v = [dv3v; dv3v(k)];
dt3v = [dt3v; dt3v(k)];

% Fourth maneuver: second half of the bielliptic
ra3v(k) = rb;
rp3 = aF*(1+eF);
a_t2v(k) = (ra3v(k) + rp3)/2;
e_t2v(k) = (ra3v(k)-rp3)/(ra3v(k)+rp3);

[dv4v(k), th4, dt4v(k)] = changeOrbitShapeALT(a_t1v(k), e_t1v(k), om3, a_t2v(k), e_t2v(k), om3, (th3-1e-7), 'apo');
dv4v = [dv4v; dv4v(k)];
dt4v = [dt4v; dt4v(k)];

% Fifth maneuver: ricircularization of the orbit to intersect the final
% orbit in its apocentre
a5 = aF*(1+eF); % raggio dell'orbita circolare = raggio dell'apocentro finale
e5 = 0;
[dv5v(k), th5, dt5v(k)] = changeOrbitShapeALT(a_t2v(k), e_t2v(k), om3, a5, e5, om3, th4, 'per');
dv5v = [dv5v; dv5v(k)];
dt5v = [dt5v; dt5v(k)];

% Sixth maneuver: decircularization into the final orbit
[dv6, th6, dt6] = changeOrbitShape(a5, e5, omF, aF, eF, omF, th5, 'apo');

% Ultimately, we reach the target point
dt7 = timeOfFlight(aF, eF, th6, thF,mu);

% Printing the total time and cost of the strategy
dt_totv(k) = sum(dt1) + dt2v(k) + dt3v(k) + dt4v(k) + dt5v(k) + sum(dt6) + dt7;
dv_totv(k) = sum(abs(dv1)) + dv2v(k) + dv3v(k) + dv4v(k) + dv5v(k) + sum(abs(dv6));

dv_totv = [dv_totv; dv_totv(k)];
dt_totv = [dt_totv; dt_totv(k)];
end

% Printing the optimal dv_tot and its iteration:
[dv_tot_opt, iter_v_opt] = min(dv_totv);

% Parameters of the optimal strategy:
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

% Final output with all the results
fprintf('\n---- Riassunto strategia biellittica Circolarizzata ----\n')
fprintf('Raggio ottimale biellittica: %2.3f km\n', rb_opt)
fprintf('Variazione di velocità totale: %2.3f km/s\n', dv_tot_opt)
fprintf('Tempo totale di manovra:\n')
fprintf('Secondi: %5.2f s\n', dt_tot)
fprintf('Minuti: %5.2f min\n', dt_tot/60)
fprintf('Ore: %5.2f h\n', dt_tot/60/60)
fprintf('Giorni: %5.2f d\n', dt_tot/60/60/24)


%% Dv_tot evolution with rb: plot
figure;
plot(rbv, dv_totv(1:10000), 'o', 'MarkerSize', 1);
title('Total velocity variation with bielliptic radius', 'FontSize', 16)
xlabel('Bielliptic radius [km]', 'FontSize', 12)
ylabel('Total velocity variation [km/s]', 'FontSize', 12)
