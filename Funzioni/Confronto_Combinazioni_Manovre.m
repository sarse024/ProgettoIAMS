%% quale combinazione di manovre risulta essere quella più efficiente?
%% strategia base
clear;
clc;
close all;

% Costant
global mu;

mu = 398600;

option = 'per'; % useful for changing start point in changeShape
precision=0.1;

% Initialize main figure
figure
Terra3d
hold on
grid on

%%%%%%%% GRUPPO B7 %%%%%%%%
dati_elaborati = [5088.9118 -3196.5659 -8222.7989 1.9090 5.6220 -1.0700 14020.0000 0.3576 1.3220 0.9764 1.8130 0.4336];

% --- Initial Orbit ---
rI = dati_elaborati(1:3)';
vI = dati_elaborati(4:6)';

% Orbital parametre of initial orbit
[aI, eI, iI, OMI, omI, thI] = car2kep(rI, vI, mu);
kepI = [aI,eI,iI,OMI,omI,thI];

% plot initial orbit
[X,Y,Z] = plotOrbit(kepI,mu,360,precision);
orbitI = plot3(X,Y,Z, 'LineWidth', 2);
% plot starting point
start = plot3(rI(1), rI(2), rI(3), 'xb', 'LineWidth', 4);

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

% plot final orbit
[X,Y,Z] = plotOrbit(kepF,mu,360,precision);
orbitF = plot3(X,Y,Z, 'LineWidth', 2);
target = plot3(rF(1), rF(2), rF(3), 'xr', 'LineWidth', 4);

%plot apocenter and pericenter first and last orbit 
%black for pericenters and cyan for apocenters
[rp1, vp1] = kep2car(aI,eI,iI,OMI,omI,0,mu);
perI=plot3(rp1(1),rp1(2),rp1(3),'ok','Linewidth',4); %pericenter first orbit
[ra1, va1] = kep2car(aI,eI,iI,OMI,omI,180,mu);
apoI=plot3(ra1(1),ra1(2),ra1(3), 'oc','Linewidth',4); %apocenter first orbit

[rp2, vp2] = kep2car(aF,eF,iF,OMF,omF,0,mu);
perF=plot3(rp2(1),rp2(2),rp2(3),'ok','Linewidth',4); %pericenter second orbit
[ra2, va2] = kep2car(aF,eF,iF,OMF,omF,180,mu);
apoF=plot3(ra2(1),ra2(2),ra2(3), 'oc','Linewidth',4); %apocenter second orbit

dv_total=[];
dt_total=[];
%% base strategy
% --- 1° Manoeuvre: change plane --- 
[dv1, om1, th1, dt1] = changeOrbitalPlane(aI, eI, iI, OMI, omI, iF, OMF, thI);

% --- 2° Manoeuvre: change periapsisArg --- 
[dv2, omPR, vec_th, dt2] = changePeriapsisArg(aI,eI,om1, (omF-om1), th1);

% Check if new om is correct (equal to omF)
if(omPR ~= omF)
    error('Error in the code')
end

th21 = vec_th(1); % th before the maneuver
th2 = vec_th(2); % th after the change of velocity

% Pericenter-to-apocenter transfer
[dv3, th3, dt3] = changeOrbitShape(aI, eI, omF, aF, eF, omF, th2, option);

% time from th3 to thF on the final orbit
dt4 = timeOfFlight(aF, eF, th3, thF,mu);

dt_tot_base = dt1 + dt2 + sum(dt3) + dt4;
dv_tot_base = dv1 + dv2 + sum(abs(dv3));

dv_total=[dv_total;dv_tot_base];
dt_total=[dt_total;dt_tot_base];

%% alternative 1 (swap manoeuvre: change shape at pericenter first and then change orbital plane and periapsisArg)

% --- 1° Manoeuvre: change plane --- pericenter-apocenter
[dv1_1, th1_1, dt1_1] = changeOrbitShape(aI, eI, omI, aF, eF, omI, thI, option);

% --- 2° Manoeuvre: change orbitalPlane --- 
[dv2_1, om1_1, th2_1, dt2_1] = changeOrbitalPlane(aF, eF, iI, OMI, omI, iF, OMF, th1_1);

% --- 3° Manoeuvre: change periapsisArg --- 
[dv3_1, omPR_1, vec_th_1, dt3_1] = changePeriapsisArg(aF,eF,om1_1, (omF-om1_1), th2_1); 

% Check if new om is correct (equal to omF)
if(omPR ~= omF) 
    error('Errore nei calcoli')
end

th3_1=vec_th_1(2);% th after the change of velocity

% time from th3_1 to thF on the final orbit
dt4_1 = timeOfFlight(aF, eF, th3_1, thF, mu);

%sum of partials dv and dt to have the total values of dv and dt
dv_tot_a1=sum(abs(dv1_1))+abs(dv2_1)+abs(dv3_1);
dt_tot_a1=sum(dt1_1)+dt2_1+dt3_1+dt4_1;

%insert these value in two vectors to compare with the values obtained from
%other strategies
dv_total=[dv_total;dv_tot_a1];
dt_total=[dt_total;dt_tot_a1];

%% alternative 2 (change shape at apocenter first and then change orbital plane and periapsisArg)

% --- 1° Manoeuvre: change shape --- apocenter-pericenter
[dv1_2, th1_2, dt1_2] = changeOrbitShape(aI, eI, omI, aF, eF, omI, thI, 'apo');

% --- 2° Manoeuvre: change orbital plane ---
[dv2_2, om1_2, th2_2, dt2_2] = changeOrbitalPlane(aF, eF, iI, OMI, omI, iF, OMF, th1_2);

% --- 3° Manoeuvre: change periapsisArg ---
[dv3_2, omPR_2, vec_th_2, dt3_2] = changePeriapsisArg(aF,eF,om1_2, (omF-om1_2), th2_2); 

% Check if new om is correct (equal to omF)
if(omPR ~= omF) 
    error('Errore nei calcoli')
end

th3_2=vec_th_2(2);% th after the change of velocity

% time from th3_2 to thF on the final orbit
dt4_2 = timeOfFlight(aF, eF, th3_2, thF, mu);

%sum of partials dv and dt to have the total values of dv and dt
dv_tot_a2=sum(abs(dv1_2))+abs(dv2_2)+abs(dv3_2);
dt_tot_a2=sum(dt1_2)+dt2_2+dt3_2+dt4_2;

%insert these value in two vectors to compare with the values obtained from
%other strategies
dv_total=[dv_total;dv_tot_a2];
dt_total=[dt_total;dt_tot_a2];

%% alternative 3 (change orbital plane, then change orbit shape and then change periapsis arg)
% --- 1° Manoeuvre: change plane --- 
[dv1_3, om1_3, th1_3, dt1_3] = changeOrbitalPlane(aI, eI, iI, OMI, omI, iF, OMF, thI);

% --- 2° Manoeuvre: change orbitShape --- pericenter-apocenter
[dv2_3, th2_3, dt2_3] = changeOrbitShape(aI, eI, om1_3, aF, eF, om1_3, th1_3, option);

% --- 3° Manoeuvre: change periapsisArg --- 
[dv3_3, omPR_3, vec_th_3, dt3_3] = changePeriapsisArg(aF,eF,om1_3, (omF-om1_3), th2_3);

% Check if new om is correct (equal to omF)
if(omPR ~= omF)
    error('Error in the code')
end

th3_3=vec_th_3(2); % th after the change of velocity

% time from th3_3 to thF on the final orbit
dt4_3 = timeOfFlight(aF, eF, th3_3, thF, mu);

%sum of partials dv and dt to have the total values of dv and dt
dv_tot_a3=abs(dv1_3)+sum(abs(dv2_3))+abs(dv3_3);
dt_tot_a3=dt1_3+sum(dt2_3)+dt3_3+dt4_3;

%insert these value in two vectors to compare with the values obtained from
%other strategies
dv_total=[dv_total;dv_tot_a3];
dt_total=[dt_total;dt_tot_a3];

%% alternativa 4 (change orbital plane, then go to the descendent node and then make the other manoeuvres)

% --- 1° Manoeuvre: change plane --- pericenter-apocenter
[dv1_4, om1_4, th1_4, dt1_4] = changeOrbitalPlane(aI, eI, iI, OMI, omI, iF, OMF, thI);
th2_4=th1_4+180;
if th2_4>360
    th2_4=th2_4-360;
end

% time from th1_4 to th2_4 on the first transfer orbit
dt2_4=timeOfFlight(aI,eI,th1_4,th2_4,mu);

% --- 2° Manoeuvre: change periapsisArg --- 
[dv3_4, omPR_4, vec_th_3_4, dt3_4] = changePeriapsisArg(aI,eI,om1_4, (omF-om1_4), th2_4);

% Check if new om is correct (equal to omF)
if(omPR ~= omF)
    error('Error in the code')
end

th3_4 = vec_th_3_4(2); % th after the change of velocity

% --- 3° Manoeuvre: change orbit shape --- Pericenter-to-apocenter transfer
[dv4_4, th4_4, dt4_4] = changeOrbitShape(aI, eI, omF, aF, eF, omF, th3_4, option);

% time from th3 to thF on the final orbit
dt5_4 = timeOfFlight(aF, eF, th4_4, thF,mu);

%sum of partials dv and dt to have the total values of dv and dt
dv_tot_a4=abs(dv1_4)+abs(dv3_4)+sum(abs(dv4_4));
dt_tot_a4=dt1_4+dt2_4+dt3_4+sum(dt4_4)+dt5_4;

%insert these value in two vectors to compare with the values obtained from
%other strategies
dv_total=[dv_total;dv_tot_a4];
dt_total=[dt_total;dt_tot_a4];

%% alternative 5 (change orbital plane first then change periapsis arg and then change orbital shape at apocenter)
% --- 1° Manoeuvre: change orbital plane --- 
[dv1_5, om1_5, th1_5, dt1_5] = changeOrbitalPlane(aI, eI, iI, OMI, omI, iF, OMF, thI);

% --- 2° Manoeuvre: change periapsisArg --- 
[dv2_5, omPR5, vec_th_5, dt2_5] = changePeriapsisArg(aI,eI,om1_5, (omF-om1_5), th1_5);

% Check if new om is correct (equal to omF)
if(omPR ~= omF)
    error('Error in the code')
end

th2_5 = vec_th_5(2); % th after the change of velocity

%  --- 2° Manoeuvre: change orbit shape --- apocenter-pericenter transfer
[dv3_5, th3_5, dt3_5] = changeOrbitShape(aI, eI, omF, aF, eF, omF, th2_5, 'apo');

% time from th3 to thF on the final orbit
dt4_5 = timeOfFlight(aF, eF, th3_5, thF,mu);

dv_tot_a5 = dv1_5 + dv2_5 + sum(abs(dv3_5));
dt_tot_a5 = dt1_5 + dt2_5 + sum(dt3_5) + dt4_5;

dv_total=[dv_total;dv_tot_a5];
dt_total=[dt_total;dt_tot_a5];

%% alternative 6 (change orbital plane, then change orbit shape at apocenter and then change periapsis arg)
% --- 1° Manoeuvre: change plane --- 
[dv1_6, om1_6, th1_6, dt1_6] = changeOrbitalPlane(aI, eI, iI, OMI, omI, iF, OMF, thI);

% --- 2° Manoeuvre: change orbitShape --- apocenter-pericenter
[dv2_6, th2_6, dt2_6] = changeOrbitShape(aI, eI, om1_3, aF, eF, om1_6, th1_6, 'apo');

% --- 3° Manoeuvre: change periapsisArg --- 
[dv3_6, omPR_6, vec_th_6, dt3_6] = changePeriapsisArg(aF,eF,om1_6, (omF-om1_6), th2_6);

% Check if new om is correct (equal to omF)
if(omPR ~= omF)
    error('Error in the code')
end

th3_6=vec_th_6(2); % th after the change of velocity

% time from th3_3 to thF on the final orbit
dt4_6 = timeOfFlight(aF, eF, th3_6, thF, mu);

%sum of partials dv and dt to have the total values of dv and dt
dv_tot_a6=abs(dv1_6)+sum(abs(dv2_6))+abs(dv3_6);
dt_tot_a6=dt1_6+sum(dt2_6)+dt3_6+dt4_6;

%insert these value in two vectors to compare with the values obtained from
%other strategies
dv_total=[dv_total;dv_tot_a6];
dt_total=[dt_total;dt_tot_a6];

%% alternativa 7 (change orbital plane, then go to the descendent node 
% and then make the other manoeuvres with change orbit shape at apocenter)

% --- 1° Manoeuvre: change plane --- pericenter-apocenter
[dv1_7, om1_7, th1_7, dt1_7] = changeOrbitalPlane(aI, eI, iI, OMI, omI, iF, OMF, thI);
th2_7=th1_7+180;
if th2_7>360
    th2_4=th2_7-360;
end

% time from th1_7 to th2_7 on the first transfer orbit
dt2_7=timeOfFlight(aI,eI,th1_7,th2_7,mu);

% --- 2° Manoeuvre: change periapsisArg --- 
[dv3_7, omPR_7, vec_th_3_7, dt3_7] = changePeriapsisArg(aI,eI,om1_7, (omF-om1_7), th2_7);

% Check if new om is correct (equal to omF)
if(omPR ~= omF)
    error('Error in the code')
end

th3_7 = vec_th_3_7(2); % th after the change of velocity

% --- 3° Manoeuvre: change orbit shape --- Pericenter-to-apocenter transfer
[dv4_7, th4_7, dt4_7] = changeOrbitShape(aI, eI, omF, aF, eF, omF, th3_7, 'apo');

% time from th3 to thF on the final orbit
dt5_7 = timeOfFlight(aF, eF, th4_7, thF,mu);

%sum of partials dv and dt to have the total values of dv and dt
dv_tot_a7=abs(dv1_7)+abs(dv3_7)+sum(abs(dv4_7));
dt_tot_a7=dt1_7+dt2_7+dt3_7+sum(dt4_7)+dt5_7;

%insert these value in two vectors to compare with the values obtained from
%other strategies
dv_total=[dv_total;dv_tot_a7];
dt_total=[dt_total;dt_tot_a7];
