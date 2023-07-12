function stampInfoManovra(kepT, th_manovra, t,dv, t_tot)
% stampInfoManovra.m - Function to print output of the maneuver
%
% PROTOTYPE:
% stampInfoManovra(kepT, th_manovra, t,dv)
%
% DESCRIPTION:
% Function to print the most important information of the maneuver:
% - Point of maneuver
% - Cost in delta_v
% - Time from point A to the point of maneuvre
% - New orbital parametre after the impulse
%
% INPUT:
% kepT          [1x6]   Orbital Parametre                                       [a,e,i,OM,om,th]
% th_manovra    [1x1]   True anomaly of the point of maneuver on kepT orbit     [deg]
% t             [1x1]   Time of flight of the maneuver                          [s]
% dv            [1x1]   Speed difference                                        [km/s]
% t             [1x1]   Total Time of flight                                    [s]

    fprintf('Punto di manovra: %2.2f gradi\n', th_manovra)
    fprintf('Costo della manovra: %2.3f m/s\n', dv*1000)
    fprintf('Tempo totale di manovra: %5.2f s\n', t_tot)
    fprintf('Tempo di manovra:\n')
    fprintf('- Secondi: %5.2f s\n', t)
    fprintf('- Minuti: %5.2f m\n', t/60)
    fprintf('- Ore: %5.2f h\n', t/60/60)
    fprintf(['Nuova orbita dopo manovra:\n' ...
        '\ta\t\t\te\t\ti\t\tOM\t\tom\t\tth\n ' ...
        '%4.2f\t %.4f\t %2.2f\t %2.2f\t %2.2f\t %2.2f\n'], kepT(1), kepT(2), kepT(3),kepT(4),kepT(5), kepT(6))
end