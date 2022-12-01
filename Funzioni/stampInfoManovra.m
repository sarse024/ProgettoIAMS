function stampInfoManovra(kepT, th_manovra, t,dv)
    fprintf('Punto di manovra: %2.2f gradi\n', th_manovra)
    fprintf('Costo della manovra: %2.3f km/s\n', dv)
    fprintf('Tempo di manovra: %4.1f s\n', t)
    fprintf(['Nuova orbita dopo manovra:\n' ...
        '\ta\t\t\te\t\ti\t\tOM\t\tom\t\tth\n ' ...
        '%4.2f\t %.4f\t %2.2f\t %2.2f\t %2.2f\t %2.2f\n'], kepT(1), kepT(2), kepT(3),kepT(4),kepT(5), kepT(6))
end