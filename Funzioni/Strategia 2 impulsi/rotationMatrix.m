function R = rotationMatrix(i,OM,om)
    
    i = deg2rad(i);
    OM = deg2rad(OM);
    om = deg2rad(om);
    
    R3_OM = [cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];
    R1_i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
    R3_om = [cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];
    
    R = R3_om*R1_i*R3_OM; %matrice di rotazione
end