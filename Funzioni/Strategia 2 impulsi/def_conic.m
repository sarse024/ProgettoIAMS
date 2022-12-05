function [f_conic,a,e] = def_conic(f1,f2,a)
    
    A = 16*a^2-4*(f1(1) - f2(1))^2;
    B = -8*(f1(1) - f2(1))*(f1(2) - f2(2));
    C = 16*a^2-4*(f1(2) - f2(2))^2;
    D = 4*(f1(1) - f2(1))*(f1(1)^2 - f2(1)^2 + f1(2)^2 - f2(2)^2) - 16*(a^2)*(f1(1) + f2(1));
    E = 4*(f1(2) - f2(2))*(f1(1)^2 - f2(1)^2 + f1(2)^2 - f2(2)^2) - 16*(a^2)*(f1(2) + f2(2));
    F = 4*(f1(1)^2 + f1(2)^2)*(f2(1)^2 + f2(2)^2) - (f1(1)^2 + f2(1)^2 + f1(2)^2 + f2(2)^2 - 4*a^2)^2;
    
    f_conic = @(x,y) A*x.^2 + B*x.*y + C*y.^2 + D*x + E*y + F;
    
    matrix = [A B/2 D/2; B/2 C E/2; D/2 E/2 F];
    Q = [A B D; B C E; D E F];
    det(Q);
    if(det(Q) ~= 0) %non degenere
        if(det(matrix)>0)
            eta = -1;
        else
            eta = 1;
        end
        e = sqrt(2*sqrt((A-C)^2 + B^2)/(eta*(A+C) + sqrt((A-C)^2+B^2)));
    else
        e = -1; %test per indicare conica degenere
    end

    if(B^2-4*A*C>0) %caso iperbole
        a = -a;
    end

end