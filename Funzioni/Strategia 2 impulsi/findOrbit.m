function [kepT, dv, th] = findOrbit(r1, r2, kep1, kep2, option)
    %voption da togliere
    %funziona anche per la direttissima
    %necessario implementare un controllo su r1 e r2 nel cicloS
    ALLERT = 6378 + 300; %Km
    DIM = 1e3;
    toll = 1e-2;
    toll2 = 10;

    a1 = kep1(1);
    e1 = kep1(2);
    i1 = kep1(3);
    OM1 = kep1(4);
    om1 = kep1(5);
    p1 = a1*(1-e1^2);
    
    a2 = kep2(1);
    e2 = kep2(2);
    i2 = kep2(3);
    OM2 = kep2(4);
    om2 = kep2(5);
    p2 = a2*(1-e2^2);

    mu = 398600;
    if(strcmp(option,'asc'))
        N = cross(r1,r2); %vettore normale al piano
    elseif(strcmp(option,'disc'))
        N = cross(r2,r1);
    else
        error('\nInput sbagliato')
    end

    N = cross(r1,r2);
    if(N(3)<0)
        N = cross(r2,r1);
    end

    iT = rad2deg(acos(N(end)/norm(N))); %inclinazione del piano

    k = [0 0 1];
    n = cross(k,N);

    %RAAN (Ascensione retta del nodo ascendente)
    if n(2)>=0
        OMT=rad2deg(acos(n(1)/norm(n)));
    else
        OMT=rad2deg(2*pi-acos(n(1)/norm(n)));
    end

    th1 = rad2deg(acos(1/e1*(p1/norm(r1)-1)));
    [new_r1,v1] = kep2car(a1,e1,i1,OM1,om1,th1,mu);
    if(norm(new_r1 - r1) > toll2)
        th1 = 360 - th1;
        [~,v1] = kep2car(a1,e1,i1,OM1,om1,th1,mu);
    end

    th2 = rad2deg(acos(1/e2*(p2/norm(r2)-1))); %sto escludendo una soluzione  
    [new_r2,v2] = kep2car(a2,e2,i2,OM2,om2, th2,mu);
    if(norm(new_r2 - r2) > toll2)
        th2 = 360 - th2;
        [~,v2] = kep2car(a2,e2,i2,OM2,om2,th2,mu);
    end

    R_I = rotationMatrix(iT,OMT,0); %definisco la rotazione rispetto ad un piano
    p1 = R_I*r1/DIM;
    p2 = R_I*r2/DIM;

    %definisco il luogo geometrico delle posizioni del secondo fuoco
    a = 1/2*abs(norm(p2) - norm(p1));
    focalHyper = def_conic(p1,p2,a);
    
    fig = figure;
    fp = fimplicit(focalHyper, [-DIM, DIM, -DIM, DIM]/DIM*100, 'MeshDensity',5000); %forse più elgante risolvere su a passando per Newton
    focal = [fp.XData; fp.YData];
    close(fig)

    orbit_shape = [];
    dv = [];
    th = [];

    for k = 1:length(focal)
        fc = focal(:, k);
        a_1 = 1/2* (norm(p1) - sqrt((p1(1) - fc(1))^2 + (p1(2) - fc(2))^2));
        a_2 = 1/2* (norm(p2) + sqrt((p2(1) - fc(1))^2 + (p2(2) - fc(2))^2));
        if(abs(a_1 - a_2) <= toll)
            a = a_1;
        else
            a = 1/2* (norm(p1) + sqrt((p1(1) - fc(1))^2 + (p1(2) - fc(2))^2));
        end
    
        [~,a, e,per] = def_conic([0,0],fc,a);

        if(e>=0 && e<1) %considero solo le ellissi 
            
            
            %om = th1 - new_th1 + om1;
            %if om < 0 
                     %om = om + 360;
            %end

            om = rad2deg(atan2(per(2),per(1)))-180;
            if(om < 0)
                om = om + 360;
            end
            
            new_th1 = rad2deg(atan2(p1(2),p1(1)))-om;
            if(new_th1 < 0)
                new_th1 = new_th1 + 360;
            end

            new_th2 = rad2deg(atan2(p2(2),p2(1)))-om;
            if(new_th2 < 0)
                new_th2 = new_th2 + 360;
            end
            
            %p = a*(1-e^2);
            %new_th1 = rad2deg(acos(1/e*(p/norm(p1)-1)));
            %new_th2 = rad2deg(acos(1/e*(p/norm(p2)-1)));
            a = a*DIM;

            %controllo se schianto a terra
            if(a*(1-e) >= ALLERT)

                [new_r1,new_v1] = kep2car(a,e,iT,OMT,om,new_th1,mu);
                
                %if(norm(new_r1 - r1) > toll2)
                    %new_th1 = 360 - new_th1;
                    %om = th1 + om1 - new_th1;
                    %if om < 0
                        %om = om + 360;
                    %end
                    %[~,new_v1] = kep2car(a,e,iT,OMT,om,new_th1,mu);
                %end
                        
                [new_r2,new_v2] = kep2car(a,e,iT,OMT,om,new_th2,mu);

               % if(norm(new_r2 - r2) > toll2)
                  %new_th2 = 360 - new_th2;
                  %[~,new_v2] = kep2car(a,e,iT,OMT,om,new_th2,mu);
                %end
                COND = norm(new_r2 - r2) < toll2 && norm(new_r1 - r1) < toll2;
                if(COND)
                %calcoliamo le velocità rispetto all'orbita appena trovata
                    orbit_shape = [orbit_shape; a e om];
                    dv1 = norm(new_v1 - v1);
                    dv2 = norm(v2 - new_v2);
                    dv = [dv; dv1 + dv2];
                    th = [th; new_th1 new_th2];
                end
            end
        end
    end
    if(isempty(dv))
        kepT = [0, 0, 0, 0, 0, 0];
        dv = 0;
        th = [th1, 0, 0, th2];
    else
        [dv, k_vmin] = min(dv);
        kepT = [orbit_shape(k_vmin,1), orbit_shape(k_vmin,2),iT, OMT, orbit_shape(k_vmin,3), th(k_vmin,1)];
        th = [th1, th(k_vmin,1), th(k_vmin,2), th2];
    end
    end