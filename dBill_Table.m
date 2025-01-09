function dz = dBill_Table(t,coefx,coefy)

    tpi = 2*pi; 
    z = zeros(2,1);
    dz = zeros(2,1);

    for kk =1: length(coefx)
        
        z(1) = z(1)+ coefx(kk)*cos(kk*tpi*t);
        z(2) = z(2) + coefy(kk)*sin(tpi*kk*t);
        dz(1) = dz(1) + -1*tpi*kk*coefx(kk)*sin(kk*tpi*t);
        dz(2) = dz(2)+ tpi*kk*coefy(kk)*cos(tpi*kk*t);

    end

    


end