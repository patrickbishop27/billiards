% INPUT: t = parameter value corresponding to point on the billiard table. 
%            Value between 0 and 1. 
%        coefx = the coefficients for the cosine terms in the x coordinate 
%        for the table boundary.
%        coefy = coefficients for the sine terms in the y coordinate for
%        the table boundary.
%OUTPUT: z = [x,y] coordinates of the point on the billiard table
%        dz = [dx,dy] tangent line at the point on the billiard table.
function [z,dz] = Bill_Table(t,coefx,coefy)

    tpi = 2*pi; 
    z = zeros(2,length(t));
    dz = zeros(2,length(t));

    for kk =1: length(coefx)
        
        z(1,:) = z(1,:)+ coefx(kk)*cos(kk*tpi*t);
        z(2,:) = z(2,:) + coefy(kk)*sin(tpi*kk*t);
        dz(1,:) = dz(1,:) + -1*tpi*kk*coefx(kk)*sin(kk*tpi*t);
        dz(2,:) = dz(2,:)+ tpi*kk*coefy(kk)*cos(tpi*kk*t);

    end

    


end