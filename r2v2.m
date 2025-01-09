function [v,flag] = r2v2(r,theta1,dz);



gamma = acos(r);
vguess = [-sin(2*pi*theta1+gamma);cos(2*pi*theta1+gamma)];

findvspec = @(v) findv(v,dz,r);
dfindvspec = @(v) dfindv(v,dz,r);

[vfin,stopv,~] = Newtons(findvspec,dfindvspec,vguess);
vfin = sign(dot(vfin,[-dz(2);dz(1)]))*vfin;
flag = 0; 
if(dot(dz,vfin)*r<0)
flag = 1;
	dz = dz(:)/norm(dz); 
	vfin = vfin - 2*dot(dz,vfin)*dz;
end

v = vfin;

end

function y = findv(v,dz,r)
    y = zeros(2,1);


    y(1) = ( dz(1)*v(1) + dz(2)*v(2))-sqrt(dz(1)^2+dz(2)^2)*r;
   	y(2) = v(1)^2+v(2)^2-1;
end

function   dy = dfindv(v,dz,r)
	dy = zeros(2,2); 

    dy(1,:) = [dz(1),dz(2)];
    dy(2,:) = 2*[v(1),v(2)]; 
end

