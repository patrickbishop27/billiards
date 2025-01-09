function r = v2r(theta,v,dz)


rguess =  dot(v,dz)/norm(dz);

findrspec = @(r) findr(r,v,dz);


rfin = Newtons(findrspec,@dfindr,rguess);


r = rfin; 

end







function y = findr(r,v,dz) 
    y = r^2 - dot(v,dz)^2/(norm(dz)^2);
end   
    
function dy = dfindr(r) 
	dy = 2*r; 
end