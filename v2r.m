% v2r converts vector v to the value r = cos(gamma).
% 
% INPUT: v = vector pointing in the direction of the path of the billiard
%        ball.
%        dz = tangent vector at the point on the boundary. 
%OUTPUT: r = cosine of the angle gamma where gamma is angle between the
%        tangent line at the point on the boundary and the vector v. 
%        theta1 = parameter value between 0 and 1 corresponding to the point 
%        on the boundary.
function r = v2r(v,dz)


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