% Finds the previous vector v that the billiard ball was directed to reach
% current position on the boundary. Uses law of reflection to find prior v.
% INPUT: v = vector pointing in the direction of the path of the billiard
%        ball.
%        dz = tangent vector at the point on the boundary.
% OUTPUT: pv = prior vector v leading to current point on the boundary. 
function pv = priorv(v,dz)
    n = [-dz(2);dz(1)];
    n = n/norm(n);
    v2 = v-2*dot(n,v)*n;
    pv = v2/norm(v2);
