function pv = priorv(v,dz)
    n = [-dz(2);dz(1)];
    n = n/norm(n);
    v2 = v-2*dot(n,v)*n;
    pv = v2/norm(v2);
