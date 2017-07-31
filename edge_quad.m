function ML = edge_quad (nodl,nodm,nod_diff,p,s,psi1,psi2,psi3,psi4,psi5,areal,aream,p_c)
xm = p(1, nodm);
ym = p(2, nodm);
xl = p(1, nodl);
yl = p(2, nodl);
x = p_c(:,1);
y = p_c(:,2);
z = p_c(:,3);
local_l = find(nodl~=nod_diff(1));
nsh_l = find(nodl==nod_diff(1));
nsh_m = find(nodm==nod_diff(2));
P1 = [xl(local_l(1)), yl(local_l(1))];
P2 = [xl(local_l(2)), yl(local_l(2))];
Bl = [P2(1)-P1(1) -P2(1)+xl(nsh_l);
      P2(2)-P1(2) -P2(2)+yl(nsh_l)];
Bm = [P2(1)-P1(1) -P2(1)+xm(nsh_m);
      P2(2)-P1(2) -P2(2)+ym(nsh_m)];
ML = ( 4*areal*aream/(4-2*s) ).*reshape(... 
      psi1*( sum( ([ones(length(x),1)  x.*z]*(Bl') - [1-x.*y  x.*(1-y)]*(Bm') ).^2, 2 ).^(-1-s) ) +...
      psi2*( sum( ([ones(length(x),1)  x]*(Bl') - [1-x.*y.*z  x.*y.*(1-z)]*(Bm') ).^2, 2 ).^(-1-s) ) +...
      psi3*( sum( ([(1-x.*y)  x.*(1-y)]*(Bl') - [ones(length(x),1) x.*y.*z]*(Bm') ).^2, 2 ).^(-1-s) ) +...
      psi4*( sum( ([1-x.*y.*z  x.*y.*(1-z)]*(Bl') - [ones(length(x),1) x]*(Bm') ).^2, 2 ).^(-1-s) ) +... 
      psi5*( sum( ([1-x.*y.*z  x.*(1-y.*z)]*(Bl') - [ones(length(x),1) x.*y]*(Bm') ).^2, 2 ).^(-1-s) )...
      , 4 , 4);
end