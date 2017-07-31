function ML = vertex_quad (nodl,nodm,sh_nod,p,s,psi1,psi2,areal,aream,p_c)
xm = p(1, nodm);
ym = p(2, nodm);
xl = p(1, nodl);
yl = p(2, nodl);
x = p_c(:,1);
y = p_c(:,2);
z = p_c(:,3);
local_l = find(nodl==sh_nod);
nsh_l = find(nodl~=sh_nod);
nsh_m = find(nodm~=sh_nod);
p_c = [xl(local_l), yl(local_l)];
Bl = [xl(nsh_l(1))-p_c(1) xl(nsh_l(2))-xl(nsh_l(1));
      yl(nsh_l(1))-p_c(2) yl(nsh_l(2))-yl(nsh_l(1))];
Bm = [xm(nsh_m(1))-p_c(1) xm(nsh_m(2))-xm(nsh_m(1));
      ym(nsh_m(1))-p_c(2) ym(nsh_m(2))-ym(nsh_m(1))];  
ML = ( 4*areal*aream/(4-2*s) ).*reshape(...
     psi1*( sum( ([ones(length(x),1)  x]*(Bl') - [y , y.*z]*(Bm') ).^2, 2 ).^(-1-s) ) +...
     psi2*( sum( ([ones(length(x),1)  x]*(Bm') - [y , y.*z]*(Bl') ).^2, 2 ).^(-1-s) )...
     , 5 , 5);
end

