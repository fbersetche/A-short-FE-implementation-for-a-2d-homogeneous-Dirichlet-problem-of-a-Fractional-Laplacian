h = 0.07; % Size of the mesh
[p,ee,tt]=initmesh('circleincircle','hmax',h);
% refine as many times as needed
% [p,ee,tt]=refinemesh('circleminuscircle',p,ee,tt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn=size(p,2);

% construction of t (triangles belonging to the square Omega) and dt (triangles 
% belonging to B \ Omega.

IndexT=find(tt(4,:)==2); 
 
IndexB=find(tt(4,:)==1); 
t=tt(:,IndexT)';
dt=tt(:,IndexB)';
nt_aux=size(dt,1);
nf=[];
for j=1:nn
    if norm(p(:,j)) <0.999
        nf=[nf; j];
    end
end

int_nod=unique(t(:,1:3));
bdrynodes=setdiff(int_nod,nf);
R=1.1; 
t=[t;dt];
t=t(:,1:3);
