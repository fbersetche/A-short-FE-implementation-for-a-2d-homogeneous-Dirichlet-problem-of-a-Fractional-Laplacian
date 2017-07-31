h = 0.1;
[p,ee,tt]=initmesh('squareincircle','hmax',h);
% refine as many times as needed
% [p,ee,tt]=refinemesh('squareincircle',p,ee,tt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn=size(p,2);

% construction of t (triangles belonging to the square Omega) and dt (triangles 
% belonging to B\ Omega

IndexT=find(tt(4,:)==2); 
IndexB=find(tt(4,:)==1); 
t=tt(:,IndexT)';
dt=tt(:,IndexB)';
nt_aux=size(dt,1);
nf=[];

for j=1:nn
    if norm(p(:,j),inf) < 1
        nf=[nf; j];
    end
end

int_nod=unique(t(:,1:3));
bdrynodes=setdiff(int_nod,nf);
R=1.6; 
t=[t;dt];
t=t(:,1:3);
