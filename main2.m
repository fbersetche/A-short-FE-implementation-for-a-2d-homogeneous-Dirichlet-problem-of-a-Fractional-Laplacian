clc
s = 0.5;
f = @(x,y) 1; 
cns = s*2^(-1+2*s)*gamma(1+s)/(pi*gamma(1-s));
load('data.mat');
nn = size(p,2); % number of nodes
nt = size(t,1) % number of elements
uh = zeros(nn,1); % discrete solution
K  = zeros(nn,nn); % stiffness matrix
b  = zeros(nn,1); % right hand side
% Compute areas
area = zeros(nt,1);
for i=1:nt
    aux = p( : , t(i,:) );
    area(i) = 0.5.*abs( det(  [ aux(:,1) - aux(:,3)  aux(:,2) - aux(:,3)]  ) );
end
% Build patches data structure
deg = zeros(nn,1);
for i=1:nt
    deg( t(i,:) ) = deg( t(i,:) ) + 1;
end
patches = cell(nn , 1);
for i=1:nn
    patches{i} = zeros( 1 , deg(i) );
end
for i=1:nt
    patches{ t(i,1) }(end - deg( t(i,1) ) + 1) = i;
    patches{ t(i,2) }(end - deg( t(i,2) ) + 1) = i;
    patches{ t(i,3) }(end - deg( t(i,3) ) + 1) = i;
    deg( t(i,:) ) = deg( t(i,:) ) - 1;
end
% Preallocate auxiliary memory
vl = zeros(6,2);
vm = zeros(6*nt,2);
norms = zeros(36,nt);
ML = zeros(6,6,nt);
empty = zeros(nt,1);   
aux_ind = reshape( repmat( 1:3:3*nt , 6 , 1 ) , [] , 1 ); 
empty_vtx = zeros(2,3*nt);
BBm = zeros(2,2*nt);
for l=1:nt-nt_aux % Main Loop
    edge =  [ patches{t(l,1)} patches{t(l,2)} patches{t(l,3)} ];
    [nonempty M N] = unique( edge , 'first' );
    edge(M) = [];
    vertex = setdiff(  nonempty , edge  );          
    ll = nt - l + 1 - sum( nonempty>=l );
    edge( edge<=l ) = []; % Elements Tm such that intersection(Tm,Tl) is an edge (and not have been computed yet) 
    vertex( vertex<=l ) = []; % Elements Tm such that intersection(Tm,Tl) is a vertex (and not have been computed yet)
    empty( 1:ll ) = setdiff_( l:nt , nonempty ); % Elements Tm such that intersection(Tm,Tl) is empty (and not have been computed yet)
    empty_vtx(: , 1:3*ll) = p( : , t( empty(1:ll) , : )' );
    nodl = t(l,:);
    xl = p(1 , nodl); yl = p(2 , nodl);
    Bl = [xl(2)-xl(1) yl(2)-yl(1); xl(3)-xl(2) yl(3)-yl(2)]'; % Bl*Tr + [xl(1) ; yl(1)] --> Tl
    b(nodl) = b(nodl) + fquad(area(l),xl,yl,f); % Volume forces
    K(nodl, nodl) = K(nodl, nodl) + triangle_quad(Bl,s,tpsi1,tpsi2,tpsi3,area(l),p_I) + comp_quad(Bl,xl(1),yl(1),s,cphi,R,area(l),p_I,w_I,p_T_12); % Compute the case Tm = Tl
    BBm(:,1:2*ll) = reshape( [ empty_vtx( : , 2:3:3*ll ) -  empty_vtx( : , 1:3:3*ll ) ,  empty_vtx( : , 3:3:3*ll ) - empty_vtx( : , 2:3:3*ll ) ] , [] , 2)' ; % BBm(:,2*m-1:2m)'*Tr + [xm(1) ; ym(1)] --> Tm
    vl = p_T_6*(Bl') + [ ones(6,1).*xl(1) ones(6,1).*yl(1) ]; % Quadrature points on Tl
    vm(1:6*ll,:) = reshape( permute( reshape( p_T_6*BBm(:,1:2*ll) , [6 1 2 ll] ) , [1 4 3 2] ) , [ 6*ll 2 ] ) + empty_vtx(: , aux_ind(1:6*ll) )'; % Quadrature points on Tm, m \in empty
    norms(:,1:ll) = reshape( distEucSq(vl,vm(1:6*ll,:)), 36 , [] ).^(-1-1*s) ; % Compute ||vl - vm||^-(2+2*s) for all m \in empty        
    ML(1:3,1:3,1:ll) =  reshape( phiA*norms(:,1:ll) , 3 , 3 , [] ); % Calculate and store the integrals on TlxTm, m \in empty
    ML(1:3,4:6,1:ll) =  reshape( phiB*norms(:,1:ll) , 3 , 3 , [] );
    ML(4:6,4:6,1:ll) =  reshape( phiD*norms(:,1:ll) , 3 , 3 , [] );
    ML(4:6,1:3,1:ll) =  permute( ML(1:3,4:6,1:ll) , [2 1 3] ) ; 
    % Assembling stiffness matrix
    for m=1:ll
        order = [nodl t( empty(m) , : )];
        K(order,order) = K(order,order) + ( 8*area(empty(m))*area(l) ).*ML(1:6,1:6,m);
    end
    for m=vertex  
        nodm = t(m,:);
        nod_com = intersect(nodl, nodm); 
        order = [nod_com nodl(nodl~=nod_com) nodm(nodm~=nod_com)];
        K(order,order) = K(order,order) + 2.*vertex_quad(nodl,nodm,nod_com,p,s,vpsi1,vpsi2,area(l),area(m),p_cube); 
    end
    for m=edge
        nodm = t(m,:);
        nod_diff = [setdiff(nodl, nodm) setdiff(nodm, nodl)]; 
        order = [ nodl( nodl~=nod_diff(1) ) nod_diff  ];
        K(order,order) = K(order,order) + 2.*edge_quad(nodl,nodm,nod_diff,p,s,epsi1,epsi2,epsi3,epsi4,epsi5,area(l),area(m),p_cube); 
    end
end 
uh(nf) = ( K(nf,nf)\b(nf) )./cns; % Solving the linear system
trimesh(t(1:nt-nt_aux , :), p(1,:),p(2,:),uh);