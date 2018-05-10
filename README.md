
This is the code exposed in the article "A short implementation for a 2D Dirichlet 
problem of a fractional laplacian" available online here: https://arxiv.org/abs/1610.05558

In that paper we provide a comprehensive and simple 2D MATLAB finite element code for a 
Dirichlet problem of a fractional laplacian. The main code is written in about 80 lines 
and can be easily modified to deal with other kernels as well as with time dependent
problems.

----------------------------------------------------------------------------------------

For a quick test in a square domain, load "example_mesh.mat" in the MATLAB workspace and
then execute "main.m".

----------------------------------------------------------------------------------------

Alternatively, we provide a couple of basic mesh generators (mesh_generator_circle.m and
mesh_generator_square.m) providing meshes with a data structure compatible with the 
main code. 

In these programs a unitary ball and an unitary square are meshed togheter with a
 suitable auxiliary ball, using the geometry files "circleincircle.m" and 
"squareincircle.m"  respectively. The user may find it easy to modify them and create
customized examples.    

In order to use these meshes, just execute either "mesh_generator_circle.m" 
or "mesh_generator_square.m" and then "main.m".  

----------------------------------------------------------------------------------------

In case the user wants to provide a custom mesh, the main code needs the
following data loaded in the MATLAB workspace before the execution: 
  
`-p` , All the mesh vertex coordinates (including those lying in the auxiliary domain).

```
    | x1 | x2 | ... |xn |
p = |-------------------|
    | y1 | y2 | ... |yn |
```
     

`-t` , All the mesh triangles (those triangles that belong to the auxiliary domain must be
listed at the end).
```
t =   

| t11 | t12 | t13 |
|-----------------|
| t21 | t22 | t23 |
|-----------------|
   .     .     .
   .     .     .  
   .     .     .
|-----------------| 
| tk1 | tk2 | tk3 |
```


`-nt_aux` , The number of triangles in the auxiliary domain.

`-nf` , An index column vector containing the free nodes. 

`-R`  , The radius of the auxiliary domain.  

`-bdrynodes` , An index column vector containing the boundary nodes.  

----------------------------------------------------------------------------------------

WARNING!!!

If you get the error message "Undefined function or method 'pdist2' for input arguments 
of type 'double'" then you are using an old MATLAB version. In such a case use "main2.m" 
instead of "main.m".

----------------------------------------------------------------------------------------   

   
