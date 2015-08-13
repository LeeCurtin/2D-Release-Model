function [u1,p,e,t,c] = CIRCLE_SQUARE_SYSTEM(tspan,gd,sf,ns)
%CIRCLE_SQUARE_SYSTEM2 - Using a finite element method to evaluate the diffusion of
%matter out of/into a circular material (PLGA) into/out of a square (water).
tic
%CREATE MESH AND BOUNDARY
g = decsg(gd,sf,ns); %Creates g using the matrices produced by PDETOOLBOX
[p,e,t] = initmesh(g, 'hmax', 0.01);

%SET INITIAL CONDITIONS

[C,A,F]=assema(p,t,0,1,'0!1'); 
u01=A\F;
[C,A,F]=assema(p,t,0,1,'0!0'); 
u02=A\F;
[C,A,F]=assema(p,t,0,1,'0.05!0'); 
u03=A\F;

u0 = [u01;u02;u03];

% CREATE A PDE ENTITY FOR THREE PDEs WITH THREE DEPENDENT VARIABLES
numberOfPDE = 3;
pb = pde(numberOfPDE);
pg = pdeGeometryFromEdges(g); % Creates a geometry entity

% BOUNDARY CONDITIONS
Square = pdeBoundaryConditions(pg.Edges([1:4]),'q', [0,0,0;0,0,0;0,0,0], 'g', [0;0;0]); %Zero flux for Circle in Middle
% Square = pdeBoundaryConditions(pg.Edges([1:5]),'q', [0,0,0;0,0,0;0,0,0], 'g', [0;0;0]); %Zero flux for Circle on Bottom
% Circle = pdeBoundaryConditions(pg.Edges([5:8]),'q', [0,0,0;0,0,0;0,0,0], 'g', [0;0;0]); %Zero flux for Circle in circle
pb.BoundaryConditions = [Square];

%SET COEFFICIENTS OF PDEs

%DIFFUSION COEFFICIENTS
D_TMZ_water = 3600*50.e-5;
D_TMZ_paste = 3600*67.e-7;
D_MTIC_water = 3600*50.e-5;
D_MTIC_paste = 3600*67.e-7;
D_pH_water = 360*50.e-5;
D_pH_paste = 360*67.e-7;

c_TMZ = sprintf('%1.3f!%1.3f',D_TMZ_water,D_TMZ_paste);
c_MTIC = sprintf('%1.3f!%1.3f',D_MTIC_water,D_MTIC_paste); 
c_pH = sprintf('%1.3f!%1.3f',D_pH_water,D_pH_paste);

c = char(c_TMZ,c_TMZ,c_MTIC,c_MTIC,c_pH,c_pH);

%COEFFICIENTS OF A
a = ['0.00!0.00';'0.00!0.00';'0.00!0.00'];

%COEFFICIENTS OF F
% f = char('-0.5*u(1).*u(3)!-0.5*u(1).*u(3)','0.5*u(1).*u(3)!0.5*u(1).*u(3)','0.0!0.0');
f = char('-10*u(1,:).*u(3,:)!-10*u(1,:).*u(3,:)','10*u(1,:).*u(3,:)!10*u(1,:).*u(3,:)','0.0!0.0');

%COEFFICIENTS OF D
d = [1];

%SOLVE SYSTEM OF PDEs
u1 = parabolic(u0,tspan,pb,p,e,t,c,a,f,d);
toc

return


