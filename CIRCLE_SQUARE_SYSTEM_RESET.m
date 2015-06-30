function [u] = CIRCLE_SQUARE_SYSTEM3(gd,sf,ns)
%CIRCLE_SQUARE_SYSTEM3 - Using a finite element method to evaluate the diffusion of
%matter out of/into a circular material (PLGA) into/out of a square (water).

%CREATE MESH AND BOUNDARY
g = decsg(gd,sf,ns); %Creates g using the matrices produced by PDETOOLBOX
[p,e,t] = initmesh(g, 'hmax', 0.05);

%SET INITIAL CONDITIONS
u01 = zeros(size(p,2),1);
u02 = zeros(size(p,2),1); %MTIC is zero everywhere
u03 = ones(size(p,2),1); %pH is 1 everywhere
% ix = find(sqrt(p(1,:).^2 + (p(2,:)).^2)<0.3); %Domain: Circle in Middle 
ix = find(sqrt(p(1,:).^2 + (p(2,:)+0.2).^2)<0.3); %Domain: Circle on Bottom
u01(ix) = ones(size(ix)); %TMZ is one in the circle
u03(ix) = zeros(size(ix)); % Make pH zero in the circle
u0 = [u01;u02;u03];
% CREATE A PDE ENTITY FOR THREE PDEs WITH THREE DEPENDENT VARIABLES
numberOfPDE = 3;
pb = pde(numberOfPDE);
pg = pdeGeometryFromEdges(g); % Creates a geometry entity

% BOUNDARY CONDITIONS
% Square = pdeBoundaryConditions(pg.Edges([1:4]),'q', [0,0,0;0,0,0;0,0,0], 'g', [0;0;0]); %Zero flux for Circle in Middle
Square = pdeBoundaryConditions(pg.Edges([1:5]),'q', [0,0,0;0,0,0;0,0,0], 'g', [0;0;0]); %Zero flux for Circle on Bottom
pb.BoundaryConditions = [Square];

%SET COEFFICIENTS OF PDEs

%DIFFUSION COEFFICIENTS
D_TMZ_water = 0.05;
D_TMZ_paste = 0.005;
D_MTIC_water = 0.05;
D_MTIC_paste = 0.005;
D_pH_water = 0.01;
D_pH_paste = 0.001;

c_TMZ = sprintf('%1.3f!%1.3f',D_TMZ_water,D_TMZ_paste);
c_MTIC = sprintf('%1.3f!%1.3f',D_MTIC_water,D_MTIC_paste); 
c_pH = sprintf('%1.3f!%1.3f',D_pH_water,D_pH_paste);

c = char(c_TMZ,c_TMZ,c_MTIC,c_MTIC,c_pH,c_pH);

%COEFFICIENTS OF A
a = ['0.00!0.00';'0.00!0.00';'0.00!0.00'];

%COEFFICIENTS OF F
k = 0.5; %Rate constant k
f_k = sprintf('%1.3f',k);

f = char('-0.5*u(1,:).*u(3,:)!-0.5*u(1,:).*u(3,:)','0.5*u(1,:).*u(3,:)!0.5*u(1,:).*u(3,:)','0.0!0.0');
% f = char('-1.0*u(1).*u(3)!-1.0*u(1).*u(3)','1.0*u(1).*u(3)!1.0*u(1).*u(3)','0.0!0.0');

%COEFFICIENTS OF D
d = [1];

%SPLIT TSPAN INTO DIFFERENT TIMES TO ACCOUNT FOR TIME POINTS OF EXPERIMENT
tspan1 = [0:0.05:2.5];
tspan2 = tspan1;
tspan3 = tspan1;
tspan4 = tspan1;

%SOLVE SYSTEM OF PDEs FOR FIRST SECTION
u1 = parabolic(u0,tspan1,pb,p,e,t,c,a,f,d);

%SPLIT SOLUTION OF SYSTEM INTO EACH PDE SOLUTION
u11 = zeros(size(p,2),length(tspan1));
u12 = u11;
u13 = u12;

for j = 1:length(tspan1)
    for i = 1:size(p,2)
        u11(i,j) = u1(i,j);
        u12(i,j) = u1(i+size(p,2),j);
        u13(i,j) = u1(i+2*size(p,2),j);
    end
end

%RESET SQUARE TO DISTILLED WATER AT TIME TSPAN1
iy = find(sqrt(p(1,:).^2 + (p(2,:)+0.2).^2)>0.3);

u11_reset = u11(:,length(tspan1));
u11_reset(iy) = zeros(size(iy));

u12_reset = u12(:,length(tspan1));
u12_reset(iy) = zeros(size(iy));

u13_reset = u13(:,length(tspan1));
u13_reset(iy) = ones(size(iy));

u1_reset = [u11_reset;u12_reset;u13_reset];

%SOLVE SYSTEM OF PDES FOR SECOND SECTION
u2 = parabolic(u1_reset,tspan2,pb,p,e,t,c,a,f,d);

%SPLIT SOLUTION OF SYSTEM INTO EACH PDE SOLUTION

u21 = zeros(size(p,2),length(tspan2));
u22 = u21;
u23 = u22;

for j = 1:length(tspan2)
    for i = 1:size(p,2)
        u21(i,j) = u2(i,j);
        u22(i,j) = u2(i+size(p,2),j);
        u23(i,j) = u2(i+2*size(p,2),j);
    end
end

%RESET SQUARE TO DISTILLED WATER AT TIME TSPAN2
iy = find(sqrt(p(1,:).^2 + (p(2,:)+0.2).^2)>0.3);

u21_reset = u21(:,length(tspan2));
u21_reset(iy) = zeros(size(iy));

u22_reset = u22(:,length(tspan2));
u22_reset(iy) = zeros(size(iy));

u23_reset = u23(:,length(tspan2));
u23_reset(iy) = ones(size(iy));

u2_reset = [u21_reset;u22_reset;u23_reset];

%SOLVE SYSTEM OF PDES FOR THIRD SECTION
u3 = parabolic(u2_reset,tspan3,pb,p,e,t,c,a,f,d);

%SPLIT SOLUTION OF SYSTEM INTO EACH PDE SOLUTION

u31 = zeros(size(p,2),length(tspan3));
u32 = u31;
u33 = u32;

for j = 1:length(tspan3)
    for i = 1:size(p,2)
        u31(i,j) = u3(i,j);
        u32(i,j) = u3(i+size(p,2),j);
        u33(i,j) = u3(i+2*size(p,2),j);
    end
end

%RESET SQUARE TO DISTILLED WATER AT TIME TSPAN3
iy = find(sqrt(p(1,:).^2 + (p(2,:)+0.2).^2)>0.3);

u31_reset = u31(:,length(tspan3));
u31_reset(iy) = zeros(size(iy));

u32_reset = u32(:,length(tspan3));
u32_reset(iy) = zeros(size(iy));

u33_reset = u33(:,length(tspan3));
u33_reset(iy) = ones(size(iy));

u3_reset = [u31_reset;u32_reset;u33_reset];

%SOLVE SYSTEM OF PDES FOR THIRD SECTION
u4 = parabolic(u3_reset,tspan4,pb,p,e,t,c,a,f,d);

%SPLIT SOLUTION OF SYSTEM INTO EACH PDE SOLUTION
 
u41 = zeros(size(p,2),length(tspan4));
u42 = u41;
u43 = u42;

for j = 1:length(tspan4)
    for i = 1:size(p,2)
        u41(i,j) = u4(i,j);
        u42(i,j) = u4(i+size(p,2),j);
        u43(i,j) = u4(i+2*size(p,2),j);
    end
end

%CREATE SOLUTION MOVIES
fig = figure(3);
u = fig.Units;
fig.Units = 'normalized';
fig.Position = [0.3 0.3 0.7 0.7];

colormap(cool);
 x = linspace(-0.5,0.5,31);
 y = x;
[~, tn, a2, a3] = tri2grid(p,t,u0,x,y);
for i = 1:length(tspan1)
    subplot(2,2,1)
umax = max(max(u11));
umin = min(min(u11));
    u = tri2grid(p,t,u11(:,i),tn,a2,a3);
    j = find(isnan(u));
    u(j) = zeros(size(j));
    surf(x,y,u,'EdgeColor','none');
    axis([-0.5 0.5 -0.5 0.5 -0.001 1.01]);
    title('Temozolomide');
    shading interp;

    subplot(2,2,2)
umax = max(max(u12));
umin = min(min(u12));
    u = tri2grid(p,t,u12(:,i),tn,a2,a3);
    j = find(isnan(u));
    u(j) = zeros(size(j));
    surf(x,y,u,'EdgeColor','none');
    caxis([umin,umax]);
    axis([-0.5 0.5 -0.5 0.5 -0.001 0.25]);
    title('MTIC (or active drug)');
    shading interp;

    subplot(2,2,3)
umax = max(max(u13));
umin = min(min(u13));
    u = tri2grid(p,t,u13(:,i),tn,a2,a3);
    j = find(isnan(u));
    u(j) = zeros(size(j));
    surf(x,y,u,'EdgeColor','none');
    caxis([umin,umax]);
    axis([-0.5 0.5 -0.5 0.5 -0.001 1.01]);
    title('pH');
    shading interp;
    pause(0.001)
%     pause
end

[~, tn, a2, a3] = tri2grid(p,t,u1_reset,x,y);
for i = 1:length(tspan2)
    subplot(2,2,1)
umax = max(max(u21));
umin = min(min(u21));
    u = tri2grid(p,t,u21(:,i),tn,a2,a3);
    j = find(isnan(u));
    u(j) = zeros(size(j));
    surf(x,y,u,'EdgeColor','none');
    caxis([umin,umax]);
    axis([-0.5 0.5 -0.5 0.5 -0.001 1.01]);
    title('Temozolomide');
    shading interp;

    subplot(2,2,2)
umax = max(max(u22));
umin = min(min(u22));
    u = tri2grid(p,t,u22(:,i),tn,a2,a3);
    j = find(isnan(u));
    u(j) = zeros(size(j));
    surf(x,y,u,'EdgeColor','none');
    caxis([umin,umax]);
    axis([-0.5 0.5 -0.5 0.5 -0.001 0.25]);
    title('MTIC (or active drug)');
    shading interp;

    subplot(2,2,3)
umax = max(max(u23));
umin = min(min(u23));
    u = tri2grid(p,t,u23(:,i),tn,a2,a3);
    j = find(isnan(u));
    u(j) = zeros(size(j));
    surf(x,y,u,'EdgeColor','none');
    caxis([umin,umax]);
    axis([-0.5 0.5 -0.5 0.5 -0.001 1.01]);
    title('pH');
    shading interp;
    pause(0.001)
%     pause
end

[~, tn, a2, a3] = tri2grid(p,t,u2_reset,x,y);
for i = 1:length(tspan3)
    subplot(2,2,1)
umax = max(max(u31));
umin = min(min(u31));
    u = tri2grid(p,t,u31(:,i),tn,a2,a3);
    j = find(isnan(u));
    u(j) = zeros(size(j));
    surf(x,y,u,'EdgeColor','none');
    caxis([umin,umax]);
    axis([-0.5 0.5 -0.5 0.5 -0.001 1.01]);
    title('Temozolomide');
    shading interp;

    subplot(2,2,2)
umax = max(max(u32));
umin = min(min(u32));
    u = tri2grid(p,t,u32(:,i),tn,a2,a3);
    j = find(isnan(u));
    u(j) = zeros(size(j));
    surf(x,y,u,'EdgeColor','none');
    caxis([umin,umax]);
    axis([-0.5 0.5 -0.5 0.5 -0.001 0.25]);
    title('MTIC (or active drug)');
    shading interp;

    subplot(2,2,3)
umax = max(max(u33));
umin = min(min(u33));
    u = tri2grid(p,t,u33(:,i),tn,a2,a3);
    j = find(isnan(u));
    u(j) = zeros(size(j));
    surf(x,y,u,'EdgeColor','none');
    caxis([umin,umax]);
    axis([-0.5 0.5 -0.5 0.5 -0.001 1.01]);
    title('pH');
    shading interp;
    pause(0.001)
%     pause
end

[~, tn, a2, a3] = tri2grid(p,t,u3_reset,x,y);
for i = 1:length(tspan3)
    subplot(2,2,1)
umax = max(max(u41));
umin = min(min(u41));
    u = tri2grid(p,t,u41(:,i),tn,a2,a3);
    j = find(isnan(u));
    u(j) = zeros(size(j));
    surf(x,y,u,'EdgeColor','none');
    caxis([umin,umax]);
    axis([-0.5 0.5 -0.5 0.5 -0.001 1.01]);
    title('Temozolomide');
    shading interp;

    subplot(2,2,2)
umax = max(max(u42));
umin = min(min(u42));
    u = tri2grid(p,t,u42(:,i),tn,a2,a3);
    j = find(isnan(u));
    u(j) = zeros(size(j));
    surf(x,y,u,'EdgeColor','none');
    caxis([umin,umax]);
    axis([-0.5 0.5 -0.5 0.5 -0.001 0.25]);
    title('MTIC (or active drug)');
    shading interp;

    subplot(2,2,3)
umax = max(max(u43));
umin = min(min(u43));
    u = tri2grid(p,t,u43(:,i),tn,a2,a3);
    j = find(isnan(u));
    u(j) = zeros(size(j));
    surf(x,y,u,'EdgeColor','none');
    caxis([umin,umax]);
    axis([-0.5 0.5 -0.5 0.5 -0.001 1.01]);
    title('pH');
    shading interp;
    pause(0.001)
%     pause
end

return

