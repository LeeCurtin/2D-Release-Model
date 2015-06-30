function [] = PLOT_CIRCLE_SQUARE_SYSTEM( u1,p,t,tspan )
%PLOT_CIRCLE_SQUARE_SYSTEM - plots the solution of the TMZ diffusion
%problem

%SPLIT SOLUTION OF SYSTEM INTO EACH PDE SOLUTION
u11 = zeros(size(p,2),length(tspan));
u12 = u11;
u13 = u12;

for j = 1:length(tspan)
    for i = 1:size(p,2)
        u11(i,j) = u1(i,j);
        u12(i,j) = u1(i+size(p,2),j);
        u13(i,j) = u1(i+2*size(p,2),j);
    end
end

%CREATE SOLUTION MOVIES
% fig = figure(2);
% u = fig.Units;
% fig.Units = 'normalized';
% % fig.Position = [0.3 0.3 0.7 0.7];
% fig.Color = [1 1 1];


F1 = pdeInterpolant(p,t,u11);
F2 = pdeInterpolant(p,t,u12);
F3 = pdeInterpolant(p,t,u13);

xgrid = -0.5:0.01:0.5;
ygrid = -0.5:0.01:0.5;
[X,Y] = meshgrid(xgrid,ygrid);

uout1 = evaluate(F1,X,Y);
uout2 = evaluate(F2,X,Y);
uout3 = evaluate(F3,X,Y);

% colormap(cool);
for i = 1:length(tspan)
    fig = figure(2);
    u = fig.Units;
    fig.Units = 'normalized';
    % fig.Position = [0.3 0.3 0.7 0.7];
    fig.Color = [1 1 1];
    colormap(cool);
%     subplot(2,2,1)
    Z1 = reshape(uout1(:,i),size(X));
    surf(xgrid,ygrid,Z1);
    axis([-0.5 0.5 -0.5 0.5 -0.01 1.01]);
    title('Temozolomide');
    zlabel('Concentration');
    shading interp;

    
    fig3 = figure(3);
    u = fig3.Units;
    fig3.Units = 'normalized';
    % fig3.Position = [0.3 0.3 0.7 0.7];
    fig3.Color = [1 1 1];
    colormap(cool);
%     subplot(2,2,2)
umax = max(max(u12));
umin = min(min(u12));
    Z2 = reshape(uout2(:,i),size(X));
    surf(xgrid,ygrid,Z2);
    axis([-0.5 0.5 -0.5 0.5 -0.01 0.3]);
    title('Active drug');
    zlabel('Concentration');
    shading interp;
    
    
    fig4 = figure(4);
    u = fig4.Units;
    fig4.Units = 'normalized';
    % fig4.Position = [0.3 0.3 0.7 0.7];
    fig4.Color = [1 1 1];
    colormap(cool);
%     subplot(2,2,3)
umax = max(max(u13));
umin = min(min(u13));
    Z3 = reshape(uout3(:,i),size(X));
    surf(xgrid,ygrid,Z3);
    axis([-0.5 0.5 -0.5 0.5 -0.001 0.051]);
    title('pH');
    zlabel('Concentration');
    shading interp;
%     pause(0.001)
    pause
end

end

