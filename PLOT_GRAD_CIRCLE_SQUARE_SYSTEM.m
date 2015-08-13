function [] = PLOT_GRAD_CIRCLE_SQUARE_SYSTEM( p,e,t,c,u )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

size_u = size(u,1)/3;

u1 = zeros(size_u,size(u,2));
u2 = u1;
u3 = u1;

for i = 1:size_u
    u1(i,:) = u(i,:);
    u2(i,:) = u(i + size_u,:);
    u3(i,:) = u(i + 2*size_u,:);
end

figure(2)

for i = 1:size(u,2)
    subplot(2,2,1)
    [ux,uy] = pdecgrad(p,t,c,u1(:,i));
    ugrad = [ux;uy];
    pdeplot(p,e,t,'flowdata',ugrad)
    
    subplot(2,2,2)
    [ux,uy] = pdecgrad(p,t,c,u2(:,i));
    ugrad = [ux;uy];
    pdeplot(p,e,t,'flowdata',ugrad)
    
    subplot(2,2,3)
    [ux,uy] = pdecgrad(p,t,c,u3(:,i));
    ugrad = [ux;uy];
    pdeplot(p,e,t,'flowdata',ugrad)
    
    pause(0.01);

end

    
end

