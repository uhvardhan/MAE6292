clear all;
close all;

x_bar=[12.3; 7.6];
P_x=diag([1.44 2.89]);

plot_gaussian_ellipsoid(x_bar,P_x,2);
axis equal;

%% Unscented Transform
kappa=30;
alpha=0.25;
beta=2;

[X_sigma W_m W_v]=UT_sigma(x_bar,P_x,kappa,alpha,beta);

n=2;
for i=1:2*n+1
    plot(X_sigma(1,i),X_sigma(2,i),'b+');
end

for i=1:2*n+1
    x=X_sigma(1,i);
    y=X_sigma(2,i);
    r=sqrt(x^2+y^2);
    theta=atan2(y,x);
    Y_sigma(:,i)=[r;theta];
end

[y_bar P_y]=UT_recover(Y_sigma,W_m,W_v);

figure;
plot_gaussian_ellipsoid(y_bar,P_y,2);
axis equal;
for i=1:2*n+1;
    plot(Y_sigma(1,i),Y_sigma(2,i),'b+');
end

plot(y_bar(1),y_bar(2),'b*');

%% Linearization

x=x_bar(1);
y=x_bar(2);
    r=sqrt(x^2+y^2);
    theta=atan2(y,x);

A=[x/r y/r;
    -y/r^2 x/r^2];
P_y_L=A*P_x*A';
y_bar_L=[r;theta];
plot_gaussian_ellipsoid(y_bar_L,P_y_L,2);
plot(y_bar_L(1),y_bar_L(2),'r*');


