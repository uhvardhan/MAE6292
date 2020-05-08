clear all;
close all;

N=201;
t=linspace(0,20,N);
h=t(2)-t(1);
sigma_x=1e-3;
sigma_y=1e-3;
sigma_t1=3*pi/180;
sigma_t2=5*pi/180;

Q=diag([0 sigma_x^2 0 sigma_y^2]);
R=diag([sigma_t1^2 sigma_t2^2]);

X0=[1 0.5 3 0.8]';
P0=diag([0.1^2 0.1^2 0.1^2 0.1^2]);

w_x=normrnd(0,sigma_x,1,N);
w_y=normrnd(0,sigma_y,1,N);
w=[zeros(1,N); w_x; zeros(1,N); w_y];
v_t1=normrnd(0,sigma_t1,1,N);
v_t2=normrnd(0,sigma_t2,1,N);
v=[v_t1;v_t2];

X_true=zeros(4,N);
X_true(:,1)=X0+[0.1 0 -0.1 0]';

Ak=[1 h 0 0;
    0 1 0 0;
    0 0 1 h;
    0 0 0 1];
B=[0;0;h^2/2;h];
u=-0.1;
