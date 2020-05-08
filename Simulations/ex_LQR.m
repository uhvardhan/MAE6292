clear all;
%close all;

A=[0 1; 0 0];
B=[0 1]';
Q=[1 0; 0 0];
R=1;

Qf=[1 0; 0 1];

N=501;
t=linspace(0,10,N);
dt=t(2)-t(1);

q=100;
Q=q*[1 0; 0 0];

P=zeros(2,2,N);
P(:,:,N)=Qf;
for k=N:-1:2
    P_dot=-P(:,:,k)*A-A'*P(:,:,k)+P(:,:,k)*B*inv(R)*B'*P(:,:,k)-Q;
    P(:,:,k-1)=P(:,:,k)-P_dot*dt;
end

x=zeros(2,N);
K=zeros(1,2,N);
u=zeros(1,N);

x(:,1)=[1 1]';
for k=1:N-1
    K(:,:,k)=inv(R)*B'*P(:,:,k);
    u(k)=-K(:,:,k)*x(:,k);
    x_dot=A*x(:,k)+B*u(:,k);mail
    x(:,k+1)=x(:,k)+x_dot*dt;
end

for k=1:N
    K1(k)=K(1,1,k);
    K2(k)=K(1,2,k);
    p1(k)=P(1,1,k);
    p2(k)=P(1,2,k);
    p3(k)=P(2,2,k);
end

figure;
subplot(2,2,1);
plot(t,x);hold on;
ylabel('x');
subplot(2,2,2);
plot(t,u);hold on;
ylabel('u');
subplot(2,2,3);
plot(t,K1,t,K2);hold on;
ylabel('K');
subplot(2,2,4);
plot(t,p1,t,p2,t,p3);hold on;
ylabel('P');



%% infinite horizon LQR

K=lqr(A,B,Q,R);
x(:,1)=[1 1]';
for k=1:N-1
    u(k)=-K*x(:,k);
    x_dot=A*x(:,k)+B*u(:,k);
    x(:,k+1)=x(:,k)+x_dot*dt;
end

figure(1);
subplot(2,2,1);
plot(t,x,'r--');hold on;
ylabel('x');
subplot(2,2,2);
plot(t,u,'r--');hold on;
ylabel('u');
subplot(2,2,3);
plot([t(1) t(end)],[K(1) K(1)],'r--');hold on;
plot([t(1) t(end)],[K(2) K(2)],'r--');hold on;
ylabel('K');



