clear all;
close all;

%% time step and sim param
N=201;
t=linspace(0,20,N);
h=t(2)-t(1);

A=[1 h;0 1];
B=[h^2/2 ; h];
H=[1 0];

sigma_p=1e-3;
sigma_v=1e-3;
sigma_mu=5e-2;

Q=diag([sigma_p^2, sigma_v^2]);
R=sigma_mu^2;

%% True trajectory
x0=[0 0]';
P0=diag([0.1^2 0.1^2]);

w_p=normrnd(0,sigma_p,1,N);
w_v=normrnd(0,sigma_v,1,N);
w=[w_p;w_v];
v=normrnd(0,sigma_mu,1,N);

u=0.01;

X_true=zeros(2,N);
X_true(:,1)=x0+[0.05 -0.05]';
for k=1:N-1
    X_true(:,k+1)=A*X_true(:,k)+B*u+w(:,k);
end
for k=1:N
    z(k)=X_true(1,k)+v(k);
end

figure(1);
plot(t,X_true(1,:),'r',t,z,'k.');
xlabel('t');
ylabel('p');

%% Kalman Filter
X_bar=zeros(2,N);
X_bar(:,1)=x0;
P=zeros(2,2,N);
P(:,:,1)=P0;

for k=1:N-1
    % prediction
    X_bar_kp=A*X_bar(:,k)+B*u;
    P_kp=A*P(:,:,k)*A'+Q;
    
    % correction
    K=P_kp*H'*inv(H*P_kp*H'+R);
    X_bar_kp=X_bar_kp+K*(z(k+1)-H*X_bar_kp);
    P_kp=(eye(2)-K*H)*P_kp;
    
    
    X_bar(:,k+1)=X_bar_kp;
    P(:,:,k+1)=P_kp;
end


figure(1);hold on;
plot(t,X_bar(1,:),'b');

figure(2);
plot(t,X_true(2,:),'r',t,X_bar(2,:),'b',t(2:end),diff(z)/h,'k');
xlabel('t');
ylabel('$\dot x$','interpreter','latex');


for k=1:N
    sigma_x(k)=sqrt(P(1,1,k));
    sigma_x_dot(k)=sqrt(P(2,2,k));
end

figure;
plot(t,X_bar(1,:)-X_true(1,:),'b',t,3*sigma_x,'k',t,-3*sigma_x,'k');
xlabel('t');
ylabel('$e_x$','interpreter','latex');

figure;
plot(t,X_bar(2,:)-X_true(2,:),'b',t,3*sigma_x_dot,'k',t,-3*sigma_x_dot,'k');
xlabel('t');
ylabel('$e_{\dot x}$','interpreter','latex');

figure;
for k=1:10:N
    plot(X_true(1,k),X_true(2,k),'r.');hold on;
    plot(X_bar(1,k),X_bar(2,k),'b*');
    plot_gaussian_ellipsoid(X_bar(:,k),P(:,:,k),3);
end
axis equal;
xlabel('$x$','interpreter','latex');
ylabel('$\dot x$','interpreter','latex');
    




