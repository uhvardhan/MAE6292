clear all;
close all;

% discretize
N=501;
t=linspace(0,1,N);
dt=t(2)-t(1);

% initial guess
u=0*ones(N,1);
alpha=0.5;
eps=1e-6;
H_u=1;

while norm(H_u) > eps
    x=zeros(N,1);
    lambda=zeros(N,1);
    
    % forward integration for x
    x(1)=10;
    for k=1:N-1
        f=-x(k)^2+u(k);
        x(k+1)=x(k)+f*dt;
    end
    
    % backward integ. for lambda
    lambda(N)=0;
    for k=N:-1:2
        H_x=x(k)-2*lambda(k)*x(k);
        lambda_dot=-H_x';
        lambda(k-1)=lambda(k)-lambda_dot*dt;
    end
    
    % update u
    H_u=zeros(N,1);
    for k=1:N
        H_u(k)=u(k)+lambda(k);
    end
    u_new=u-alpha*H_u;
    
    subplot(2,2,1);
    plot(t,x);hold on;
    ylabel('x');
    subplot(2,2,2);
    plot(t,lambda);hold on;
    ylabel('\lambda');
    subplot(2,2,3);
    plot(t,u);hold on;
    ylabel('u');
    subplot(2,2,4);
    plot(t,H_u);hold on;
    ylabel('H_u');
    
    drawnow;
    u=u_new;
    pause;
end


figure;
subplot(2,2,1);
plot(t,x);hold on;
ylabel('x');
subplot(2,2,2);
plot(t,lambda);hold on;
ylabel('\lambda');
subplot(2,2,3);
plot(t,u);hold on;
ylabel('u');
subplot(2,2,4);
plot(t,H_u);hold on;
ylabel('H_u');
