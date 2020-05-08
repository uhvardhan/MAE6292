function ex0317_QL
close all;

N=501;
t=linspace(0,1,N);
dt=t(2)-t(1);

% init guess
x_i=0*ones(N,1);
lambda_i=0*ones(N,1);

eps=1e-6;
delta=1;

while delta > eps
    
    % STM and particular sol
    Phi=zeros(2,2,N);
    p=zeros(2,N);
    Phi(:,:,1)=eye(2);
    p(:,1)=[0 0]';
    for k=1:N-1
        [A e]=eomLin(x_i(k),lambda_i(k));
        Phi_dot=A*Phi(:,:,k);
        p_dot=A*p(:,k)+e;
        Phi(:,:,k+1)=Phi(:,:,k)+Phi_dot*dt;
        p(:,k+1)=p(:,k)+p_dot*dt;
    end
    
    % find lambda_0
    x0=10;
    Phi_lx=Phi(2,1,N);
    Phi_ll=Phi(2,2,N);
    p_l=p(2,N);
    lambda0=inv(Phi_ll)*(-Phi_lx*x0-p_l);
    
    % compute x, lambda
    x=zeros(N,1);
    lambda=zeros(N,1);
    x(1)=x0;
    lambda(1)=lambda0;
    z0=[x0;lambda0];
    for k=2:N
        zk=Phi(:,:,k)*z0+p(:,k);
        x(k)=zk(1);
        lambda(k)=zk(2);
    end
    u=-lambda;
    
    subplot(2,2,1);
    plot(t,x);hold on;
    ylabel('x');
    subplot(2,2,2); 
    plot(t,lambda);hold on; 
    ylabel('\lambda');
    subplot(2,2,3);
    plot(t,u);hold on;
    ylabel('u');
    drawnow;
    pause;
    
    
    % update x^i and lambda^i
    delta=max(abs(x-x_i))+max(abs(lambda-lambda_i));
    
    x_i=x;
    lambda_i=lambda;
    
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


end


function [A e]=eomLin(x,lambda)

Hx=x-2*x*lambda;
Hl=-lambda-x^2;
Hxx=1-2*lambda;
Hll=-1;
Hlx=-2*x;
Hxl=-2*x;

A=[Hlx Hll;
    -Hxx -Hxl];
e=-A*[x;lambda]+[Hl; -Hx];
end


