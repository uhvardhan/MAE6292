function ex0317_Shooting
close all;

N=501;
t=linspace(0,1,N);

% init cond
x0=10;
lam0=0.4; % guess
P_x0=0;
P_l0=1;

eps=1e-6;
del_lambda0=1
alpha=0.1;

while norm(del_lambda0) > eps
    
    % forward integration
    X0=[x0;lam0;P_x0;P_l0];
    [t X]=ode45(@optcond,t,X0);
    x=X(:,1);
    lambda=X(:,2);
    P_x=X(:,3);
    P_l=X(:,4);
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
    subplot(2,2,4);
    plot(t,P_l);hold on;
    ylabel('P_l');
    
    drawnow;
   
    del_lambda0=-alpha*inv(P_l(N))*lambda(N);
    lam0=lam0+del_lambda0;
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
plot(t,P_l);hold on;
ylabel('P_l');

filename='ex0317_Shooting';
save(filename);
evalin('base',['load ' filename]);

end

function X_dot=optcond(t,X)
x=X(1);
lambda=X(2);
P_x=X(3);
P_l=X(4);

Hx=x-2*x*lambda;
Hl=-lambda-x^2;
Hxx=1-2*lambda;
Hll=-1;
Hlx=-2*x;
Hxl=-2*x;

x_dot=Hl;
lambda_dot=-Hx;
P_dot=[Hlx Hll; -Hxx -Hxl]*[P_x;P_l];

X_dot=[x_dot; lambda_dot; P_dot];


end

