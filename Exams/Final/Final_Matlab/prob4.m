clc;
clear all;
close all;

N=201;
t=linspace(0,20,N);
h=t(2)-t(1);
sigma_wV=1e-2;
sigma_wt=1e-2;
sigma_vr=1e-1;
sigma_vt=3*pi/180;

Nref=2;
ref=[0 10; 5 -5];
n=3;

Q=diag([sigma_wV^2 sigma_wt^2]);
R=diag([sigma_vr^2 sigma_vt^2 sigma_vr^2 sigma_vt^2]);

X0=[0.05 0.1 3*pi/180]';
P0=diag([0.2^2 0.2^2 (5*pi/180)^2]);

load prob1_wv;

X_true=zeros(3,N);
X_true(:,1)=[0 0 0]';

V=pi;

% True Trajectory
for k=1:N-1
    if k <= (N-1)/2
        u(k)=pi/5;
    else
        u(k)=-pi/5;
    end
    theta_k=X_true(3,k);
    X_true(:,k+1)=X_true(:,k)...
        +[h*(V+w_V(k))*cos(theta_k);h*(V+w_V(k))*sin(theta_k);h*(u(k)+w_t(k))];
end


% Measurement
for k=1:N
    x=X_true(1,k);
    y=X_true(2,k);
    theta=X_true(3,k);
    
    for j=1:Nref
        rx=ref(1,j);
        ry=ref(2,j);
        
        delx=rx-x;
        dely=ry-y;
        
        dis=sqrt(delx^2+dely^2);
        
        delP=[delx; dely];
        delP_body=[cos(theta) sin(theta);
            -sin(theta) cos(theta)]*delP;
        angle=atan2(delP_body(2),delP_body(1));
        
        z(2*j-1:2*j,k)=[dis; angle];
    end
    
    z(:,k)=z(:,k)+v(:,k);
end

%%
% Particle Filter

NN=1000;
XX_k=zeros(n,NN);
WW_k=zeros(NN,1);
XX_kp=zeros(n,NN);
WW_kp=zeros(NN,1);
X_bar=zeros(n,N);

for i=1:NN
    XX_k(:,i)=mvnrnd(X0,P0);
    WW_k(i)=1;
end
%
for i=1:NN
    X_bar(:,k)=X_bar(:,k)+WW_k(i)*XX_k(:,i);
end
X_bar(:,k)=X_bar(:,k)./sum(WW_k);
%
%%
for k=1:N-1
    
    % prediction
    for i=1:NN
        theta_k=XX_k(3,i);
        XX_kp_bar=XX_k(:,i) + [h*V*cos(theta_k);h*V*sin(theta_k);h*u(k)];
        Gk=[h*cos(theta_k) 0;h*sin(theta_k) 0;0 h];
        PP_kp=Gk*Q*Gk';
        XX_kp(:,i)=mvnrnd(XX_kp_bar,PP_kp);
    end
    
    % correction
    for i=1:NN
        curr_x = XX_kp(1,i);
        curr_y = XX_kp(2,i);
        curr_theta = XX_kp(3,i);
        
        for n=1:Nref
            refx = ref(1,n);
            refy = ref(2,n);
            
            delx = refx - curr_x;
            dely = refy - curr_y;
            
            dist = sqrt(delx^2 + dely^2);
            
            delP = [delx;dely];
            
            delP_body = [cos(curr_theta) sin(curr_theta);
                -sin(curr_theta) cos(curr_theta)]*delP;
            
            z_theta = atan2(delP_body(2), delP_body(1));
            
            Z_bar(2*n-1:2*n,i) = [dist z_theta];
        end
        
        WW_kp(i)=WW_k(i)*mvnpdf(z(:,k+1), Z_bar(:,i),R);
    end
    
    for i=1:NN
        X_bar(:,k+1)=X_bar(:,k+1)+WW_kp(i)*XX_kp(:,i);
    end
    X_bar(:,k+1)=X_bar(:,k+1)./sum(WW_kp);
    
    XX_k=XX_kp;
    WW_k=WW_kp;
    
    if rem(k,10)==0
        disp(k/N);
    end
    
end

figure(1);
plot(X_true(1,:),X_true(2,:),'r'); hold on;
plot(ref(1,:),ref(2,:),'b*');
plot(X_bar(1,:),X_bar(2,:),'b');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
axis equal;

for i=1:N
    x = X_bar(1,i);
    y = X_bar(2,i);
    theta = X_bar(3,i);
    
    for n=1:Nref
        refx = ref(1,n);
        refy = ref(2,n);
        
        delx = refx - x;
        dely = refy - y;
        
        dist = sqrt(delx^2 + dely^2);
        
        delP = [delx;dely];
        
        delP_body = [cos(theta) sin(theta);
            -sin(theta) cos(theta)]*delP;
        
        z_theta = atan2(delP_body(2), delP_body(1));
        
        Z_est(2*n-1:2*n,i) = [dist; z_theta];
    end
end
figure(2);
plot(t,z); hold on
plot(t,Z_est,'-');
xlabel('$t$','interpreter','latex');
ylabel('$z$','interpreter','latex');

figure(3);
plot(t,X_true,'r',t,X_bar,'b--');

