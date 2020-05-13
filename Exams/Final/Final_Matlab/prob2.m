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
% UKF
kappa=50;
alpha=0.25;
beta=2;
n=3;


X_bar=zeros(3,N);
P=zeros(3,3,N);
X_bar(:,1)=X0;
P(:,:,1)=P0;
Z_bar=zeros(2*Nref,N);
%%
for k=1:N-1
    % prediction
    % 1: choose sigma points
    % Output:
    % theta = X_bar(3,k);
    % X_k_sigma = 3x7 matrix
    [X_k_sigma, W_m, W_v]=UT_sigma(X_bar(:,k), P(:,:,k),kappa,alpha,beta);
    
    % 2: Pass sigma points through motion model and compute Gaussian
    % statistics
    X_kp_sigma = zeros(n,2*n+1);
    for i=1:2*n+1
        x = X_k_sigma(1,i);
        y = X_k_sigma(2,i);
        theta = X_k_sigma(3,i);
        
        curr_x = x + h*V*cos(theta);
        curr_y = y + h*V*sin(theta);
        curr_theta = theta + h*u(k);
        
        Ak = [1 0 -h*V*sin(theta);
            0 1 h*V*cos(theta);
            0 0 1];
        Bk = [0; 0; h];
        
        X_kp_sigma(:,i) = [curr_x; curr_y; curr_theta];
    end
    
    % 3: Recover mean and variance
    [X_bar_kp, P_kp] = UT_recover(X_kp_sigma, W_m, W_v);
    
    % 4: Add process noise
    %theta = X_bar_kp(3,1);
    Gk = [h*cos(theta) 0; h*sin(theta) 0; 0 h];
    P_kp = P_kp + Gk*Q*Gk';
    
    % Correction
    % 1: Get sigma points
    [X_sigma, W_m, W_v] = UT_sigma(X_bar_kp, P_kp, kappa, alpha, beta);
    % 2: Generate measurement for each sigma point
    Z_sigma = zeros(4,2*n+1);
    for i = 1:2*n+1
        curr_x = X_sigma(1,i);
        curr_y = X_sigma(2,i);
        curr_theta = X_sigma(3,i);
        
        for j = 1:Nref
            refx = ref(1,j);
            refy = ref(2,j);
            
            delX = refx - curr_x;
            delY = refy - curr_y;
            
            dist = sqrt(delX^2 + delY^2);
            delP = [delX;delY];
            
            delP_body = [cos(curr_theta) sin(curr_theta);
                -sin(curr_theta) cos(curr_theta)]*delP;
            
            z_theta = atan2(delP_body(2), delP_body(1));
            % Finally gives out a 4x7 matrix
            Z_sigma(2*j-1:2*j,i) = [dist; z_theta];
        end
    end
    % 3: recover mean and variance
    [Z_bar(:,k), P_z] = UT_recover(Z_sigma, W_m, W_v);
    P_z = P_z+R;
    
    P_xz = zeros(n,4);
    for i=1:2*n+1
        P_xz = P_xz + W_v(i) * (X_sigma(:,i) - X_bar_kp)*(Z_sigma(:,i)-Z_bar(:,k))';
    end
    
    % 4: Update
    K = P_xz * inv(P_z);
    X_bar_kp = X_bar_kp + K*(z(:,k+1) - Z_bar(:,k));
    P_kp = P_kp - K*P_z*K';
    
    P(:,:,k+1) = P_kp;
    X_bar(:,k+1) = X_bar_kp;
end
%%

figure;
plot(X_true(1,:),X_true(2,:),'r'); hold on;
plot(ref(1,:),ref(2,:),'b*');
plot(X_bar(1,:),X_bar(2,:),'b.')
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
axis equal;
legend('True', 'reference points', 'Estimated')
hold off;
%%
figure;
plot(t,z); hold on;
plot(t,Z_bar);
xlabel('$t$','interpreter','latex');
ylabel('$z$','interpreter','latex');
hold off

% Estimation Error
for m = 1:N
    sigma_x(m) = sqrt(P(1,1,m));
    sigma_y(m) = sqrt(P(2,2,m));
    sigma_theta(m) = sqrt(P(3,3,m));
end

diffX = X_bar(1,:) - X_true(1,:);
% 3-Sigma plots
figure
plot(t,diffX, 'b',t,3*sigma_x,'k',t,-3*sigma_x,'k');
xlabel('t');
ylabel('$e_x$','interpreter','latex');

diffY = X_bar(2,:) - X_true(2,:);
figure
plot(t,diffY, 'b',t,3*sigma_y,'k',t,-3*sigma_y,'k');
xlabel('t');
ylabel('$e_y$','interpreter','latex');

diffTheta = X_bar(3,:) - X_true(3,:);
figure
plot(t,diffTheta, 'b',t,3*sigma_theta,'k',t,-3*sigma_theta,'k');
xlabel('t');
ylabel('$e_\theta$','interpreter','latex');

figure
for a=1:10:N
    %plot(X_true(1,a), X_true(2,a),'r.');hold on;
    plot(X_bar(1,a), X_bar(2,a),'b*');
    plot_gaussian_ellipsoid(X_bar([1,2],a),P([1 2],[1,2],a),3);
end
axis equal;
xlabel('x');
ylabel('y');
