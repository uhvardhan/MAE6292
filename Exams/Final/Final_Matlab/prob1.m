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
% EKF

X_bar=zeros(3,N);
P=zeros(3,3,N);
X_bar(:,1)=X0;
P(:,:,1)=P0;
Z_bar=zeros(2*Nref,N);

for l=1:N-1
    % Prediction
    prev_theta = X_bar(3,l);
    prev_x = X_bar(1,l);
    prev_y = X_bar(2,l);
    
    curr_x = prev_x + h*V*cos(prev_theta);
    curr_y = prev_y + h*V*sin(prev_theta);
    curr_theta = prev_theta + h*u(l);
    
    Ak = [1 0 -h*V*sin(prev_theta);
        0 1 h*V*cos(prev_theta);
        0 0 1];
    
    Bk = [0; 0; h];
    
    Gk = [h*cos(prev_theta) 0; h*sin(prev_theta) 0; 0 h];
    
    X_bar_kp = [curr_x curr_y curr_theta]';
    P_kp = Ak*P(:,:,l)*Ak' + Gk*Q*Gk';
    
    % Correction:
    for n=1:Nref
        refx = ref(1,n);
        refy = ref(2,n);
        
        delx = refx - curr_x;
        dely = refy - curr_y;
        
        dist = delx^2 + dely^2;
        
        % 2x1 matrix
        delP = [delx;dely];
        
        delP_body = [cos(curr_theta) sin(curr_theta);
            -sin(curr_theta) cos(curr_theta)]*delP;
        
        z_theta = atan2(delP_body(2), delP_body(1));
        
        % For each measurement, we have a 2x1 matrix
        Z_bar(2*n-1:2*n,l) = [sqrt(dist);z_theta] + v(2*n-1:2*n,l);
        
        % 2x3 matrix
        H = [-delx/sqrt(dist) -dely/sqrt(dist) 0;
            dely/dist -delx/dist -1];
        
        S = H*P_kp*H' + R(2*n-1:2*n,2*n-1:2*n);
        
        K = P_kp * H'/S;
        
        X_bar_kp = X_bar_kp + K*(z(2*n-1:2*n,l+1) - Z_bar(2*n-1:2*n,l));
        
        P_kp = (eye(3) - K*H)*P_kp;
    end
    
    P(:,:,l+1) = P_kp;
    
    X_bar(:,l+1) = X_bar_kp;
    
end
%%

figure(1);
plot(X_true(1,:),X_true(2,:),'r'); hold on;
plot(ref(1,:),ref(2,:),'b*');
plot(X_bar(1,:),X_bar(2,:),'b.')
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
axis equal;
legend('True', 'reference points', 'Estimated')
hold off;

figure(2);
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
