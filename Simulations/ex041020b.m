clear all;
close all;

% measurement
l=[0 500 1000]';
Z=[30.1 45 73.6]';
R=diag([0.001 0.001 1]);

for i=1:3
    Bi_p=[l(i) 0];
    Bi_q=Bi_p+[cosd(Z(i)) sind(Z(i))]*1500;
    plot([Bi_p(1) Bi_q(1)],[Bi_p(2) Bi_q(2)]);hold on;
end
axis equal;

%%
xlim([1150 1250]);

%%
X=[1000 700]';
dX=1;
eps=1e-6;

while norm(dX) > eps
    x=X(1);
    y=X(2);
    
    Z_bar=zeros(3,1);
    for i=1:3
        Z_bar(i)=atand(y/(x-l(i)));
    end
    dZ=Z-Z_bar;
    
    H=zeros(3,2);
    for i=1:3
        H(i,:)=[-y, x-l(i)]/((x-l(i))^2+y^2)*180/pi;    
    end

    P=inv(H'*inv(R)*H);
    K = P*H'*inv(R);
    
    dX=K*dZ
    X=X+dX;
    
    plot(X(1),X(2),'r*');
end

%%
plot_gaussian_ellipsoid(X,P,2);
axis equal;
xlim([1185 1225]);
