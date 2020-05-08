function ex0427
close all;

N=201;
t=linspace(0,20,N);
h=t(2)-t(1);
sigma_w=1e-2;
sigma_v=2e-1;
Qk=sigma_w^2;
Rk=sigma_v^2;

w=normrnd(0,sigma_w,1,N);
v=normrnd(0,sigma_v,1,N);

X_true=zeros(N,1);
X_true(1)=0;
uk=50;

% True Trajectory
for k=1:N-1
    X_true(k+1)=X_true(k)+h*uk+w(k);
end

for k=1:N
    z(k)=X_true(k)+v(k);
end

%plot(t,X_true,t,z);

NN=100;
XX=zeros(NN,N);
WW=zeros(NN,N);

X0=1;
sigma_x0=1;
XX(:,1)=normrnd(X0,sigma_x0,NN,1);
for i=1:NN
    WW(i,1)=1/NN;
end

figure(1);
plot(X_true(1),0,'r*');hold on;
plot_pdf_particle(XX(:,1),WW(:,1),'b-');
pause;

for k=1:4%N-1
    
    % prediction
    for i=1:NN
        XX(i,k+1)=normrnd(XX(i,k)+h*uk,sigma_w);
    end
    
    plot(X_true(k+1),0,'r*');hold on;
    plot_pdf_particle(XX(:,k+1),WW(:,k),'b-');
     pause; 
    
    % correction
    for i=1:NN
        WW(i,k+1)=WW(i,k)*normpdf(z(k+1),XX(i,k+1),sigma_v);
    end
    sumWW =sum(WW(:,k+1));
    WW(:,k+1)= WW(:,k+1)/sumWW;
    
    plot(z(k+1),0,'k*');hold on;
    plot_pdf_particle(XX(:,k+1),WW(:,k+1),'k-');
     pause 
    
end


X_bar=zeros(N,1);
for k=1:N
    for i=1:NN
        X_bar(k)=X_bar(k)+WW(i,k)*XX(i,k);
    end
    X_bar(k)=X_bar(k)/sum(WW(:,k));
end

figure;
plot(t,X_true,'r',t,X_bar,'b:');

figure;
plot_pdf_particle(XX(:,N),WW(:,N),'k-');

save ex0427;
evalin('base','clear all');
evalin('base','load ex0427');

end

function plot_pdf_particle(XX,WW,linetype)
NN=length(XX);
for i=1:NN
    plot([XX(i) XX(i)],[0 WW(i)],linetype);hold on;
end
xlim([-5 25 ]);
end

