function [X_sigma, W_m, W_v]=UT_sigma(x_bar, P,kappa,alpha,beta)
% Generate sigma points from (x_bar,P)

n=length(x_bar);
X_sigma=zeros(n,2*n+1);

X_sigma(:,1)=x_bar;
lambda=alpha^2*(n+kappa)-n;

sqrt_P=sqrtm((n+lambda)*P);
for i=1:n
    X_sigma(:,i+1)=x_bar+sqrt_P(:,i);
    X_sigma(:,i+n+1)=x_bar-sqrt_P(:,i);
end

W_m=zeros(2*n+1,1);
W_v=zeros(2*n+1,1);

W_m(1)=lambda/(n+lambda);
W_v(1)=W_m(1)+(1-alpha^2+beta);
for i=2:2*n+1
    W_m(i)=1/2/(n+lambda);
    W_v(i)=W_m(i);
end

end