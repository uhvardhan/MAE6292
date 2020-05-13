function [y_bar P_y]=UT_recover(Y_sigma,W_m,W_v)
% Recover (y_bar,P_y) from the transformed sigma points

n=size(Y_sigma,1);
y_bar=zeros(n,1);
for i=1:length(W_m)
    y_bar=y_bar+W_m(i)*Y_sigma(:,i);
end
P_y=zeros(n,n);
for i=1:length(W_m)
    P_y=P_y+W_v(i)*(Y_sigma(:,i)-y_bar)*(Y_sigma(:,i)-y_bar)';
end

end