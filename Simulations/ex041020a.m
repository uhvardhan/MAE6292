clear all;
close all;

% prior belief
x_bar = 1;
M = 0.2^2; %with 95% prob,   abs (x - x_bar ) < 2-sigma : 2 meter
n = 1;

% sensor
z = 4;
R = 2^2; %with 95% prob,   error < 2-sigma : 1 meter

H = 1;

% estimation
S = H*M*H'+R;
K = M*H'*inv(S);
P = (eye(n) - K*H)*M*(eye(n)-K*H)' + K*R*K';

x_hat = x_bar + K*(z-H*x_bar);

% visualization
figure;
N = 1000;
x_vec = linspace(-2,8,N);
p_x = normpdf(x_vec, x_bar, M);
p_z = normpdf(x_vec, z, R);
p_x_z = normpdf(x_vec, x_hat, P);


plot(x_vec,p_x);
hold on;
plot(x_vec,p_z,'r');
xlabel('$x$','interpreter','latex');
ylabel('$p$','interpreter','latex');
plot(x_vec,p_x_z,'b','LineWidth',2);


