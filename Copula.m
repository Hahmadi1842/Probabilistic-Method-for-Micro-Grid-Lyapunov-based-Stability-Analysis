clc
clear all
close all


n = 1000;
rho = 0.4;
nu = 2;

rng('default') % for reproducibility
U = copularnd('t',[1 rho; rho 1],nu,n);
X = [gaminv(U(:,1),5,3) betainv(U(:,2),7,10)];

figure()
scatterhist(X(:,1),X(:,2),'Direction','out')


% 
% n = 1000;
% Rho = [1 .4 .2; .4 1 -.8];
% rng('default') % for reproducibility
% U = copularnd('Gaussian',Rho,n);
% X = [gaminv(U(:,1),2,1) betainv(U(:,2),2,2)];
% figure()
% scatterhist(X(:,1),X(:,2),'Direction','out')
% % subplot(1,1,1)
% plot3(X(:,1),X(:,2),X(:,3),'.')
% grid on
% view([-55, 15])
% xlabel('X1')
% ylabel('X2')
% zlabel('X3')