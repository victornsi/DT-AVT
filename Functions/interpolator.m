function [l,dl] = interpolator(x,g,X)
% About: Interpolator or basis functions for mesh interior and boundaries
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% /////////////////////////////////////////////////////////////////////////
% Calculate Monomials
d = size(x,2);
ind = fullfact(1+repmat(g,d,1));
ind = ind - 1;
ind(sum(ind,2)>g,:) = [];

l = 1;
% Monomials
for ii = d:-1:1
    l = l.*(x(:,ii)+1e-12).^(ind(:,ii)');
    dl(:,:,ii) = (ind(:,ii)').*(x(:,ii)+1e-12).^(-1);
end
dl = dl.*l;

% Gaussian RBF's
% [l,dl] = kernel(x,X,2,1,0,[]);
% 
% l = l./sum(l,2);
% dl = dl./sum(l,2) - l./sum(l,2).^2.*(sum(dl,2));
end

% clc
% L = 1;
% x = linspace(-1,1,200)'*L;
% 
% N = 21;
% ui = linspace(-1,1,N)'*L;
% 
% % Specify function
% e = 1e-9;
% f = @(y) sin(10*y);
% 
% % Specify Interpolator
% [l,dl] = interpolator(x,ui,1);
% [L] = interpolator(ui,ui,1);
% L = (L + e*eye(size(L)))\f(ui)
% g = l*L;
% % g = l*f(ui);
% 
% % Specify Kernel
% gg = 4;
% K = kernel(ui,ui,gg,1,0,[]);
% k = kernel(x,ui,gg,1,0,[]);
% K = (K+e*eye(size(K)))\f(ui)
% h = k*K;
% 
% figure(2)
% plot(x,f(x),'.-'), hold all
% plot(x,g,'-o')
% plot(x,h,'-^')
% plot(ui,0*ui,'*')
% hold off
% grid on
% %%
% clear all
% clc
% x = [1 2]/2;
% y = [1 2 3 4]'/4;
% y = [y,y];
% v = 0;
% g = 3.123;
% flg = 1;
% 
% dx = 1e-6;
% [k,dk] = interpolator(x,y,g);
% for ii = 1:size(x,2)
%     xp = x;
%     xp(ii) = x(ii) + dx;
%     [kp,dkp] =  interpolator(xp,y,g);
%     dk2(:,:,ii) = (kp-k)/dx;
%     
% end
% dk
% dk2

% clc
% X = 1;
% x = linspace(-1,1,200)'*X;
% 
% N = 20;
% gg = 4;
% % ui = linspace(-1,1,N)'*X;
% ui = 2*(rand(N,1)-0.5)
% 
% % Specify function
% e = 1e-9;
% f = @(y) sin(10*y);
% 
% % Specify Interpolator
% [l,dl] = interpolator(x,ui,4);
% [L,R] = kernel2(ui,ui,4,0);
% l = l/sqrt(N);
% L = L/sqrt(N);
% R = R/N;
% L = (L'*L + R + e*eye(size(L)))\(L'*f(ui));
% g = l*L;
% % Specify Kernel
% 
% K = kernel(ui,ui,gg,1,0,[]);
% k = kernel(x,ui,gg,1,0,[]);
% K = (K+e*eye(size(K)))\f(ui);
% h = k*K;
% 
% figure(2)
% plot(x,f(x),'.-'), hold all
% plot(x,g,'-ob')
% plot(x,h,'-^r')
% plot(ui,0*ui,'*')
% hold off
% grid on
% axis([-X X -2 2])