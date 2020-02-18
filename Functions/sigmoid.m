function [s,ds,d2s] = sigmoid(p)
% Compute sigmoid
s = 1./(1+exp(-p));
ds = s.^2.*exp(-p);
d2s = (2*s.*ds - s.^2).*exp(-p);

end

% clc
% clear all
% x = [1 2 3.2]'/2;
% 
% dx = 1e-6;
% [k,dk,ddk] = sigmoid(x);
% for ii = 1:size(x,1)
%     xp = x;
%     xp(ii) = x(ii) + dx;
%     [kp,dkp] =  sigmoid(xp);
%     dk2(:,ii) = (kp-k)/dx;
%     ddk2(:,ii) = (dkp-dk)/dx;
%     
% end
% 
% [dk,dkp]
% num2cell([ddk,ddk2])