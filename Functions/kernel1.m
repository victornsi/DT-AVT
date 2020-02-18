function [k,dk,ddk] = kernel1(x,y,g,v,flg,z)
% About: Compute Gaussian RBF kernel and its derivatives
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% /////////////////////////////////////////////////////////////////////////
% Form Arrays
Nx = size(x,1);
Ny = size(y,1);
d = size(x,2);
X = reshape(x,Nx,1,d);
Y = reshape(y,1,Ny,d);

% Calculate displacements
dX = X-Y;
dX2 = sum(dX.^2,3);

% Pre multiply vector
if isempty(z)
    z = 1;
end

% Squared Exponential Kernel
k1 = exp(-g*dX2);
dk1 = (-2*g*k1).*dX;

% Exponential Kernel
if v ~= 1
    dX1 = sum(abs(dX),3);
    k2 = exp(-g*dX1);
    dk2 = -g*bsxfun(@times,k2,sign(dX));
else
    k2 = 0;
    dk2 = 0;
end

if length(z) > 1
    % Premultiply
    dk1 = transpose(transpose(z)*squeeze(dk1));
    dk2 = 0;%transpose(z)*squeeze(dk2);
%     dk1 = sum(dk1.*transpose(z),2);
%     dk2 = sum(dk2.*transpose(z),2);

end

% Second Derivatives   
if flg
    % Compute Second Derivative
    if length(z)==1
        for ii = d:-1:1
            for jj = d:-1:1
                % Squared Exponential Kernel
                ddk1(:,:,ii,jj) = (4*g^2.*(dX(:,:,ii).*dX(:,:,jj).*k1) - 2*g*(ii==jj)*(k1));
                
                % Only subgradients away from zero are defined for this. Otherwise,
                % add a large diagonal term
                ddk2(:,:,ii,jj) = (g^2*sign(dX(:,:,ii)).*sign(dX(:,:,jj)).*k2 - 1e6*g*k2.*(dX(:,:,ii)==0)*(ii==jj));
            end
        end
    elseif length(z)>1
        % Squared Exponential Kernel
        dXs = squeeze(dX);
        ddk1 = 4*g.^2.*(transpose(dXs)*((transpose(k1).*z).*dXs)) - 2*g*eye(d)*(k1*z);
        ddk2 = 0;
    end   
else
    ddk1 = 0;
    ddk2 = 0;
end

% Combine
k = v*k1 + (1-v)*k2;
dk = v*dk1 + (1-v)*dk2;
ddk = v*ddk1 + (1-v)*ddk2;

end

% clc
% x = [1 2]/2;
% y = [1 2 3 4]'/4;
% y = [y,y];
% v = 0;
% g = 3;
% flg = 1;
% 
% dx = 1e-6;
% [k,dk,ddk] = kernel(x,y,g,v,flg);
% for ii = 1:size(x,2)
%     xp = x;
%     xp(ii) = x(ii) + dx;
%     [kp,dkp] =  kernel(xp,y,g,v,flg);
%     dk2(:,:,ii) = (kp-k)/dx;
%     ddk2(:,:,:,ii) = (dkp-dk)/dx;
%     
% end
% num2cell([ddk-ddk2])