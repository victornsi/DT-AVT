function [k,dk,d2k] = kernel(x,y,g,hessFlg,z,o)
% About: Kernel routine with derivatives (up to 2nd order). Current kernel
% is squared exponential
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
dX2 = sum(dX.^2,3) + o;

% Pre multiply vector
if isempty(z)
    z = 1;
end

% Compute Kernel
k = exp(-g*dX2);

% Compute Derivatives
dk = (-2*g*k).*dX;

% Premultiply
if length(z) > 1 && isnumeric(z)
    % Premultiply
    dk = transpose(transpose(z)*permute(dk,[2 3 1]));
end

% Second Derivatives   
if hessFlg
    if length(z)==1
        for ii = d:-1:1
            for jj = d:-1:1
                % Compute Terms
                d2k(:,:,ii,jj) = (4*g^2.*(dX(:,:,ii).*dX(:,:,jj).*k) - 2*g*(ii==jj)*(k));                                            
            end
        end
    elseif length(z)>1
        % Compute Second Derivative with Dimension Collapse
        dXs = squeeze(dX);
        d2k = 4*g.^2.*(transpose(dXs)*((transpose(k).*z).*dXs)) - 2*g*speye(d)*(k*z);
    end   
else
    d2k = 0;
end

end

% Debugging Code
% clc
% clear all
% x = [1 2]/2;
% y = [1 2 3 5]'/4;
% y = [y,y];
% g = 1;
% flg = 1;
% 
% dx = 1e-8;
% [k,dk,ddk] = kernel(x,y,g,flg,[],1.23);
% for ii = 1:size(x,2)
%     xp = x;
%     xp(ii) = x(ii) + dx;
%     [kp,dkp] =  kernel(xp,y,g,flg,[],1.23);
%     dk2(:,:,ii) = (kp-k)/dx;
%     ddk2(:,:,:,ii) = (dkp-dk)/dx;
%     
% end
% disp('Gradient')
% num2cell(dk-dk2)
% disp('Hessian')
% num2cell([ddk-ddk2])




% clc
% clear all
% x = [1 2 3]/3;
% y = [1 2 3 5]'/5;
% y = [y,y,y];
% g = 1;
% z = (1:4)';
% flg = 1;
% 
% dx = 1e-8;
% [k,dk,ddk] = kernel(x,y,g,flg,z,1.2);
% 
% for ii = 1:size(x,2)
%     xp = x;
%     xp(ii) = x(ii) + dx;
%     [kp,dkp] =  kernel(xp,y,g,flg,z,1.2);
%     dk2(ii) = (kp-k)/dx*z;
%     ddk2(:,ii) = (dkp-dk)/dx;
%     
% end
% disp('Gradient')
% num2cell(dk-dk2(:))
% disp('Hessian')
% num2cell([ddk-ddk2])