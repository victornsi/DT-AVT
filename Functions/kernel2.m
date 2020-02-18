function [k] = kernel2(x,y,g,o)
% Modified from kernel to reduce memory overhead
% High Dimensional Version
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% ///////////////////////////////////////////////////////////////////////// 
% Form Arrays
X2 = sum(x.^2,2);
Y2 = sum(y.^2,2);
XY = x*transpose(y);

% Compute dX2
dX2 = X2 + transpose(Y2) - 2*XY + o;

% Compute kernel
k = exp(-g.*dX2);

end
