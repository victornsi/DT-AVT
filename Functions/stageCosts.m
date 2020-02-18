function [r, dr, d2r] = stageCosts(digitalThread,U,n,derivativeFlg,varargin)
% About: Stage Cost r_t(D_t,u_t,y_t) wrapper
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% ///////////////////////////////////////////////////////////////////////// 
% Get Design and Operation
design = digitalThread.design;
operation = digitalThread.operation;

% Get Control
u = U.u;
ub = U.ub;

% Input override
if ~isempty(varargin)
    P = varargin{1}.manufacturingParam;
else
    P = design.inputDisturbances.samples.P(:,n);
end

% Only consider full cost for design only
switch ub
    case 'D'
        % As is
        % Calculate current stage cost of design
        [r,dr,d2r] = computeCost(design.region,U,P,derivativeFlg);

    case 'E'
        r = 0.005;
        dr = 0*u;
        d2r = zeros(size(u,1));
end

% Calculate operational costs for components in operation
co = 0;
for ii = 1:size(operation,2)
    co = co + 0.05*operation(ii).design.region.complexity.r;
end
r = r + co;

end

function [r,dr,d2r] = computeCost(region,U,Pinput,derivativeFlg)
% Create a copy of region
obj = region.copy();

% Parse input variables 
[a,z,s,p] = parseInput(U.u,obj);

% Update with input u
obj.designVariables.a = a;
obj.designVariables.z = z;
obj.designVariables.s = s;
obj.designVariables.p = p;

% Call complexity costs
G = obj.buildGlobalBasis(1);
[r,dr,d2r] = computeComplexity(obj,G,Pinput,derivativeFlg);
% num2cell([max(abs(a)),max(abs(z)),max(abs(s))])

end

function [a,z,s,p] = parseInput(u,r)
% Get Dimensions
dV = r.designVariables;
Na = size(dV.a,1);
Nz = size(dV.z,1);
Ns = size(dV.s(:),1);
Np = size(dV.p,1);

% Parse u
a = (u((1:Na)));
z = u(Na+(1:Nz));
s = reshape(u(Na+Nz+(1:Ns)),[],2);
p = u(Na + Nz + Ns + (1:Np));

end