function Dp = transitionModel(D,ub,u,inputs)
% About: Transition model for Digital Thread
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% ///////////////////////////////////////////////////////////////////////// 
% Initialize a new Digital Thread by inheriting previous
Dp = D.copy();
Dp.stage = Dp.stage + 1;

% Create a dummy component to access design
component = Component(Dp.design,inputs);

% Update design parameters with u
% Parse input variables
[a,z,s,p] = parseInput(u,Dp.design.region);

% Update with input u
Dp.design.region.designVariables.a = a;
Dp.design.region.designVariables.z = z;
Dp.design.region.designVariables.s = s;
Dp.design.region.designVariables.p = p;

% Update Basis
G = Dp.design.region.buildGlobalBasis(0);

% Update Dp with action u
switch ub
    case 'D'                
        % Get manufacturing times and FEM fields
        computeComplexity(Dp.design.region,G,inputs.manufacturingParam,'');
        computeMargins(Dp.design.region,G,inputs.loads,inputs.matAllow,0);
        
        % Collect manufacturing data from current design    
        d = component.measure('manufacture',G);    
        
        % Update design process with d
        Dp.design.inputDisturbances.update(d);
        
        % Manufacture and Deploy
        Dp.manufactureAndDeploy(inputs);
                
    case 'E'
        % Collect experimental data from experiments in current design
        d = component.measure('experiment',G);
        
        % Update design process with d
        Dp.design.inputDisturbances.update(d);

end

% Collect operational data from components in operation and update input disturbances
for ii = 1:size(Dp.operation,2)
    % Get Measurements
    d = Dp.operation(ii).measure('operate',G);
    
    % Update disturbance models with new measurements
    Dp.operation(ii).design.inputDisturbances.update(d);
    
    % Update design with latest
    if ii == size(Dp.operation,2)
        Dp.design.inputDisturbances.update(d);
    end
end

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