function [g,dg,margins] = stageConstraints(digitalThread,U,varargin)
% About: Stage Constraint g_t(D_t,u_t,y_t) wrapper
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% ///////////////////////////////////////////////////////////////////////// 
% Get Design
design = digitalThread.design;

% Get Control
u = U.u;
ub = U.ub;

% Input override
if ~isempty(varargin)
    samples.L = varargin{1}.loads;
    samples.A = varargin{1}.matAllow;
else
    samples = design.inputDisturbances.samples;
end

% Calculate current stage constraint of design
switch ub
	case 'D'
		[g,dg,margins] = computeConstraint(design.region,u,samples);
	case 'E'
		g = [-0.01;-0.01];
		dg = [u*0,u*0];
        margins = [];
end

end

function [g,dg,margins] = computeConstraint(region,u,samples)
% Create a copy of region
obj = region.copy();

% Get Input Disturbance samples
loads = samples.L;
matAllow = samples.A;

% Parse input variables 
[a,z,s,p] = parseInput(u,obj);

% Update with input u
obj.designVariables.a = a;
obj.designVariables.z = z;
obj.designVariables.s = s;
obj.designVariables.p = p;

% Run Margins
G = obj.buildGlobalBasis(0);
margins = computeMargins(obj,G,loads,matAllow,1);

% Compute Margins
g = (margins.f-1);
dg = [margins.ga;margins.gz];

% Compute Interior Function
ds = 1e-6;
f = 1;
[F] = obj.distanceFunction(obj.sensorLocations(s));
dFds1 = distanceFunction(obj,obj.sensorLocations(s + [ds 0]));
dFds2 = distanceFunction(obj,obj.sensorLocations(s + [0 ds]));
dFds = ([dFds1-F;dFds2-F])/(ds);
F = F*f;
dFds = dFds*f;

% Aggregate
[d,dds] = softmaximum(F,1);
dds = [dds(:);dds(:)].*dFds;

% Construct final constraint and gradient
g = [g;d];
dg = [[dg;0*dds;p*0],[dg*0;dds;p*0]];

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

function [fM,g1] = softmaximum(f,d)
% Maximization along d dimension fM,fMdf,fMdff
% Normalize
fn = sqrt(sum(f.^2,d));
fdn = bsxfun(@times,f,1./fn);

% Compute Soft Maximum
a = 500;%100*sqrt(size(f,d));
e = sum(exp(a*fdn),d);
lg = log(e)/a;
fM = bsxfun(@times,fn,lg);

% Compute Derivatives
dfM.dfn = lg;
dfM.dlg = fn;

dlg.de = 1./e./a;
dfn.dfn = 1;

de.dfdn = exp(a*fdn)*a;

dfdn.df = 1./fn;
dfdn.dfn = -bsxfun(@times,f,1./fn.^2);
dfn.df = fdn;
df.df = 1;

% Construct Derivatives
% Parents ('lg','fn') Children ('fn','e')
dfM.dfn =  dfM.dfn.*dfn.dfn;
dfM.de =  dfM.dlg.*dlg.de;

% Parents ('fn','e') Children ('fn','fdn')
dfM.dfdn =  bsxfun(@times,dfM.de,de.dfdn);
dfM.dfn =  dfM.dfn.*dfn.dfn;

% Parents ('fn','fdn') Childre ('f','fn')
dfM.df =  bsxfun(@times,dfM.dfdn,dfdn.df);
dfM.dfn = sum(dfM.dfdn.*dfdn.dfn,d) + dfM.dfn.*dfn.dfn;

% Parents ('fn','fdn') Children ('f')
dfM.df = dfM.df.*df.df + bsxfun(@times,dfM.dfn,dfn.df);

% Store
g1 = transpose(dfM.df);

end