function [Margins,f,ga,gz] = computeMargins(region,Gb,loads,matAllow,FM)
% About: Margins routine
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% /////////////////////////////////////////////////////////////////////////
% Load Materials
load Materials

% Get Field
[region,~,~,~,sol] = runFEMRoutine(region,Gb,loads,[],[]);

if FM == 0
    % Only care about FEM 
    Margins = [];
    f = 0;
    ga = [];
    gz = [];
    return
end

% Margins are Evaluated on Quadrature Points of plate 
% Get Angle and Thickness Values
[a,da] = region.mappingA(Gb.plate.volume);
[z,dz] = region.mappingZ(Gb.thickness.planar);

% Replicate z values for thickness
Ntq = size(region.mesh.plate.thicknessElem.q,1);
z = repmat(z,Ntq,1);
dz.dZ = repmat(dz.dZ,Ntq,1);
t = Gb.plate.volume.t; % Thickness at quadrature in reference domain

% Get Midplane Strains
% Gradient Order: query points, components, design Variables
S = computeStrain(region,Gb,[],0,0);
S = repmat(S,Ntq,1);

% Derivative Info
dff.Ntq = Ntq;
dff.I = 1;

% % Margins are evaluated on both nodal and quadrature points
% % Get Angle and Thickness Values
% [a,da] = region.mappingA(Gb.margin.volume);
% [z,dz] = region.mappingZ(Gb.margin.planar);
% 
% % Replicate z values for thickness
% Ntq = size(region.mesh.plate.thicknessElem.q,1);
% z = repmat(z,Ntq,1);
% dz.dZ = repmat(dz.dZ,Ntq,1);
% t = Gb.margin.volume.t;
% 
% % Get Midplane Strains
% % Gradient Order: query points, components, design Variables
% S = computeStrain(region,Gb,[],0,0);
% 
% % Project onto nodes
% N = Gb.thickness.planar.N;
% w = Gb.thickness.planar.w;
% I = ((w.*N)'*N)\((w.*N)'); 
% Sn = reshape(I*reshape(S,size(S,1),[]),[],size(S,2),size(S,3));
% S = [S;Sn];
% 
% % Replicate Strain through thickness
% S = repmat(S,Ntq,1);
% 
% % Derivative Info
% dff.Ntq = Ntq;
% dff.I = [speye(size(w,1),size(w,1));I];
 
% Get Material Properties
Q = Materials.Q;
G = Materials.G;
Q3 = Materials.Q3;
Q3(3) = 0; % (e33) is zero for Mindlin-Reissner

% Reshape Material Allowables
matAllow = reshape(matAllow,size(matAllow,1),1,[]);

% Compute Ply Properties
% Need to invert and transpose T because in-plane strain transformation rule is different than for the stress
[TT,~,TTda,~] = plyTransformations(-a);
TT = permute(TT,[2 1 3]);
TTda = -permute(TTda,[2 1 3]);

% Transverse shear transformation is the same
[~,UU,~,UUda] = plyTransformations(a);

% Synthesize matrices
TTUU = [TT 0*TT 0*TT(:,[1 2],:)
        0*TT TT 0*TT(:,[1 2],:)
        0*TT([1 2],:,:) 0*TT([1 2],:,:) UU];
TTUUda = [TTda 0*TTda 0*TTda(:,[1 2],:)
          0*TTda TTda 0*TTda(:,[1 2],:)
          0*TTda([1 2],:,:) 0*TTda([1 2],:,:) UUda];

% Loop through Loadcases
NL = size(S,3);
fTL = 0;
for l = NL:-1:1
    % Compute Total Inplane Strain in Material Axis
    ee = sum(bsxfun(@times,permute(TT,[3 1 2]),permute(S(:,[1 2 3],l),[1 3 2])),3);
    kk = sum(bsxfun(@times,permute(TT,[3 1 2]),permute(S(:,[4 5 6],l),[1 3 2])),3);
    yy = sum(bsxfun(@times,permute(UU,[3 1 2]),permute(S(:,[7 8],l),[1 3 2])),3);
   
    % Compute Total Inplane Strain
    eT = ee + t.*bsxfun(@times,z/2,kk);

    % Compute Stresses and Strains
    % Order is [e11 e22 e12 e13 e23 e33] (6 components) for e = {e,s}
    [sQ,dsQ] = cSC3(Q,eT);
    [sG,dsG] = cSC2(G,yy);
    [sQ3,dsQ3] = cSC3(Q3,eT); % Last column is actually e33 = 0. Enforce instead with Q3(3) = 0 from above.
    s = [sQ,sG,sQ3];
    e = [eT,yy,0*yy(:,1,:)];
    
    % Compute Failure Criteria
    [f,df.ds] = calculateFailureCriteria(s,e,matAllow(:,1,l));
    
    % Compute Stresses and Strain Derivatives   
    ds.deT = [dsQ;zeros(2,3);dsQ3];
    ds.dyy = [zeros(3,2);dsG;zeros(1,2)];
    
    deT.dee = eye(3);
    deT.dz = t/2.*kk;
    deT.dkk = t/2.*z;
    dyy.dyy = 1;
    deekkyy.da = sum(bsxfun(@times,permute(TTUUda,[3 1 2]),permute(S(:,:,l),[1 3 2])),3);
    deekkyy.dS = TTUU;
    dz.dz = 1;
    
    % Construct Derivatives
    df.deT = df.ds*ds.deT;
    df.dyy = df.ds*ds.dyy;
    
    df.dee = df.deT*deT.dee;
    df.dkk = bsxfun(@times,df.deT,deT.dkk);
    df.dz  = sum(df.deT.*deT.dz,2);
    df.dyy = df.dyy*dyy.dyy;
    df.deekkyy = [df.dee,df.dkk,df.dyy];
    
    df.da = sum(df.deekkyy.*deekkyy.da,2);
    df.dS = sum(bsxfun(@times,permute(deekkyy.dS,[3 2 1]),permute(df.deekkyy,[1 3 2])),3);
    df.dz = df.dz.*dz.dz;
       
    % Compute Soft Max through Quadrature
    [fN,dfN.df] = softmaximum(f,1);
    
    % Store for Loadcases
    dff.da(:,l) = dfN.df(:).*df.da;
    dff.dz(:,l) = dfN.df(:).*df.dz;
    dff.dS(:,:,l) = dfN.df(:).*df.dS;
    fNL(l) = fN;
    fTL = fTL + 1/NL*f;
     
end

% Compute Gradient through strain
[~,dff] = computeStrain(region,Gb,dff,0,0);
[region,~,~,~] = runFEMRoutine(region,Gb,loads,dff.dN,sol);

% Compute Final Gradients - Not divided by NL in runFEMroutine
dN.dA = region.field.da;
dN.dZ = region.field.dz;

df.dA = (sum(transpose(da.dA)*dff.da,2) + dN.dA)/NL;
df.dZ = (sum(transpose(dz.dZ)*dff.dz,2) + dN.dZ)/NL;
f = sum(fNL)/NL;

% Store
[~,ind] = max(fNL);
Margins.fTL = fTL;
Margins.f = f;
Margins.ga = df.dA;
Margins.gz = df.dZ;
Margins.loads = loads(:,ind);
Margins.matAllow = matAllow(:,ind);
Margins.ind = ind;
ga = df.dA;
gz = df.dZ;

% Store into region
region.margins = Margins;

end

function [T,U,Tda,Uda] = plyTransformations(a)
% Ply transformations (in radians)
a = transpose(a);

% Specify cos and sin
c = cos(a);
s = sin(a);

% Calculate T for in plane
T = [c.^2
     s.^2
    -s.*c
     s.^2
     c.^2
     s.*c
     2*s.*c
    -2*s.*c
     c.^2-s.^2];
Tda = [-2.*c.*s
        2.*s.*c
       -c.*c + s.*s
       2.*s.*c
       -2.*c.*s
        c.*c - s.*s
       2*(c.*c - s.*s)
      -2*(c.*c - s.*s)
      -4.*c.*s];

% Calculate U for shear
U = [c
    -s
     s
     c];
Uda = [-s
       -c
        c
       -s];

% Reshape
T = reshape(T,3,3,[]);
U = reshape(U,2,2,[]);
Tda = reshape(Tda,3,3,[]);
Uda = reshape(Uda,2,2,[]);
end

function [fM,g1] = softmaximum(f,d)
% Maximization along d dimension fM,fMdf
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

function [s,sd] = cSC3(Q,strain)
N = size(strain);
N(2) = 1;
s = zeros(N);
s(:,1,:,:) = Q(1,1)*strain(:,1,:,:) + Q(1,2)*strain(:,2,:,:) + Q(1,3)*strain(:,3,:,:);
if size(Q,1) > 1
    s(:,2,:,:) = Q(2,1)*strain(:,1,:,:) + Q(2,2)*strain(:,2,:,:) + Q(2,3)*strain(:,3,:,:);
    s(:,3,:,:) = Q(3,1)*strain(:,1,:,:) + Q(3,2)*strain(:,2,:,:) + Q(3,3)*strain(:,3,:,:);
end

% Compute Adjoint
sd = Q;

end
function [s,sd] = cSC2(Q,strain)
N = size(strain);
N(2) = 1;
s = zeros(N);
s(:,1,:,:) = Q(1,1)*strain(:,1,:,:) + Q(1,2)*strain(:,2,:,:);
s(:,2,:,:) = Q(2,1)*strain(:,1,:,:) + Q(2,2)*strain(:,2,:,:);

% Compute Adjoint
sd = Q;
end

% % Get Initial
% [~,F,dA,dZ,sol] = runFEMRoutine(region,Gb,loads,[],[]);
% 
% a = region.designVariables.z;
% da = 1e-5;
% for ii = 1:5
%     [ii]
%     % Perturb
%     ap = a;
%     ap(ii) = ap(ii) + da;
%     region.designVariables.z = ap;
%     
%     % Evaluate
%     [~,Fp] = runFEMRoutine(region,Gb,loads,[],[]);
%     dAn(ii,1) = (Fp-F)/(da);
%     
% end
% 
% dZ(1:5)
% dAn
% 
% pause

