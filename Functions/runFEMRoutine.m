function [region,f,ga,gz,sol] = runFEMRoutine(region,G,loads,dfdN,sol)
% About: Finite Element Routine
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% ///////////////////////////////////////////////////////////////////////// 
% Load Materials Data
load Materials

% Assemble As and bs
[bs,cs,ds,B] = assembleBs(region,loads,G);

% Check if dfdN is provided to collapse dimension
if isempty(dfdN)
    y = repmat(1:size(bs,1),size(bs,2),1)';
else
    y = dfdN;
end

% Solve System
[x,da,dz,Lf,sol] = solveSystem(region,bs,cs,ds,B,Materials,y,sol,G);

% Compute Gradient
f = sum(sum(y.*x));
ga = da;
gz = dz;

% Compute Field and Store
region.field.x = x;
region.field.f = f;
region.field.da = da;
region.field.dz = dz;
region.field.Lf = Lf;
region.field.A = sol.A;

end
function [bs,cs,ds,B] = assembleBs(region,loads,G)
% Assemble Bs
% Get Coordinates
iFill = find(strcmp({region.curves(:).type},'Fill'));
tbq = G.plate.boundary.tq;
ibq = G.plate.boundary.cq;
t = G.plate.boundary.ti;
in = G.plate.boundary.ci;

% Get Boundary Basis Function
M = G.plate.boundary.Nq;

% Get Quadrature Weights on boundary
vq = G.plate.boundary.wq;

% Determine Outer Boundary Inputs. Assume zero traction on all inner
% boundaries
L = loads;

% Transform normal and tangent loads to cartesian basis
[~,xt,xn] = region.computeTN(tbq(ibq==iFill),iFill);
sL = size(L);
r = [0 -1
     1  0];
xtr = xt*r;
xnr = xn*r;
L = reshape(L,size(xt,1),5,[]);
L = [bsxfun(@times,L(:,1,:),xn) + bsxfun(@times,L(:,2,:),xt),...
    (bsxfun(@times,L(:,3,:),xtr) - bsxfun(@times,L(:,4,:),xnr)),...
     L(:,5,:)];
L = reshape(L,sL);

o = sparse(size(xt,1),size(xt,1));
l = ones(size(xt,1),1);
Sp = @(x) spdiags(x,0,size(x,1),size(x,1));
T = [Sp(xn(:,1)) Sp(xt(:,1)) o         o      o
     Sp(xn(:,2)) Sp(xt(:,2)) o         o      o
     o       o                   Sp(xtr(:,1))  -Sp(xnr(:,1)) o
     o       o                   Sp(xtr(:,2))  -Sp(xnr(:,2)) o
     o       o                   o         o                  Sp(l)];

% Calculate Transformation to b
[~,dqs] = region.calculateBoundaryConditions(tbq(ibq==iFill));
S = repmat([1 1 -1 -1 1],size(L,1)/5,1).*(dqs(:,1:5)~=1); % Inner Product of coordinate basis unit vectors
S = S(:);
Vq = repmat(vq(ibq==iFill),5,1);%kron(eye(5),diag(vq(ibq==iFill)));
Mq = M(ibq==iFill,:);
B = reshape(Mq'*reshape(Vq.*S.*T,[],5*size(T,2)),[],size(T,2));
bs = reshape(Mq'*reshape(Vq.*S.*L,[],5*size(L,2)),[],size(L,2));

% Calculate Dirichlet Term
ds = zeros(size(t,1),10);
[~,xt,xn] = region.computeTN(t,iFill);
[~,ds(in==iFill,:)] = region.calculateBoundaryConditions(t(in==iFill));
ds1 = ds(:,1:5);
ds2 = ds(:,6:10);
ds = [ds1(:),ds2(:)];

% Compute Dirichlet Term
di = ds(:,1) > 0;
M = G.plate.boundary.Ni;
O = sparse(size(M,1),size(M,1));

% Compute x1 and x2 in n coordinates
x1 = [xn*[1;0],xt*[1;0]];
x2 = [xn*[0;1],xt*[0;1]];
p1 = [xn*[0;-1],xt*[0;-1]];
p2 = x1;
cs = [x1(:,1).*M x2(:,1).*M O O O
      x1(:,2).*M x2(:,2).*M O O O
      O O -p1(:,2).*M -p2(:,2).*M O
      O O p1(:,1).*M p2(:,1).*M O
      O O O O M];
cs = cs(di,:);
ds = repmat(ds(di,2),1,size(L,2));

end

function [x,da,dz,Lf,sol] = solveSystem(region,bs,cs,ds,B,Materials,dfdN,sol,G)
% Get System Matrix
if isempty(sol)
    y = repmat(1:size(bs,1),size(bs,2),1)';
    [As,~,~,ABDA] = assembleAs(region,Materials,reshape(y,[],5,size(bs,2)),reshape(y,[],5,size(bs,2)),[],G);
    
    % Solve System for x
    [x,Lf,A] = solveS(As,bs,cs,ds,B);
    
    % Store solution temporarily
    sol.x = x;   
    sol.Lf = Lf;
    sol.A = A;
    sol.ABDA = ABDA;
else
    x = sol.x;
    A = sol.A;
    Lf = sol.Lf;
    ABDA = sol.ABDA;
end

% Compute System Gradients
if ~isempty(dfdN)
    N = size(x,1);
    z = A\[dfdN;zeros(size(A,1)-N,size(dfdN,2))];
    y = -z(1:N,:);
    [~,da,dz] = assembleAs(region,Materials,reshape(y,[],5,size(bs,2)),reshape(x,[],5,size(bs,2)),ABDA,G);
else
    da = 0;
    dz = 0;
end
end

function [As,da,dz,ABDA] = assembleAs(region,Materials,y,x,ABDA,G)
% Get Basis Function and Basis Function Gradients on Quadrature Points
% Get Configuration Quadrature
wq = G.plate.planar.w;
JJinv = G.plate.planar.Jinv;

% Jinv
Jinv11 = reshape(JJinv(1,1,:),[],1);
Jinv12 = reshape(JJinv(1,2,:),[],1);
Jinv21 = reshape(JJinv(2,1,:),[],1);
Jinv22 = reshape(JJinv(2,2,:),[],1);

% Get Thicknesses and Ply Orientations
% Get Parameter Values
[a,da] = region.mappingA(G.plate.volume);
[z,dz] = region.mappingZ(G.thickness.planar);

% Resize z and dz for thickness
tq = region.mesh.plate.thicknessElem.q;
yq = region.mesh.plate.thicknessElem.w;
Nq = size(tq,1);
z = repmat(z,Nq,1);
dz.dZ = repmat(dz.dZ,Nq,1);

% Reshape a and z
a = reshape(a,size(wq,1),[]);
z = reshape(z,size(wq,1),[]);

% Compute Ply Transformations
ABDFlag = ~isempty(ABDA);
if ~ABDFlag
    [TT,UU,TTda,UUda] = plyTransformations(transpose(a(:)),size(a,1),size(a,2)); 
    ABD = laminateMatrices(z,tq,yq,Materials,TT,UU,TTda,UUda); 
    ABDA = ABD;
else
    ABD = ABDA;
end

% Compute Basis Functions
N = G.plate.planar.NT;
dN1 = G.plate.planar.dNT(:,:,1);
dN2 = G.plate.planar.dNT(:,:,2);
Nx1 = bsxfun(@times,Jinv11,dN1) + bsxfun(@times,Jinv21,dN2);
Nx2 = bsxfun(@times,Jinv12,dN1) + bsxfun(@times,Jinv22,dN2);
N = permute(N,[3 2 1]);
Nx1 = permute(Nx1,[3 2 1]);
Nx2 = permute(Nx2,[3 2 1]);
wq = reshape(wq,1,1,[]);

% Compute Strains
Nbar1 = [Nx1 0*Nx2];
Nbar2 = [0*Nx1 Nx2];
Nbar3 = [Nx2,Nx1];

Mbar1 = [-[N,N*0],Nx1];
Mbar2 = [-[N*0,N],Nx2];
Nbar = [Nbar1;Nbar2;Nbar3];
Mbar = [Mbar1;Mbar2];

% Compute ee,ek,kk,gg
Nbara = multiprod3(ABD.A,Nbar);
Nbarb = multiprod3(ABD.B,Nbar);
Nbard = multiprod3(ABD.D,Nbar);
Mbars = multiprod3(ABD.S,Mbar)*Materials.kappa;

ee = wq.*multiprod3(permute(Nbar,[2 1 3]),Nbara);
ek = wq.*multiprod3(permute(Nbar,[2 1 3]),Nbarb);
kk = wq.*multiprod3(permute(Nbar,[2 1 3]),Nbard);
gg = wq.*multiprod3(permute(Mbar,[2 1 3]),Mbars);

% Sum through quadrature
E = region.mesh.plate.elements;
Mq = size(region.mesh.plate.planarElem.w,1);
ee = sum(reshape(ee,size(ee,1),size(ee,2),[],Mq),4);
ek = sum(reshape(ek,size(ek,1),size(ek,2),[],Mq),4);
kk = sum(reshape(kk,size(kk,1),size(kk,2),[],Mq),4);
gg = sum(reshape(gg,size(gg,1),size(gg,2),[],Mq),4);

% Form Global Matrices
n = @(ii) E + max(E(:))*(ii-1);
r = @(ii) repmat(permute(n(ii),[2 3 1]),1,size(E,2),1);
c = @(ii) repmat(permute(n(ii),[3 2 1]),size(E,2),1,1);
R = @(ii,jj) reshape([r(ii) r(ii)
                      r(jj) r(jj)],[],1);
C = @(ii,jj) reshape([c(ii) c(jj)
                      c(ii) c(jj)],[],1);
Rg = @(ii,jj,kk) reshape([r(ii) r(ii) r(ii)
                          r(jj) r(jj) r(jj)
                          r(kk) r(kk) r(kk)],[],1);
Cg = @(ii,jj,kk) reshape([c(ii) c(jj) c(kk)
                          c(ii) c(jj) c(kk)
                          c(ii) c(jj) c(kk)],[],1);

% Sparse sums through repeated rows and columns
ee = sparse(R(1,2),C(1,2),ee(:));
ek = sparse(R(1,2),C(1,2),ek(:));
kk = sparse(R(1,2),C(1,2),kk(:));
gg = sparse(Rg(1,2,3),Cg(1,2,3),gg(:));

% Compute Derivatives
EE = repmat(E,Mq,1);
U1 = squeeze(x(:,1,:));
U2 = squeeze(x(:,2,:));
U3 = squeeze(x(:,3,:));
U4 = squeeze(x(:,4,:));
U5 = squeeze(x(:,5,:));

V1 = squeeze(y(:,1,:));
V2 = squeeze(y(:,2,:));
V3 = squeeze(y(:,3,:));
V4 = squeeze(y(:,4,:));
V5 = squeeze(y(:,5,:));

% Loop over quadrature points
for ii = size(a,1):-1:1
    % Get matrices
    Nbarq = Nbar(:,:,ii);
    Mbarq = Mbar(:,:,ii);
    ind = EE(ii,:);
    
    % Get Nodal Values
    u12 = [U1(ind,:);U2(ind,:)];
    u34 = [U3(ind,:);U4(ind,:)];
    u345 = [U3(ind,:);U4(ind,:);U5(ind,:)];
    
    v12 = [V1(ind,:);V2(ind,:)];
    v34 = [V3(ind,:);V4(ind,:)];
    v345 = [V3(ind,:);V4(ind,:);V5(ind,:)];
    
    % Compute (and sum over loads)
    eed(:,:,ii) = wq(ii)*(Nbarq*v12)*transpose(Nbarq*u12);   
    ekd(:,:,ii) = wq(ii)*(Nbarq*v12)*transpose(Nbarq*u34);   
    ked(:,:,ii) = wq(ii)*(Nbarq*v34)*transpose(Nbarq*u12);   
    kkd(:,:,ii) = wq(ii)*(Nbarq*v34)*transpose(Nbarq*u34);    
    ggd(:,:,ii) = wq(ii)*(Mbarq*v345)*transpose(Mbarq*u345);
end

% Compute Derivative wrt a
eeda = sum(ABD.Ada.*reshape(eed,[],1,size(eed,3)),1);
ekda = sum(ABD.Bda.*reshape(ekd,[],1,size(ekd,3)),1);
keda = sum(ABD.Bda.*reshape(ked,[],1,size(ked,3)),1);
kkda = sum(ABD.Dda.*reshape(kkd,[],1,size(kkd,3)),1);
ggda = sum(ABD.Sda.*reshape(ggd,[],1,size(ggd,3)),1)*Materials.kappa;

% Compute Derivative wrt z
eedz = sum(ABD.Adz.*reshape(eed,[],1,size(eed,3)),1);
ekdz = sum(ABD.Bdz.*reshape(ekd,[],1,size(ekd,3)),1);
kedz = sum(ABD.Bdz.*reshape(ked,[],1,size(ked,3)),1);
kkdz = sum(ABD.Ddz.*reshape(kkd,[],1,size(kkd,3)),1);
ggdz = sum(ABD.Sdz.*reshape(ggd,[],1,size(ggd,3)),1)*Materials.kappa;

% Assemble Global Stiffness Matrix
[As] = assembleStiffnessMatrix(ee,ek,kk,gg);

% Assemble Global Stiffness Matrices for derivatives
da = assembleStiffnessMatrix2(eeda,ekda,keda,kkda,ggda,da.dA);
dz = assembleStiffnessMatrix2(eedz,ekdz,kedz,kkdz,ggdz,dz.dZ);

end

function [As] = assembleStiffnessMatrix(ee,ek,kk,gg)
% Assemble Global Stiffness Matrix
% Order: u1,u2,phi1,phi2,w
Nu = size(ee,1)/2;
o = sparse(Nu,Nu);
As = [ee -ek [o;o]
     -transpose(ek) kk [o;o]
      repmat(o,1,5)];
Ng = size(gg,1);
As((end-Ng+1):end,(end-Ng+1):end) = As((end-Ng+1):end,(end-Ng+1):end) + gg;

end

function [d] = assembleStiffnessMatrix2(ee,ek,ke,kk,gg,k)
% Sum
ee = ee - ek - ke + kk + gg;
ee = transpose(squeeze(ee));

% Calculate gradient
d = transpose(transpose(ee(:))*k);

end


function [x,Lf,A] = solveS(As,bs,cs,ds,B)
% Assemble System
o = sparse(size(cs,1),size(cs,1));

% Constrained System
A = [As cs'
     cs o];
b = [bs
     ds];
 
% Solve system
N = size(As,1);
x = A\b;
x = x(1:N,:);

% Parse out inputs for filter
o = sparse(size(cs,1),size(B,2));
oo = sparse(size(bs,1),size(bs,2));
Lf.A = [B;o];%A\[B;o];
Lf.b = [oo;ds];%A\[oo;ds];
% Lf.A = Lf.A(1:N,:);
% Lf.b = Lf.b(1:N,:);

end
function [ABD] = laminateMatrices(z,tq,wq,Materials,TT,UU,TTda,UUda)
% Classical Laminate Plate Matrices
% Construct through thickness quadrature points
zq = 1/2*z.*tq';
dzq = 1/2*diag(tq);

% Construct through thickness quadrature weights
Iw = ones(size(tq,1),1)./size(tq,1);
zo = z*Iw;
wq1 = wq;
wq = 1/2*zo*wq1';
dwq = 1/2*wq1*Iw';

% Speecify angle gradient coefficients
I = eye(size(tq,1));

% Initialize
ABD.A=0;ABD.B=0;ABD.D=0;ABD.S=0;
ABD.Ada=0;ABD.Bda=0;ABD.Dda=0;ABD.Sda=0;
ABD.Adz=0;ABD.Bdz=0;ABD.Ddz=0;ABD.Sdz=0;

% Determine Q and G
MQ = repmat(Materials.Q,1,1,numel(z));
MG = repmat(Materials.G,1,1,numel(z));

% Precalculate
QQ = multiprod3(TT,MQ,permute(TT,[2 1 3 4]));
QQda = multiprod3(TT,MQ,permute(TTda,[2 1 3 4]));

GG = multiprod3(UU,MG,permute(UU,[2 1 3 4]));
GGda = multiprod3(UU,MG,permute(UUda,[2 1 3 4]));

% Loop Through Thickness
for ii = size(z,2):-1:1
    % Calculate Q and G
    Q = QQ(:,:,:,ii);%multiprod3(TT(:,:,:,ii),MQ,permute(TT(:,:,:,ii),[2 1 3]));%
    G = GG(:,:,:,ii);%multiprod3(UU(:,:,:,ii),MG,permute(UU(:,:,:,ii),[2 1 3]));%
    
    % Get Thickness and Quadrature
    pw = permute(wq(:,ii),[2 3 1]);
    pz = permute(zq(:,ii),[2 3 1]);
    
    % Calculate A, B, D, and S
    ABD.A = ABD.A + Q.*(pw);
    ABD.B = ABD.B + Q.*(pz.*pw);
    ABD.D = ABD.D + Q.*(pz.^2.*pw);
    ABD.S = ABD.S + G.*(pw);
    
    % Compute Gradient
    Qda = QQda(:,:,:,ii);%multiprod3(TT(:,:,:,ii),MQ,permute(TTda(:,:,:,ii),[2 1 3]));%
    Gda = GGda(:,:,:,ii);%multiprod3(UU(:,:,:,ii),MG,permute(UUda(:,:,:,ii),[2 1 3]));%
    Qda = Qda + permute(Qda,[2 1 3]);
    Gda = Gda + permute(Gda,[2 1 3]);
    
    % Gradient wrt angle
    pQda = reshape(Qda,9,1,[]);
    pGda = reshape(Gda,4,1,[]);    
    ABD.Ada = ABD.Ada + pQda.*(pw).*I(ii,:);
    ABD.Bda = ABD.Bda + pQda.*(pz.*pw).*I(ii,:);
    ABD.Dda = ABD.Dda + pQda.*(pz.^2.*pw).*I(ii,:);
    ABD.Sda = ABD.Sda + pGda.*(pw).*I(ii,:);
    
    % Gradient wrt thickness
    pQ = reshape(Q,9,1,[]);
    pG = reshape(G,4,1,[]);
    ABD.Adz = ABD.Adz + pQ.*dwq(ii,:);
    ABD.Bdz = ABD.Bdz + pQ.*pw.*dzq(ii,:) + pQ.*pz.*dwq(ii,:);
    ABD.Ddz = ABD.Ddz + 2*pQ.*pz.*pw.*dzq(ii,:) + pQ.*pz.^2.*dwq(ii,:);
    ABD.Sdz = ABD.Sdz + pG.*dwq(ii,:);
    
end

end

function [T,U,Tda,Uda] = plyTransformations(a,N,M)
% Ply transformations (in radians)
a = -a; % Need to invert to go from material to global axis

% Specify cos and sin
c = cos(a);
s = sin(a);

% Calculate T for inplane
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
T = reshape(T,3,3,N,M);
U = reshape(U,2,2,N,M);
Tda = -reshape(Tda,3,3,N,M);
Uda = -reshape(Uda,2,2,N,M);
end

