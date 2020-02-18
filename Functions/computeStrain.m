function [S,df,Lf,Sinplane] = computeStrain(region,G,df,xi,F)
% About: Compute strains from displacement fields
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% /////////////////////////////////////////////////////////////////////////
% Get Basis Function and Basis Function Gradients
Field = region.field.x;
Lf = region.field.Lf;

% Jinv
JJinv = G.plate.planar.Jinv;
N = G.plate.planar.N;
dN1 = G.plate.planar.dN1;
dN2 = G.plate.planar.dN2;
Jinv11 = reshape(JJinv(1,1,:),[],1);
Jinv12 = reshape(JJinv(1,2,:),[],1);
Jinv21 = reshape(JJinv(2,1,:),[],1);
Jinv22 = reshape(JJinv(2,2,:),[],1);

% Compute Nx
Nx1 = bsxfun(@times,Jinv11,dN1) + bsxfun(@times,Jinv21,dN2);
Nx2 = bsxfun(@times,Jinv12,dN1) + bsxfun(@times,Jinv22,dN2);

% Compute Transformation
O = sparse(size(Nx1,1),size(Nx1,2));

if ~isempty(df) || F == 1
    Tinp.e = [Nx1 O O O O
              O Nx2 O O O
              Nx2 Nx1 O O O];
    Tinp.b = [O O -Nx1 O O
              O O O -Nx2 O
              O O -Nx2 -Nx1 O];
else
    Tinp = [];
end

% Compute Strains (voigt note the factor of 2 on shears)
z = region.mappingZ(G.thickness.planar);
[S,Sinplane,Lf] = computeStrainComponents(region,Tinp,N,Nx1,Nx2,Lf,z,xi,F,G.plate.planar.w);
S = reshape(S,[],8,size(Field,2));

% Compute Strain Derivatives
if ~isempty(df)
    T = [Tinp.e
         Tinp.b
         O O -N O Nx1
         O O O -N Nx2];
    df.dS = permute(reshape(df.dS,[],df.Ntq,8,size(df.dS,3)),[1 3 4 2]); % Sum through thickness
    temp = transpose(df.I)*reshape(sum(df.dS,4),size(df.dS,1),[]);
    df.dN = transpose(T)*reshape(temp,[],size(df.dS,3));
end

end

function [S,Sinplane,Lf] = computeStrainComponents(region,Tinp,N,Nx1,Nx2,Lf,z,xi,F,w)
% Compute Strains (voigt note the factor of 2 on shears)
% Parse input
Field = region.field.x;
I = reshape(1:(size(Field,1)),[],5);

U1 = Field(I(:,1),:);
U2 = Field(I(:,2),:);
P1 = Field(I(:,3),:);
P2 = Field(I(:,4),:);
W = Field(I(:,5),:);

e11 = Nx1*U1;
e12 = (Nx1*U2 + Nx2*U1);
e22 = Nx2*U2;
k11 = -Nx1*P1;
k12 = -(Nx1*P2 + Nx2*P1);
k22 = -Nx2*P2;
y1 = (Nx1*W - N*P1);
y2 = (Nx2*W - N*P2);

% Assemble strain components
S = [e11; e22; e12; k11; k22; k12 ;y1; y2];

% Check if filter information is required
if F == 0
    Lf = [];
    Sinplane = [];
    return
end

% Get Transformation Matrices and strains for sensors at sensor locations
Sp = @(x) spdiags(x,0,size(x,1),size(x,1));
Ti = region.sensorInterp(xi);
Ti = blkdiag(Ti,Ti,Ti);
Sinplane = [e11 + (z/2).*k11
            e22 + (z/2).*k22
            e12 + (z/2).*k12];
Tinplane = Tinp.e + Sp([z;z;z]/2)*Tinp.b;

% Interpolate at sensor locations
Sinplane = Ti*Sinplane;
Tinplane = Ti*Tinplane;

Lf.A = Tinplane*Lf.A;
Lf.b = Tinplane*Lf.b;

end