function Materials = createMaterials()
% About: Generate Material Properties for lamina
% Each lamina has the same properties for now and is treated as tranversely
% isotropic (1,2=3)
% E22 = E33
% v12 = v13
% v21 = v31
% v23 = v32
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% /////////////////////////////////////////////////////////////////////////

% Specify Material
% All material properties in Pa where applicable
% Material Data Taken from Table 3.2 Page 36 of Spec Sheet NCP-RP-2008-04 N/C
% RTD - Mean - Normalized Data (where applicable)
ii = 1;
si2pa = 6894.76; % psi to Pa
M(ii).kappa = 5/6;
M(ii).Name = 'MTM45-1/AS4';
M(ii).E11c = 17.02e6*si2pa;
M(ii).E11t = 18.51e6*si2pa;
M(ii).E22c = 1.25e6*si2pa;
M(ii).E22t = 1.15e6*si2pa;
M(ii).v12c = 0.31;
M(ii).v21c = 0.02;
M(ii).G12 = 0.53e6*si2pa;

% Specify Fiber Material Properties
M(ii).E11f = 231e9*si2pa;
M(ii).v12f = 0.2;

% Specify Material Stress Strengths
M(ii).F11c = 202.80e3*si2pa;
M(ii).F11t = 270.76e3*si2pa;
M(ii).F22c = 26.81e3*si2pa;
M(ii).F22t = 6.92e3*si2pa;
M(ii).F12 = 9.36e3*si2pa;% Use 5% strain value 6.67

% Specify Material Strain Strengths
M(ii).S11c = 11107e-6;
M(ii).S11t = 14239e-6;
M(ii).S22c = 10899e-6;
M(ii).S22t = 14136e-6;
M(ii).S12 = 18679e-6;% Use 5% strain value 6.67

% Compute Stiffness and Shear Matrices in Material Coordinates (average
% compression and tension values for now)
E11 = 0.5*(M(ii).E11c + M(ii).E11t);
E22 = 0.5*(M(ii).E22c + M(ii).E22t);
E33 = E22;
v12 = M(ii).v12c;
v21 = M(ii).v21c;
v13 = v12;
v31 = v21;
v23 = 0;
v32 = v23;
G12 = M(ii).G12;
G13 = G12;
G23 = E22/(2*(1+v23));

% Isotropic Case
% E11 = 0.5*(M(ii).E11c + M(ii).E11t);
% E22 = E11;
% E33 = E22;
% v12 = M(ii).v12c;
% v21 = M(ii).v12c;
% v13 = v12;
% v31 = v21;
% v23 = v12;
% v32 = v23;
% G12 = E11/(2*(1+v12));
% G13 = G12;
% G23 = E22/(2*(1+v23));

% Compute Compliance Matrix
S = [1/E11 -v21/E22 -v31/E33 0 0 0
     -v12/E11 1/E22 -v32/E33 0 0 0
     -v13/E11 -v23/E22 1/E33 0 0 0
     0 0 0 1/G23 0 0
     0 0 0 0 1/G13 0
     0 0 0 0 0 1/G12];
 
% Symmetrize since Modulus and Poisson constraints aren't exact
S = 0.5*(S+S');

% Compute Stiffness Matrix
C = S\eye(6);

% Extract Q
Q = blkdiag(C(1:2,1:2),C(6,6));

% Extract Q3 (33 direction)
Q3 = C(3,1:3);%[C(3,1:2),0];

% Extract G
G = diag([C(5,5),C(4,4)]);

% Store Data
M(ii).Q = Q;
M(ii).G = G;
M(ii).Q3 = Q3;
Materials = M;

save('Materials\Materials.mat','Materials')
end