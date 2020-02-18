function [f,g] = calculateFailureCriteria(s,e,matAllow)
% About: Calculate failure criteria
% (intra-laminar/phenomenological/fiber/matrix) using a variety of popular
% methods: 'Tsai-Wu','Max-Strain','Max-Stress','Hashin','Christensen','Puck','LaRC03'
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% /////////////////////////////////////////////////////////////////////////
% Calculate Streses
% Specify Labels
names = {'Tsai-Wu','Max-Strain','Max-Stress','Hashin','Christensen','Puck','LaRC03'};

% Load Materials Data
% load Materials
Materials = [];

% Convert to Voigt
[sv] = voigt(s);
ev = voigt(e);

% Calculate Margins
[f,g] = TsaiWu(sv,Materials,matAllow);

% Reverse Transformation
[g] = unvoigt(g);
% [f(2,:),l{2}] = maxStrain(ev,Materials);
% [f(3,:),l{3}] = maxStress(sv,Materials);
% [f(4,:),l{4}] = Hashin(sv,Materials);
% [f(5,:),l{5}] = Christensen(sv,Materials);
% [f(6,:),l{6}] = Puck(sv,ev,Materials);
% [f(7,:),l{7}] = LaRC03(sv,ev,Materials);

% Stitch Strings
% l = reshape([l{:}],size(f'))';

end

% Transform Stresses to Voigt Notation
function [sv] = voigt(s)
% sv1 - s1
% sv2 - s2
% sv3 - s6 (normal 33)
% sv4 - s5 (shear 23)
% sv5 - s4 (shear 13)
% sv6 - s3 (shear 12)
% A = [1 0 0 0 0 0
%      0 1 0 0 0 0
%      0 0 0 0 0 1
%      0 0 0 0 1 0
%      0 0 0 1 0 0
%      0 0 1 0 0 0];
% sv = A*s;
% sv(:,1,:) = s(:,1,:);
% sv(:,2,:) = s(:,2,:);
% sv(:,3,:) = s(:,6,:);
% sv(:,4,:) = s(:,5,:);
% sv(:,5,:) = s(:,4,:);
% sv(:,6,:) = s(:,3,:);

sv = s(:,[1 2 6 5 4 3],:);

end
function [s] = unvoigt(sv)
% sv1 - s1
% sv2 - s2
% sv3 - s6 (shear 12)
% sv4 - s5 (shear 13)
% sv5 - s4 (shear 23)
% sv6 - s3 (normal 33)
% A = [1 0 0 0 0 0
%      0 1 0 0 0 0
%      0 0 0 0 0 1
%      0 0 0 0 1 0
%      0 0 0 1 0 0
%      0 0 1 0 0 0];
% sv = A*s;
% s(:,1,:) = sv(:,1,:);
% s(:,2,:) = sv(:,2,:);
% s(:,6,:) = sv(:,3,:);
% s(:,5,:) = sv(:,4,:);
% s(:,4,:) = sv(:,5,:);
% s(:,3,:) = sv(:,6,:);

s = sv(:,[1 2 6 5 4 3],:);
end

% Calculate Failure Criterias
function [f,g] = TsaiWu(s,Materials,matAllow)
% Calculate Tsai-Wu Failure Criteria
% This is the standard implementation without the special factor and
% biaxial specialization for F12
% Get Material Properties: Assume Transverse Isotropy
M = Materials;
N = size(s,1);
F11t = repmat(matAllow(1,:,:),N,1,1);%M.F11t;
F11c = repmat(matAllow(2,:,:),N,1,1);%M.F11c;
F22t = repmat(matAllow(3,:,:),N,1,1);%M.F22t;
F22c = repmat(matAllow(4,:,:),N,1,1);%M.F22c;
F33t = F22t;
F33c = F22c;
t12 = repmat(matAllow(5,:,:),N,1,1);%M.F12;

% Calculate Factors
F1 = 1./F11t - 1./F11c;
F2 = 1./F22t - 1./F22c;
F3 = 1./F33t - 1./F33c;
F11 = 1./(F11c.*F11t);
F22 = 1./(F22c.*F22t);
F33 = 1./(F33c.*F33t);
F44 = 1./t12.^2;
F55 = 1./t12.^2;
F66 = 1./t12.^2;

% Calculate Interaction Terms
F12 = -0.5*sqrt(abs(F11.*F22));
F13 = -0.5*sqrt(abs(F11.*F33));
F23 = -0.5*sqrt(abs(F22.*F33));

% Compute Failure Index
s1 = s(:,1,:);
s2 = s(:,2,:);
s3 = s(:,3,:);
s4 = s(:,4,:);
s5 = s(:,5,:);
s6 = s(:,6,:);

f = F1.*s1 + F2.*s2 + F3.*s3 + ...
    F11.*s1.^2 + F22.*s2.^2 + F33.*s3.^2 + F44.*s4.^2 + F55.*s5.^2 + F66.*s6.^2 + ...
    2*F12.*s1.*s2 + 2*F13.*s1.*s3 + 2*F23.*s2.*s3;


% A = [F11 F12 F13
%      F13 F22 F23
%      F13 F23 F33]
% eig(A)
% [F1 F2 F3]

% Compute Gradient of Failure Index
g = [F1 + 2*F11.*s1 + 2*F12.*s2 + 2*F13.*s3,...
     F2 + 2*F22.*s2 + 2*F12.*s1 + 2*F23.*s3,...
     F3 + 2*F33.*s3 + 2*F13.*s1 + 2*F23.*s2,...
     2*F44.*s4,...
     2*F55.*s5,...
     2*F66.*s6];

% Return label
% v = f;
% l = repmat({'Tsai-Wu - 1 Mode'},1,size(f,2));

end

function [f,l,v] = maxStrain(s,Materials)
% Calculate Max Strain Failure Criteria
% Get Parameters
% Get Material Properties: Assume Transverse Isotropy
M = Materials;
S11c = M.S11c;
S11t = M.S11t;
S22c = M.S22c;
S22t = M.S22t;
S33c = S22c;
S33t = S22t;
S12 = M.S12;
S13 = S12;
S23 = S12;

% Get Strains
s1 = s(1,:);
s2 = s(2,:);
s3 = s(3,:);
s23 = s(4,:);
s13 = s(5,:);
s12 = s(6,:);

% Calculate criteria for s1
ind = s1 >= 0;
F(1,ind)  = abs(s1(ind))./S11t;
F(1,~ind) = abs(s1(~ind))./S11c;

% Calculate criteria for s2
ind = s2 >= 0;
F(2,ind)  = abs(s2(ind))./S22t;
F(2,~ind) = abs(s2(~ind))./S22c;

% Calculate criteria for s3
ind = s3 >= 0;
F(3,ind)  = abs(s3(ind))./S33t;
F(3,~ind) = abs(s3(~ind))./S33c;

% Calculate criteria for s12
ind = s12 >= 0;
F(4,ind)  = abs(s12(ind))./S12;

% Calculate criteria for s13
ind = s13 >= 0;
F(5,ind)  = abs(s13(ind))./S13;

% Calculate criteria for s23
ind = s23 >= 0;
F(6,ind)  = abs(s23(ind))./S23;

% Determine Failure
[f,iF] = max(F,[],1);
v = F;

% Return Label
l = {'Fiber Tension','Fiber Compression','Matrix Tension','Matrix Compression','Transverse Tension','Transverse Compression',...
     'In-Plane Shear','Transverse Shear In Fiber','Transverse Shear in Matrix'};
l = l(iF);
end

function [f,l,v] = maxStress(s,Materials)
% Calculate Max Stress Failure Criteria
% Get Material Properties: Assume Transverse Isotropy
M = Materials;
S11c = M.F11c;
S11t = M.F11t;
S22c = M.F22c;
S22t = M.F22t;
S33c = S22c;
S33t = S22t;
S12 = M.F12;
S13 = S12;
S23 = S12;

% Get Strains
s1 = s(1,:);
s2 = s(2,:);
s3 = s(3,:);
s23 = s(4,:);
s13 = s(5,:);
s12 = s(6,:);

% Calculate criteria for s1
ind = s1 >= 0;
F(1,ind)  = abs(s1(ind))./S11t;
F(1,~ind) = abs(s1(~ind))./S11c;

% Calculate criteria for s2
ind = s2 >= 0;
F(2,ind)  = abs(s2(ind))./S22t;
F(2,~ind) = abs(s2(~ind))./S22c;

% Calculate criteria for s3
ind = s3 >= 0;
F(3,ind)  = abs(s3(ind))./S33t;
F(3,~ind) = abs(s3(~ind))./S33c;

% Calculate criteria for s12
ind = s12 >= 0;
F(4,ind)  = abs(s12(ind))./S12;

% Calculate criteria for s13
ind = s13 >= 0;
F(5,ind)  = abs(s13(ind))./S13;

% Calculate criteria for s23
ind = s23 >= 0;
F(6,ind)  = abs(s23(ind))./S23;

% Determine Failure
[f,iF] = max(F,[],1);
v = F;

% Return Label
l = {'Fiber Tension','Fiber Compression','Matrix Tension','Matrix Compression','Transverse Tension','Transverse Compression',...
     'In-Plane Shear','Transverse Shear In Fiber','Transverse Shear in Matrix'};
l = l(iF);
end

function [f,l,v] = Hashin(s,Materials)
% Calculate Hashin Failure Criteria
% Get Material Properties
% Assume Transverse Isotropy
M = Materials;
F11t = M.F11t;
F11c = M.F11c;
F22t = M.F22t;
F22c = M.F22c;
F33t = M.F22t;
F33c = M.F22c;
F12 = M.F12;
F23 = F12;
F13 = F12;

% Get Stresses
s1 = s(1,:);
s2 = s(2,:);
s3 = s(3,:);
s12 = s(6,:);
s13 = s(5,:);
s23 = s(4,:);

% Initialize Failure Indices
F = zeros(6,size(s,2));

% Calculate Tensile/Compressive Fiber Failure
ind = s1>=0;
f1 = (s1./F11t).^2 + (s12.^2 + s13.^2)./F12^2;
f2 = (s1./F11c).^2;
F(1,ind) = f1(ind);
F(2,~ind) = f2(~ind);

% Calculate Tensile/Compressive Matrix Failure
ind = s2 + s3 >= 0;
f1 = (s2 + s3).^2./F22t^2 + (s23.^2 - s2.*s3)./F23^2 + (s12.^2 + s13.^2)./F12^2;
f2 = ((F22c./(2*F23)).^2 - 1)*(s2 + s3)./F22c + (s2 + s3).^2./(4*F23^2) + (s23.^2 - s2.*s3)./F23^2 + (s12.^2 + s13.^2)./F12^2;
F(3,ind) = f1(ind);
F(4,~ind) = f2(~ind);

% Calculate Interlaminar Tensile/Compressive Matrix Failure
ind = s3 >= 0;
f1 = (s3./F33t).^2;
f2 = (s3./F33c).^2;
F(5,ind) = f1(ind);
F(6,~ind) = f2(~ind);

% Take Maximums
[f,iF] = max(F,[],1);
v = F;

% Return Labels
l = {'Fiber Failure in Tension','Fiber Failure in Compresion',...
     'Matrix Failure in Tension','Matrix Failure in Compression',...
     'Interlaminar Tensile Failure','Interlaminar Compressive Failure'};
l = l(iF);
end

function [f,l,v] = Christensen(s,Materials)
% Calculate Hashin Failure Criteria
% Get Material Properties
% Assume Transverse Isotropy
M = Materials;
F11t = M.F11t;
F11c = M.F11c;
F22t = M.F22t;
F22c = M.F22c;
F33t = M.F22t;
F33c = M.F22c;
F12 = M.F12;
F23 = F12;
F13 = F12;

% Get Stresses
s1 = s(1,:);
s2 = s(2,:);
s3 = s(3,:);
s12 = s(6,:);
s13 = s(5,:);
s23 = s(4,:);

% Calculate Matrix Controlled Failure
f1 = (1./F22t - 1./F22c)*(s2 + s3) + 1./(F22t*F22c)*(s2 + s3).^2 + 1./F23.^2*(s23.^2 - s2.*s3) + 1./F12^2*(s12.^2 + s13.^2);

% Calculate Fiber Controlled Failure
f2 = (1./F11t - 1./F11c)*s1 + 1./(F11t*F11c)*s1.^2;

% Take Maximums
[f,iF] = max([f1;f2],[],1);
v = [f1;f2];

% Return Labels
l = {'Matrix Controlled Failure','Fiber Controlled Failure'};
l = l(iF);
end

function [f,l,v] = Puck(s,e,Materials)
% Calculate Puck Failure Criteria
% Get Material Properties
% Assume Transverse Isotropy
M = Materials;
F11t = M.F11t;
F11c = M.F11c;
F22t = M.F22t;
F22c = M.F22c;
F33t = M.F22t;
F33c = M.F22c;
F12 = M.F12;
F23 = F12;
F13 = F12;
S11c = M.S11c;
S11t = M.S11t;
E11f = M.E11f;
v12f = M.v12f;

% Get Stresses and Strains
s1 = s(1,:);
s2 = s(2,:);
s3 = s(3,:);
s12 = s(6,:);
s13 = s(5,:);
s23 = s(4,:);
e1 = e(1,:);
e2 = e(2,:);
e3 = e(3,:);
e12 = e(6,:);
e13 = e(5,:);
e23 = e(4,:);

% Specify Puck Inclination Parameters and Magnification Parameters
% Assumed for CFRP for now
pvpp = 0.35;
pvpm = 0.3;
ms1 = 1.1;

% Initialize Failure Modes
F = zeros(5,size(e,2));

% Calculate Tensile/Compressive Fiber Failure Mode
S = e1 + v12f./E11f*ms1*(s2 + s3);
ind = S >= 0;
f1 = S./S11t;
f2 = -S./S11c + (10*e12).^2;
F(1,ind) = f1(ind);
F(2,~ind) = f2(~ind);

% Calculate Inter Fiber Failure Mode Parameters
RvvA = F12./(2*pvpm)*(sqrt(1 + 2*pvpm*F22c./F12)-1);
pvvm = pvpm*RvvA./F12;
t12c = F12*sqrt(1 + 2*pvvm);

% Calculate Interfiber Mode Conditions
indA = s2 >= 0;
indB = s2 < 0 & abs(s2./s12) <= RvvA./t12c;
indC = s2 < 0 & abs(s2./s12) > RvvA./t12c;

% Calculate Mode A: Matrix Failure in Transverse Tension
fA = sqrt((s12./F12).^2 + (1 - pvpp*F22t./F12)^2*(s2./F22t).^2) + pvpp*s2./F12;

% Calculate Mode B: Matrix Failure in Moderate Transverse Compression
fB = 1./F12*(sqrt(s12.^2 + (pvpm*s2).^2) + pvpm*s2);

% Calculate Mode C: Matrix Failure in Large Transverse Compression
fC = ((s12./(2*(1 + pvvm)*F12)).^2 + (s2./F22c).^2).*(-F22c./s2);

% Store
F(3,indA) = fA(indA);
F(4,indB) = fB(indB);
F(5,indC) = fC(indC);

% Take Maximums
[f,iF] = max(F,[],1);
v = F;

% Return Labels
l = {'Tensile Fiber Failure Mode','Compressive Fiber Failure Mode','Mode A: Matrix Failure in Transverse Tension',...
     'Mode B: Matrix Failure in Moderate Transverse Compression','Mode C: Matrix Failure in Large Transverese Compression'};
l = l(iF);
end

function [f,l,v] = LaRC03(s,e,Materials)
% Calculate LaRC03 Failure Criteria
% Get Material Properties
% Assume Transverse Isotropy
M = Materials;
F11t = M.F11t;
F11c = M.F11c;
F22t = M.F22t;
F22c = M.F22c;
F33t = M.F22t;
F33c = M.F22c;
F12 = M.F12;
F23 = F12;
F13 = F12;
S11c = M.S11c;
S11t = M.S11t;
E11f = M.E11f;
v12f = M.v12f;
G12 = M.G12;
v21 = M.v21c;
E11 = 1/2*(M.E11c + M.E11t);
E22 = 1/2*(M.E22c + M.E22t);

% Get Stresses and Strains
s1 = s(1,:);
s2 = s(2,:);
s3 = s(3,:);
s12 = s(6,:);
s13 = s(5,:);
s23 = s(4,:);
e1 = e(1,:);
e2 = e(2,:);
e3 = e(3,:);
e12 = e(6,:);
e13 = e(5,:);
e23 = e(4,:);

% Calculate Method Parameters
XT = F11t;
XC = F11c;
YT = F22t;
YC = F22c;
SL = F12;
ao = 53*pi/180;
nL = -SL*cos(2*ao)/(YC*cos(ao)^2);
nT = -1/tan(2*ao);
ST = YC*cos(ao)*(sin(ao) + cos(ao)/tan(2*ao));

% Calculate Fracture Toughness Properties: Assume Thickply
L220 = 2*(1/E22 - v21^2/E11);
L440 = 1/G12;
YisT = 1.12*sqrt(2)*YT;
SisL = sqrt(2)*SL;
g = 1.12^2*L220/L440*(YT/SL)^2;

% Calculate Angles
phiC = atan((1 - sqrt(1 - 4*(SisL./XC + nL)*(SisL./XC)))./(2*(SisL./XC + nL)));
phi = (abs(s12) + (G12 - XC)*phiC)./(G12 + s1 - s2);

% Calculate stresses in the fiber misalignment frame
s1m = cos(phi).^2.*s1 + sin(phi).^2.*s2 + 2*sin(phi).*cos(phi).*s12;
s2m = sin(phi).^2.*s1 + cos(phi).^2.*s2 - 2*sin(phi).*cos(phi).*s12;
s12m = -sin(phi).*cos(phi).*s1 + sin(phi).*cos(phi).*s2 + (cos(phi).^2 - sin(phi).^2).*s12;

% Calculate Effective Shear Stresses. Do Grid Point Maximization over a
a = linspace(-10,80,50)'*pi/180;
ls = ones(1,size(s,2));
la = ones(size(a,1),1);
aa = a*ls;
ss2 = la*s2;
ss12 = la*s12;
ss2m = la*s2m;
ss12m = la*s12m;
tauT = ReLU(-ss2.*cos(aa).*(sin(aa) - nT*cos(aa)));
tauL = ReLU(cos(aa).*(abs(ss12) + nL*ss2.*cos(aa)));
tauTm = ReLU(-ss2m.*cos(aa).*(sin(aa) - nT*cos(aa)));
tauLm = ReLU(cos(aa).*(abs(ss12m) + nL*ss2m.*cos(aa)));

% Calculate Fiber Failure In Tension (s1 > 0)
ind1 = s1 >= 0;
f1 = e1./S11t; % Max Strain Criterion

% Calculate Fiber Failure in Compression 1 (s1 < 0 and s2m < 0)
ind2 = s1 < 0 & s2m < 0;
f2 = ReLU((abs(s12m) + nL*s2m)./SisL);

% Calculate Fiber Failure in Compression 2 (s1 < 0 and s2m >= 0)
ind3 = s1 < 0 & s2m >= 0;
f3 = (1-g)*(s2m./YisT) + g*(s2m./YisT).^2 + (s12m./SisL).^2;

% Calculate Matrix Cracking in Tension (s2 > 0)
ind4 = s2 > 0;
f4 = (s2./YT).^2 + (s12./SL).^2;

% Matrix Cracking in Biaxial Compression (s2 < 0 and s1 < -YC)
ind5 = s2 < 0 & s1 < -YC;
f5 = (tauTm./ST).^2 + (tauLm./SL).^2;
f5 = max(f5,[],1); % Grid point max over a

% Matrix Cracking in Compression (s2 < 0 and s1 >= -YC)
ind6 = s2 < 0 & s1 >= -YC;
f6 = (tauT./ST).^2 + (tauL./SL).^2;
f6 = max(f6,[],1); % Grid point max over a

% Store Failure Criteria
F = zeros(6,size(s,2));
F(1,ind1) = f1(ind1);
F(2,ind2) = f2(ind2);
F(3,ind3) = f3(ind3);
F(4,ind4) = f4(ind4);
F(5,ind5) = f5(ind5);
F(6,ind6) = f6(ind6);

% Take Maximums
[f,iF] = max(F,[],1);
v = F;

% Return Labels
l = {'Fiber Failure in Tension'
     'Fiber Failure in Compression with Misalignment Normal Compression'
     'Fiber Failure in Compression with Misalignment Normal Tension'
     'Matrix Cracking in Transverse Tension'
     'Matrix Cracking in Compression with Biaxial Compression'
     'Matrix Cracking in Compression with Transverse Compression'}';
l = l(iF);

end
% Rectified Linear Unit for LarC03
function y = ReLU(x)
y = x.*(x>=0);
end