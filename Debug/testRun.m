% Get application data
function testRun()
% Single test run (pair with "Run and Time" to diagnos performance
% bottlenecks

% Get application data
DATA = getappdata(1,'DATA');

% Initialize Design and Digital Thread
D0 = initializeDigitalThread();

% Do for each policy
tic
for ii = [2]%1:size(DATA.Policies,2)
    P = DATA.Policies(ii);
    bellmanBackup(D0.copy,P.policy,ii);
end
toc
disp('Policies Optimized')


% % Copy Policies
% DATA = getappdata(1,'DATA');
% 
% DATA.Policies(2).policy = DATA.Policies(1).policy.copy;
% 
% setappdata(1,'DATA',DATA) 
end

function D0 = initializeDigitalThread()
DATA = getappdata(1,'DATA');
G = DATA.GeometryMaker;
design = Design(G.region,G.inputDisturbances);
D0 = DigitalThread(0,0,design);
D0.getStatistics();
end

function DATA = initializeData(flg)
% Initialize a Digital Thread Trajectory
for jj = 1:2
for ii = 1:1:1
    Trajectories(jj).Trajectory(ii) = DigitalThread(ii-1,ii-1,Design(Region(),InputDisturbances(1,1,1,1,1)));
end
end
% Specify Policies
% p = {{'E','E','D','D'}
%      {'E','D','E','D'}
%      {'E','D','D','E'}
%      {'D','E','E','D'}
%      {'D','D','E','E'}
%      {'D','E','D','E'}};
% p = {{'E','D','D'}
%      {'D','E','D'}};

p = {{'E','D','D'}
     {'E','D','D'}
     {'E','D','D'}
     {'D','E','D'}
     {'D','E','D'}
     {'D','E','D'}};
m = {'G','P','S','G','P','S'};
d = [0 1.0 1.0 0 1.0 1.0];
ind = [3 6];
sp = 1 + 0*ind;

% sp = (-0.25*log(1./linspace(0.1,0.9,length(d))-1));
% ind = d*0+1;

p = p(ind);
m = m(ind);
d = d(ind);

for ii = 1:size(p,1)
    Policies(ii).Trajectories = Trajectories;
    Policies(ii).rules = p{ii};
    Policies(ii).discounts = d(ii);
    Policies(ii).modes = m{ii};
    Policies(ii).Costs = [];
    Policies(ii).MConstraints = [];
    Policies(ii).SConstraints = [];
    Policies(ii).ESU = [];
    Policies(ii).SensorProbability = sp(ii);
end
DATA.Policies = Policies;

% Get Reference to Main Panel
if flg
    mP = findobj('tag','mainPanel');
    DATA.PanelID = mP;
    DATA.CurrentMenuOption = {'Geometry'};
end
end
% 
function TEMP = initializeTEMP()
TEMP.Margins.Z.k = [];
TEMP.Margins.Z.dk = [];
TEMP.Margins.Z.ddk = [];
TEMP.Margins.A(100).k = [];
TEMP.Margins.A(100).dk = [];
TEMP.Margins.A(100).ddk = [];
TEMP.Margins.A(100).dfk = [];
TEMP.Margins.A(100).fk = [];

TEMP.Complexity.Z.k = [];
TEMP.Complexity.Z.dk = [];
TEMP.Complexity.Z.ddk = [];
TEMP.Complexity.A.dfk = [];
TEMP.Complexity.A.fk = [];
TEMP.Complexity.A.k = [];
TEMP.Complexity.A.dk = [];
TEMP.Complexity.A.ddk = [];

TEMP.FEM.Z.k = [];
TEMP.FEM.Z.dk = [];
TEMP.FEM.Z.ddk = [];
TEMP.FEM.A.k = [];
TEMP.FEM.A.dk = [];
TEMP.FEM.A.ddk = [];

TEMP.Quad.k = [];
TEMP.Quad.dk = [];
TEMP.Quad.ddk = [];

end

function kernelCheck()
clc
x = randn(1,2);

y = randn(5,2);
v = 1;
yL = 1;
z = randn(5,1);
dx = 1e-4;

[k,dk,ddk] = kernel(x,y,yL,v,1,z);

for ii = 1:size(x,2)
    xp = x;
    xp(ii) = xp(ii) + dx;
    [kp,dkp] = kernel(xp,y,yL,v,1,z);
    g(:,ii) = (kp-k)*z/dx;
    H(:,:,ii) = (dkp-dk)/dx;  
end

[dk',g']
ddk
H
end

function policyCheck()

% [U] = computeControl(obj,digitalThread,v)
clc
% % rng default
DATA = getappdata(1,'DATA') 
DT = DATA.Policies(3).Trajectories(1).Trajectory(2);
policy = DATA.Policies(3).policy;

% a = randn(size(policy.coefficients(2).A,1),1);
t = 2; 
Co = policy.coefficients(t).C;
policy.coefficients(t).C = rand(size(Co))*1e6;

% [rank(policy.coefficients(t).X),size(policy.coefficients(t).X)]
% 
% [rank(policy.coefficients(t).Y),size(policy.coefficients(t).Y)]

U = policy.computeControl(DT,[],1,0);
a = U.a+randn(size(U.a));
U = policy.computeControl(DT,a,1,1);

dv = 1e-5*1i;
ind = 1:3;
ind = ind' + [0,cumsum(U.No(1:3))];
ind = [ind(:)]';

for ii = ind
    vp = a;
    vp(ii) = vp(ii) + dv;
    
    Up = policy.computeControl(DT,vp,1,1);
    dZ(ii) = imag(Up.Z - U.Z)./imag(dv);
    dZZ(:,ii) = imag(Up.dZ-U.dZ)./imag(dv);
        
end
B = U.dZZ;%transpose(U.kx)*U.dZZ*U.kx;
num2cell([U.dZ(ind),dZ(ind)'])
% num2cell(U.dZZ(ind,ind))
% num2cell(dZZ(ind,ind))
num2cell((B(ind,ind))./(dZZ(ind,ind)'+1e-12))



policy.coefficients(t).C = Co;


end