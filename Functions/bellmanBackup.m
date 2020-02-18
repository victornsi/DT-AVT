function bellmanBackup(D0,policy,pInd)
% About: Approximate Dynamic Programming Algorithm - Bellman Backup
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% /////////////////////////////////////////////////////////////////////////
% Specify number of iterations
NIter = 1;
R = 100;

% Specify Optimization Options
options = optimset('fminunc');
options.Algorithm = 'trust-region';
options.Display = 'iter-detailed';
options.Hessian = 'user-supplied';
options.GradObj = 'on';
options.MaxFunEvals = 50;
options.MaxIter = 5;

% Test Derivative flg
% Use with an EDD policy
flg = 0;

% Do for each iteration
T = size(policy.coefficients,2)-flg;

% Specify cuttoff and sampling parameters
policy.temp.cutoff = 0.25;
policy.temp.N = 3;
policy.temp.q = 0.1;
policy.temp.r = 0.9;

% Determine input variations
iD0 = D0.design.inputDisturbances;

% Store Initial Inputs
initialInputs.iD0 = iD0.copy;
initialInputs.loads = iD0.loads.u;
initialInputs.matAllow = iD0.matAllow.u;
initialInputs.manufacturingParam = iD0.manufacturingParam.u;

for ii = 1:NIter
    % Do for each stage
    for t = T:-1:1
        disp('Status')
        disp([pInd ii NIter t T])
                                      
        % Regenerate full new trajectory based on current policy
        if t == T
            [trajectory] = generateTrajectory(D0,policy,initialInputs);
        end         
        D = trajectory(t);
                        
        % Get Control
        F = size(policy.coefficients(t).VTT,2) <= R*5;
        U = policy.computeControl(D,[],1,0);
        U.NIter = NIter;
        U.R = R;
        U.flg = flg;
        
        % Initialize Cost, Constraint, and Hessian Functions
        f = @(a) costFunction(a,D,policy,U,F,initialInputs);
        c = @(a) constraintFunction(a,D,policy);
%         options.OutputFcn = @(a,optimValues,state) outfun(a,optimValues,state,D,policy); % Uncomment for stage diagnostics and visualization
        
        % Initalize control
        at = U.a;
                                             
        % Test derivatives
        if flg
            testGradientc(at,c,U)
            testGradientf(at,trajectory(t),policy,U)
            pause
        end
                    
        % Run Optimization for control
        counter = 0;
        SF = 1;
        ct = 0;
        while counter < 3 && SF~=-1       
            % Run optimizer
            [at,~,SF] = fminunc(f,at,options);
            ct = c(at);
            if max(ct(1)) > policy.temp.cutoff
                counter = counter + 1;
            else
                counter = 3;
            end
        end
        
        % Update Function Approximation
        updateFunctionApproximation(at,SF,ct,D,policy,R)
               
    end
    
end

end

function Vt = fastObjective(D,U)
% Evaluate objective 
Ns = D.design.inputDisturbances.samples.N;
r_t = stageCosts(D,U,1:Ns,'');
[Vt] = r_t + U.y*exp(U.Z) + U.Y;
end

function [stop] = outfun(x,optimValues,state,D,policy)
stop = false;
switch state
    case 'init'
        if optimValues.fval > 10000
            stop = true;
            disp('Bad Sample')
        end
    case {'iter','done'}       
        a = x;
        U = policy.computeControl(D,a,1,0);        
        if U.t == 2 && optimValues.iteration > -1%&& strcmp('done',state)%
            % Plot Angles and Thickness
            figure(10)
            clf(10)
            mP = uipanel(10,'units','normalized','position',[0.13 0.01 0.865 0.98],...
                'tag','mainPanel');

            iD = D.design.inputDisturbances;
            ind = randi(50);
            inputs.loads = iD.samples.L;
            inputs.matAllow = iD.samples.A;
            inputs.manufacturingParam = iD.samples.P(:,ind);

            % Update
            [an,z,s,p] = parseInput(U.u,D.design.region);
            D.design.region.designVariables.a = an;
            D.design.region.designVariables.z = z;
            D.design.region.designVariables.s = s;
            D.design.region.designVariables.p = p;

            CV = ComponentVisualizer(mP);
            CV.plot(Component(D.design,inputs),{'Input Disturbances'})
            colormap(linspecer)
            drawnow
            pause(1)

        end

end
end

% Update function approximations
function updateFunctionApproximation(at,SF,ct,D,policy,R)
t = D.stage + 1;
if SF ~= -1 && (max(ct(1))) < policy.temp.cutoff    
    % Get Iteration
    n = policy.coefficients(t).n;

    % Get costs
    U = policy.computeControl(D,at,1,0);
    Vt = fastObjective(D,U);
    
    % Get parameters
    kj = U.kj;
    A = policy.coefficients(t).A;
    B = policy.coefficients(t).B;
    H = policy.coefficients(t).H;
    
    % Update Training Data
    if size(policy.coefficients(t).VT,2) <= R
        policy.coefficients(t).KJ = [policy.coefficients(t).KJ, kj];
        policy.coefficients(t).VT = [policy.coefficients(t).VT, Vt];
        policy.coefficients(t).AT = [policy.coefficients(t).AT, at];
        policy.coefficients(t).PT = [policy.coefficients(t).PT, U.phit];
%         policy.coefficients(t).TT = [policy.coefficients(t).TT, U.thet];
        policy.coefficients(t).ZP = [policy.coefficients(t).ZP, exp(B*kj + U.bo)];
    end
    
    % Get Errors
    da = at - (A*kj);
    dV = log(Vt) - (B*kj + U.bo);

    % Compute Preliminary quantities
    p = H*kj;
    H = H - (p*p')./(1 + kj'*p);
    v = H*kj;
    
    % Update coefficients
    policy.coefficients(t).A = A + da*v';
    policy.coefficients(t).B = B + dV*v';
    policy.coefficients(t).H = H;
    policy.coefficients(t).n = n + 1;
    
end

end

function [f,df,d2f] = costFunction(a,D,policy,Uo,flg,initialInputs)
% Get Control
hessFlg = 1;
getFutureValueFunction(D,a,policy,Uo,policy.temp.N*(strcmp(Uo.ub,'E')+1),flg,initialInputs);
U = policy.computeControl(D,a,1,hessFlg);

% Get parameters
y = U.y;

% Compute stage costs
Ns = D.design.inputDisturbances.samples.N;
[r_t,dr_t,d2r_t] = stageCosts(D,U,1:Ns,'derivative');

% Compute Stage Constraint
[c,~,dc,~,margins] = constraintFunction(a,D,policy);
D.design.region.margins = margins;

% Specify weights
l = [10,0];
if ~Uo.flg
    % Ease into the discount (for good initialization)
    n = policy.coefficients(U.t).n;
    y = y*(n > 0);
    
    % Do for each constraint
    for ii = 1:length(c)
        % Specify Constraint Smoother
        cs = 0.5*max(0,c(ii)).^2;
        dcs = max(0,c(ii))*dc(:,ii);
        d2cs = dc(:,ii)*max(0,sign(c(ii)))*transpose(dc(:,ii));
        
        % Add constraint
        r_t = r_t + l(ii)*cs;
        dr_t = dr_t + l(ii)*dcs;
        d2r_t = d2r_t + l(ii)*d2cs;     
    end
end

% Store values
fR = r_t + (y>0)*U.Y;
df1 = dr_t + (y>0)*U.dY;
d2f = d2r_t + (y>0)*U.dYY;

% Add future value and barrier
f = fR + y*exp(U.Z) + U.P;
df = df1 + y*exp(U.Z)*(U.dZ) + U.dP;
d2f = d2f + y*exp(U.Z)*(U.dZ*transpose(U.dZ) + U.dZZ)*hessFlg + U.dPP;

% Verify derivative magnitudes
if y > 0
%     z = exp(U.Z)*(U.dZ);
%     [da1,dz1,ds1,dp1] = parseInput(df1 + U.dP,D.design.region);
%     [za,zz,zs,zp] = parseInput(z,D.design.region);
%     num2cell([da1'*da1,za'*za
%               dz1'*dz1,zz'*zz
%               ds1(:)'*ds1(:),zs(:)'*zs(:)
%               dp1'*dp1,zp'*zp])
        
end

% Transform back into a coordinates
df = transpose(U.kx.*df);
d2f = (U.kx.*(d2f.*transpose(U.kx)));

% Parse control
[~,~,s,p] = parseInput(U.a,D.design.region);
s = sqrt(sum(s.^2,2));
s = [max(s),min(s)];
p = [min(p),max(p)];

% Display costs
disp(['Cost      : ',sprintf('%8.4f ',fR,y,exp(U.Z),U.Y,U.P,s(:),U.p,p(:))])
disp('----------------------------------------')

end

% Get Future Value Function
function getFutureValueFunction(Di,at,policy,U,Ns,flg,initialInputs)
% Get Time index
t = Di.stage+1;

% Only do for specific stages or settings
if t < policy.horizonLength && U.y > 0 && ~U.flg
    % Initialize
    V = zeros(1,Ns);
    
    % Get Coefficients
    C = policy.coefficients(t).C;
    I = policy.coefficients(t).I;
        
    % Calculate Control Variations
    A = perturbControl(at,U,policy,Ns,1);
    A(:,1) = at;
        
    % Loop through samples
    HJ = [];
    TT = [];
    for jj = 1:Ns
        % Compute Control
        a = A(:,jj);
        
        % Get Thread
        D = Di.copy;
        
        % Sample under Digital Thread and Control
        W = policy.computeControl(D,a,1,0);

        % Compute Transition
        Dp = oneStepForward(D,policy,a,initialInputs);
        
        % Compute Forward Value Function Sample
        S = policy.computeControl(Dp,[],0,0);
             
        % Compute Value Function Sample
        V(jj) = (S.V);

        % Get Basis
        hj = W.hj;
           
        % Get Errors
        dZ = V(jj) - (C*hj + W.co);
                
        % Compute Preliminary Quantities
        p = I*hj;
        I = I - (p*p')./(1 + hj'*p);
        v = I*hj;
        HJ = [HJ,hj];
        TT = [TT,W.thet];
        
        % Update Coefficients
        C = C + dZ*v';      
    end

     % Store updates
     policy.coefficients(t).C = C;
     policy.coefficients(t).I = I;
          
     % Store Points
     if flg
         policy.coefficients(t).VTT = [policy.coefficients(t).VTT exp(V(:))'];
         policy.coefficients(t).TT = [policy.coefficients(t).TT, TT];
         policy.coefficients(t).HJ = [policy.coefficients(t).HJ, HJ];
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
a = u(1:Na,:);
z = u(Na+(1:Nz),:);
s = reshape(u(Na+Nz+(1:Ns),:),[],2,size(u,2));
p = u(Na + Nz + Ns + (1:Np),:);

end

function s = weightFunction(x)
% Weight function for learning rates
% N = 10;
% bn = 0.05; % Initial Slope
% yn = 0.9; % Target value at N
% 
% an = (1./(1-yn) - bn*N - 1)/N.^2;
% s = 1 - 1./(an.*x.^2 + bn.*x + 1);
s = 1./(x + 1);

end

function [c,ceq,dc,dceq,margins] = constraintFunction(a,D,policy)
% Get Control
U = policy.computeControl(D,a,0,0);

% Compute stage constraints
[g_t,dg_t,margins] = stageConstraints(D,U);

% Build c and dc
c = g_t;
dc = dg_t;
ceq = 0;
dceq = 0*a;

% Transform back into b coordinates
disp('----------------------------------------')
disp(['Constraint: ',sprintf('%8.4f ',c(:))])

end

% Trajectory generation
function [D] = generateTrajectory(D0,policy,initialInputs)
% Get horizon length
N = size(policy.coefficients,2)-1;
   
% Initialize Trajectory
D(1) = D0.copy;
for ii = 1:N
    % Vary initial condition
    if ii == 1
        D(ii) = perturbThread(D(ii),policy);         
    end  
        
    % Perturb Sensor Selection for better sampling
    U = policy.computeControl(D(ii),[],1,0);
    at = U.a;
    if U.us > 0
        at = perturbControl(at,U,policy,1,0); 
%         U = policy.computeControl(D(ii),at,1,0);
%         policy.coefficients(ii).TT = [policy.coefficients(ii).TT, U.thet];
    else
        at = [];
    end
    
    % Transition state
    D(ii+1) = oneStepForward(D(ii),policy,at,initialInputs);
            
end

end

function D = perturbThread(D,policy)
% Determine input variations
iD = D.design.inputDisturbances;

% Sample Initial State
if policy.coefficients(1).n > 0
    % Adjust Means - Periodically Sample Near Origin
    ii = D.transform.i;
    o = rand(1,1);
    q = 0.2*o;
    L = D.transform.L;
    du = L*(2*rand(size(L,2),1)-1);
    
    iD.loads.u = iD.loads.u + q*du(ii{1},1);
    iD.matAllow.u = iD.matAllow.u + q*du(ii{2},1);
    iD.manufacturingParam.u = iD.manufacturingParam.u + q*du(ii{3},1);
    
    % Adjust Variances (method 1)
%     L = 1 + q*(2*rand(size(L,1),1)-1);
%     L1 = L(ii{1});
%     L2 = L(ii{2});
%     L3 = L(ii{3});
%     
%     iD.loads.S = (L1*L1').*iD.loads.S;
%     iD.matAllow.S = (L2*L2').*iD.matAllow.S;
%     iD.manufacturingParam.S = (L3*L3').*iD.manufacturingParam.S;
    
    % Adjust Variances (method 2)
    ii = D.transform.i;
    V = D.transform.V;
    L = 1 + q*(2*rand(size(L,1),1)-1);
    S = V*diag((L).^2.*D.transform.l)*V';
    
    iD.loads.S = S(ii{1},ii{1});
    iD.matAllow.S = S(ii{2},ii{2});
    iD.manufacturingParam.S = S(ii{3},ii{3});
       
end

% Resample
iD.sample(iD.samples.N,'all')

end

function at = perturbControl(at,U,policy,Ns,flg)
% Perturb Control
% Local perturbation plus global for sensor selection
q = policy.temp.q;
r = policy.temp.r;
dU = 2*rand(size(at,1),Ns)-1;
dUs = 2*rand(1,Ns)-1;

% Perturb Sensor Selection (global)
v = [zeros(U.No(1),1)
     zeros(U.No(2),1)
     zeros(U.No(3),1)
     ones(U.No(4),1)];
if flg
    % Perturb all control components
    w = [ones(U.No(1),1)
         ones(U.No(2),1)
         ones(U.No(3),1)
         ones(U.No(4),1)];
else
    % Only perturb sensor selection
    w = v;
end
q = q.*w;
r = r.*v;
at = (1 - r).*at + r.*((U.kx.*U.L).*dUs); 
at = (1 - q).*at + q.*((U.kx.*U.L).*dU);

end

function Dp = oneStepForward(D,policy,a,initialInputs)
% Get action
U = policy.computeControl(D,a,0,0);

% Sample an input to use
iD = D.design.inputDisturbances;
ind = randi(iD.samples.N);
w = weightFunction(policy.coefficients(U.t).n);
s = (1-w);

% Choose an easy case to ease into general inputs
inputs.loads = iD.samples.L(:,ind)*s + (1-s)*initialInputs.loads;
inputs.matAllow = iD.samples.A(:,ind)*s + (1-s)*initialInputs.matAllow*1.2;
inputs.manufacturingParam = iD.samples.P(:,ind)*s + (1-s)*initialInputs.manufacturingParam;
    
% Transition dynamics
Dp = transitionModel(D,U.ub,U.u,inputs);

end

% Test Gradient
function testGradientf(b,D,policy,Uo)
% Get initial values
t = D.stage+1;
Co = policy.coefficients(t).C;
policy.coefficients(t).C = Co + randn(size(Co));
U = policy.computeControl(D,b,0,0);

f = @(a) costFunction(a,D,policy,Uo,0);
[g,dg,d2g] = f(b);
db = 1e-6;

L = 1:2;
L = L' + [0,cumsum(U.No(1:end-1))];
L = [L(:)]';

for ii = L
    % Get new direction
    bp = b;
    bp(ii) = bp(ii)+db;
    
    % Get transitions
    [gp,dgp] = f(bp);
    dgn(ii,1) = (gp-g)/(db);
    d2gn(:,ii) = (dgp-dg)/(db);
end

% Compare
disp('Cost Gradient')
num2cell([transpose(dg(L)),dgn(L)])

disp('Cost Hessian')
num2cell(d2g(L,L))
num2cell(d2gn(L,L))
num2cell((d2g(L,L)./(d2gn(L,L)+1e-12)))
f(b);
policy.coefficients(t).C = Co;
end

% Test Gradient
function testGradientc(b,f,U)
% Get initial values
[g,~,dg] = f(b);
db = 1e-6;
dg = U.kx(:).*dg;
L = [(1:2),(1:2)+U.No(1),(1:2)+U.No(1)+U.No(2),(1:2)+U.No(1)+U.No(2)+U.No(3)];
for ii = L
    % Get new direction
    bp = b;
    bp(ii) = bp(ii)+db;

    % Get transitions
    gp = f(bp);
    dgn(ii,:) = (gp-g)/(db);
end

% Compare
disp('Constraint Gradient')
num2cell([dg(L,:),dgn(L,:)])
f(b);
end