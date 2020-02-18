
function tester()
% Function Approximation Parameterization behavior tester
rng default
run1DN()

end

function run1DN()
clc
% Get Grid Points points
D = 3000;
Ngrid = 100;
Xgrid = linspace(-1,1,Ngrid)';
ind = fullfact(Ngrid*ones(1,2));
Xgrid = [Xgrid(ind),zeros(size(ind,1),D-2)];

% Specify Objective
sa = 2;
% f = @(X) sa*prod(cos(1*pi*X(:,[1 2])),2);
f = @(X) sa*cos(1.0*pi*sqrt(sum(X.^2,2)));

% Specify model parameters
Nmod = 2000;
Lskip = 1000;
hX = haltonset(D);
hX = scramble(hX,'RR2');
as = 2*(net(hX,Nmod+Lskip)-0.5);
as = as((Lskip+1):end,:);
g = 1/D;

% Basis Function
KX = @(x) [kernel2(x,as,g,0),ones(size(x,1),1)];

% Specify A matrix
ao = 0;
A = [zeros(1,Nmod+1)];

H = blkdiag(1e-2*speye(Nmod),1e-9);
Hinv = H\speye(size(H));

% Precompute Quantities
N = 1000;
KXgrid = KX(Xgrid);
% Xtrain = (rand(N,D)-0.5)*2.*rand(N,1);%.*[ones(1,2),zeros(1,D-2)];
% Xtest = (rand(N,D)-0.5);
Xtrain = randn(N,D).*rand(N,1);%.*[ones(1,2),zeros(1,D-2)];
Xtest = randn(N,D);

% Xtrain(1,:) = Xtrain(1,:)*0;
% Xtest(1,:) = Xtest(1,:)*0;
% Xtest = Xtrain;
Utrain = f(Xtrain);
Utest = f(Xtest);
KXtrain = KX(Xtrain);
KXtest = KX(Xtest);

% Loop through each sample
flg = 0;
figure(3)
if flg
    % Do Incremental Least Squaers
    for ii = 1:N
        % Get sample point
        u = Utrain(ii,:);

        % Calculate Update
        kj = KXtrain(ii,:);
        dx = u - A*kj' - ao;

        % Update Matrices
        s = 1;%./(1 - 0.7./(ii));
        p = s*Hinv*kj';
        Hinv = s*Hinv - (p*p')./(1 + kj*p);
        v = Hinv*kj';
        A = A + (dx*v');

        % Plot
        plotResults(f,A,ao,D,Xgrid,Xtrain,Utrain,Utest,KXtrain,KXtest,KXgrid,sa,ii,N)
        [max(A),min(A),A*KX(as(1,:)*0)'+ao]
        pause(0.01)
    end
else
    % Just Do Least Squares
    A = (KXtrain'*KXtrain + H)\(KXtrain'*Utrain);
    A = transpose(A);
    plotResults(f,A,ao,D,Xgrid,Xtrain,Utrain,Utest,KXtrain,KXtest,KXgrid,sa,N,N)
end


end

function plotResults(f,A,ao,D,Xgrid,Xtrain,Utrain,Utest,KXtrain,KXtest,KXgrid,sa,ii,N)
switch D
    case 1
        % Plot Test Data
        plot(Xgrid,f(Xgrid),'-'), hold all
        plot(Xgrid,A*KXgrid'+ao,'r.-')

        % Plot Training Points
        plot(Xtrain,f(Xtrain),'o')
        plot(Xtrain,A*KXtrain'+ao,'rs')
        plot(as,as*0,'o'), hold off
        grid on
        axis([-1 1 -10 10])
        
    otherwise
        subplot(4,1,[1 2])
        % Plot Train Data
        R = @(x) reshape(x,sqrt(numel(x))*ones(1,2));
        surf(R(Xgrid(:,1)),R(Xgrid(:,2)),R(A*KXgrid'+ao)), hold all
        
        % Plot Training Points
        surf(R(Xgrid(:,1)),R(Xgrid(:,2)),R(f(Xgrid)),'facealpha',0.1)
        v = A*KXtrain'+ao;
        plot3(Xtrain(:,1),Xtrain(:,2),v,'rs'), hold off
        grid on
        axis([-1 1 -1 1 -sa-2 sa+2])
        
        subplot(4,1,3)
        plot(1:ii,Utrain(1:ii,:),'*-',1:ii,v(:,1:ii),'o-')
        axis([0 N -sa-.5 sa+0.5]), grid on
        
        subplot(4,1,4)
        v = A*KXtest' + ao;
        plot(1:ii,Utest(1:ii,:),'*-',1:ii,v(:,1:ii),'o-')
        ax = axis;
        axis([0 100 ax(3) ax(4)]), grid on
        
        drawnow
end


end


function [d,dd,U] = distanceFunction2(region,V)
% Calculate area collision
U = [];
D = [];
for ii = 1:length(region.curves)
    % Get Points
    u = region.curves(ii).u;
    t = cumsum(sqrt(sum((u(2:end,:)-u(1:end-1,:)).^2,2)));
    t = [0;t/max(t)];
       
    T = linspace(0,1,length(t)*2.5);
    u = interp1(t,u,T);    
    u = u(1:end-1,:);
    
    % Assign conditions
    if strcmp(region.curves(ii).type,'Hole')
        u = [u;mean(u,1)];
        d = [zeros(size(u,1)-1,1);1];
    else
        d = zeros(size(u,1),1);
    end
       
    % Store
    U = [U;u];
    D = [D;d];
    
end

% Build System
y = 0.01;
v = 0;
K = kernel2(U,U,y,v,0,[],[]);
[k,dk] = kernel2(V,U,y,v,0,[],[]);
m = @(x) sum(x.^2,2);

% Calculate distance
a = K\(D-m(U));
d = m(V) + k*a;
dd = 2*V + permute(sum(dk.*transpose(a),2),[1 3 2]);

% Truncate
d(d<0) = 0;
dd(d<0,:) = 0;
% d = log1p(d);
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
