function policyTest
% Load Policy
DATA = getappdata(1,'DATA');


% Train
for ii = [2]
    figure(1+1)
    policy = DATA.Policies(ii).policy;
    D = policy.coefficients;
    train(policy,D);
end


% setappdata(1,'DATA',DATA)

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

function train(policy,D)
% Do for each stage
for t = 1:3
    % Get Coefficients
    KJ = D(t).KJ;
    HJ = D(t).HJ;
    AT = D(t).AT;
%     CT = D(t).CT;
    VTT = D(t).VTT;
    VT = D(t).VT;
    PT = D(t).PT;
    TT = D(t).TT;
    ZP = D(t).ZP;
    
    % Truncate
    ind = [1:size(KJ,2)];
    KJT = KJ(:,ind);
%     KJT(1:end-1,:) = KJT(1:end-1,:).^(1)*sqrt(policy.Nb).^(1-1);
        
%     AT(:,ind) = [];
%     VT(:,ind) = [];  
%     ZT(:,ind) = [];
    
    % Train
%     M = blkdiag(1e-2*eye(size(KJ,1)-1),1e-9);
% 
%     % Least Squares
%     policy.coefficients(t).A = ((KJT*KJT' + M)\(KJT*AT(:,ind)'))';
%     policy.coefficients(t).B = ((KJT*KJT' + M)\(KJT*(log(VT(ind)) - policy.coefficients(t).bo)'))';

    H = max(abs(policy.coefficients(t).A*KJT-AT(:,ind)),[],1);
    G = [VT(ind);exp(policy.coefficients(t).B*KJT + policy.coefficients(t).bo)];
    if size(HJ,2) > 0
        I = [VTT;exp(policy.coefficients(t).C*HJ + policy.coefficients(t).co)]';
    else
       I = 0; 
    end
   
    N = 5;
    subplot(N,3,t)
    plot(1:length(G),G,'.-'), grid on
    
    subplot(N,3,3+t)
    plot(1:size(PT(:,:),1),PT(:,:)), grid on
    
    subplot(N,3,6+t)
    plot(1:size(AT,1),AT'), grid on
    
    
    subplot(N,3,9+t)
    plot(1:size(I,1),I,'.-'), grid on
    
    subplot(N,3,12+t)
    plot(1:size(TT(:,:),1),TT(:,:),'-'), grid on
    
    
%     subplot(4,3,9+t)
%     plot(1:length(H),max(abs(AT),[],1),'.-'), grid on
    drawnow
    
    
end

end

% M = 1;
% x = linspace(-M,M,100);
% y = [x/3
%     sigmoid(4*x)-0.5
%      asinh(2*x)/3];
% 
% figure(2)
% plot(x,y,'.-')
% axis equal