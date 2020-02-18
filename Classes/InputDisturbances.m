classdef InputDisturbances < handle
% About: Class definition for Input Disturbances
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% /////////////////////////////////////////////////////////////////////////
    properties
        loads
        matAllow
        manufacturingParam
        samples
    end
    
    methods
        % Class constructor
        function obj = InputDisturbances(t,L,A,P,MA)           
           % Initialize Gaussian Parameters for loads                      
           obj.loads.u = L;
           obj.loads.S = 1/16*1e12*0.00015*kron(diag([1 1 1 1 1]),kernel1(t,t,2,0,0,[]))*(16);%(0.001*diag(1+0*abs(L))).^2;%
           obj.loads.t = t;
           obj.loads.names = {'N_n (N/m)','N_t (N/m)','M_n (N)','M_t (N)','Q (N/m)'};
           
           % Initialize Gaussian Parameters for Materials
           obj.matAllow.u = A;
           obj.matAllow.S = 1000*diag([100 100 0.05 1.5 0.05])/3*1e12;
           obj.matAllow.names = {'F_{11t} (N/m^2)','F_{11c} (N/m^2)','F_{22t} (N/m^2)','F_{22c} (N/m^2)','F_{12} (N/m^2)'};
           
           % Initialize Gaussian Parameters for Manufacturing Parameters
           obj.manufacturingParam.u = P;
           obj.manufacturingParam.S = (0.10*diag(P)).^2 + 1e-9*eye(size(P,1));
           obj.manufacturingParam.A = MA;

           % Sample Values
           obj.samples.N = 50;
           if length(A) > 1
               obj.sample(obj.samples.N,'all');
           end
           
        end
        
        % Class Copy
        function inputDisturbances = copy(obj)
           % Copy parameters
           inputDisturbances = InputDisturbances(1,1,1,1,1);
           inputDisturbances.loads = obj.loads;
           inputDisturbances.matAllow = obj.matAllow;
           inputDisturbances.manufacturingParam = obj.manufacturingParam;
           inputDisturbances.samples = obj.samples;            
        end
        
        % Sampler
        function sample(obj,M,d)
            % Draw Samples from Gaussian Process
            if strcmp(d,'all') || isfield(d,'strains')
                u = obj.loads.u;
                S = obj.loads.S;
                try
                    cholS = chol(S+1e-9*eye(size(S,1)),'lower');
                catch
                    [V,D] = eig(S);
                    l = diag(D);
                    li = l > 0;
                    S = V(:,li)*D(li,li)*V(:,li)';                    
                    cholS = chol(real(S)+1e-9*eye(size(S,1)),'lower');
                end
                obj.samples.L = bsxfun(@plus,u,cholS*randn(size(S,1),M));
            end
            
            % Draw Samples from Multivariate Gaussian for Materials
            if strcmp(d,'all') || isfield(d,'stresses')
                u = obj.matAllow.u;
                S = obj.matAllow.S;
                cholS = chol(S,'lower');
                obj.samples.A = bsxfun(@plus,u,cholS*randn(size(S,1),M));
                               
            end
            
            % Draw Samples from Multivariate Gaussian for Manufacturing
            % Parameters
            if strcmp(d,'all') || isfield(d,'times')
                u = obj.manufacturingParam.u;
                S = obj.manufacturingParam.S;
                try
                    cholS = chol(S+1e-9*eye(size(S,1)),'lower');
                catch
                    [V,D] = eig(S);
                    l = diag(D);
                    li = l > 0;
                    S = V(:,li)*D(li,li)*V(:,li)';                    
                    cholS = chol(S+1e-9*eye(size(S,1)),'lower');
                end
                obj.samples.P = bsxfun(@plus,u,cholS*randn(size(S,1),M));
                
                % Restrict samples to physical values
                obj.samples.P = (obj.samples.P).*sign(obj.samples.P);
            end
                 
        end
        
        % Update Posterior
        function update(obj,d)            
            % Update Loads
            if isfield(d,'strains') && size(d.strains,1) > 0
                strains = d.strains;
                Lf = d.Lf;
                Qstrain = d.Qstrain;
                u = obj.loads.u;
                S = obj.loads.S;      
                
                K = S'*Lf.A'*((Lf.A*S*Lf.A' + Qstrain)\eye(size(Qstrain,1)));
                obj.loads.u = u + K*(strains-Lf.A*u-Lf.b);
                obj.loads.S = S - K*(Lf.A*S);
                
            end
            
            % Update Material Allowables
            if isfield(d,'stresses')
                stresses = d.stresses;
                Qstress = d.Qstress;
                u = obj.matAllow.u;
                S = obj.matAllow.S;   
                
                K = S'*((S + Qstress)\eye(size(Qstress,1)));
                obj.matAllow.u = u + K*(stresses-u);
                obj.matAllow.S = S - K*(S);
                
            end
            
            % Update Manufacuturing Parameters
            if isfield(d,'times')
                time = d.times;    
                Mf = d.Mf;
                Qtime = d.Qtime;
                u = obj.manufacturingParam.u;
                S = obj.manufacturingParam.S;  

                K = S'*Mf.A'*((Mf.A*S*Mf.A' + Qtime)\eye(size(time,1)));
                obj.manufacturingParam.u = u + K*(time-Mf.A*u);
                obj.manufacturingParam.S = S - K*(Mf.A*S);
                obj.manufacturingParam.A = Mf.A;
            end
            
            % Sample New
            obj.sample(obj.samples.N,d)
        end
        
        % Get Statistics
        function [loads,matAllow,mParam] = getStatistics(obj)
            % Get Mean and Variance for loads
            loads.u = reshape(obj.loads.u,[],5);
            loads.S = sqrt(diag(obj.loads.S));
            loads.SS = obj.loads.S;
            
            % Get upper and lower bounds for loads
            loads.U = reshape(bsxfun(@plus,loads.u(:),loads.S*2),[],5);
            loads.L = reshape(bsxfun(@plus,loads.u(:),-loads.S*2),[],5);
            
            % Return mean and variance for material allowables
            matAllow.u = obj.matAllow.u;
            matAllow.S = sqrt(diag(obj.matAllow.S));
            matAllow.SS = obj.matAllow.S;
            
            % Get upper and lower and bounds for material allowables
            matAllow.U = matAllow.u+matAllow.S*2;
            matAllow.L = matAllow.u-matAllow.S*2;
            
            % Return Mean and Variance for manufacturing parameters
            A = obj.manufacturingParam.A;
            mParam.u = A*obj.manufacturingParam.u;
            mParam.S = sqrt(diag(A*obj.manufacturingParam.S*A'));
            mParam.SS = obj.manufacturingParam.S; 
            mParam.uu = obj.manufacturingParam.u;
            
            % Get upper and lower and bounds for manufacturing parameters
            mParam.U = mParam.u+mParam.S*2;
            mParam.L = mParam.u-mParam.S*2;
            
        end
            
        
        
    end
    
end