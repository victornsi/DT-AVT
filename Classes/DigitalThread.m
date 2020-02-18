classdef DigitalThread < handle
% About: Class definition for Digital Thread
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% /////////////////////////////////////////////////////////////////////////
    properties
        stage                 % Current Stage
        operation             % Components in operation
        design                % Current Design
        transform             % Feature Vector Transformation 
    end
    
    methods
        % Copy
        function DT = copy(obj)
            % Initialize
            DT = DigitalThread(0,0,obj.design.copy());
            
            % Copy stage
            DT.stage = obj.stage;
            
            % Copy operation
            for ii = 1:size(obj.operation,2)
                operation(ii) = obj.operation(ii).copy();
            end
            if size(obj.operation) == 0
                operation = [];
            end
            DT.operation = operation;
            
            % Copy transform
            DT.transform = obj.transform;
            
        end
        
        % Class constructor
        function obj = DigitalThread(t,n,design) 
            obj.stage = t;
            for ii = 1:n
                operation(ii) = Component(design,ii);
            end
            if n == 0
                operation = [];
            end
            obj.operation = operation;
            obj.design = design;
            obj.transform = [];
        end 
        
        % Methods
        function manufactureAndDeploy(obj,inputs)
            % Put current design into operation
            if isempty(obj.operation)
                obj.operation = Component(obj.design,inputs);
            else
                obj.operation(end+1) = Component(obj.design,inputs);
            end
            
            % Update new design with latest knowledge
            obj.design = obj.operation(end).design.copy();
            
        end
        
        % Feature super vector of design and operation
        function [phi,ao,bo,co,No] = getFeatureVector(obj)            
            % Get statistics from design
            [loads,matAllow,mParam] = obj.design.inputDisturbances.getStatistics;
            
            % Compute Feature Vector
            phi = [loads.u(:);matAllow.u(:);mParam.uu(:)];
            
            % Get Input Transform
            if isempty(obj.transform)
                % Generate Transform
                % Compute Mean
                N = [length(loads.u(:)),length(matAllow.u(:)),length(mParam.uu(:))];
                obj.transform.u = phi;
                obj.transform.N = N;
                obj.transform.i = {(1:N(1)),...
                                   (1:N(2)) + N(1),...
                                   (1:N(3)) + N(1) + N(2)};
                
                % Extract Cholesky Decomposition of Covariance
                SS = blkdiag(loads.SS,matAllow.SS,mParam.SS);
%                 LL = chol(loads.SS,'lower');
%                 LA = chol(matAllow.SS,'lower');
%                 LM = chol(mParam.SS,'lower');
                
                % Extract Diagonals of Covariance
                LL = diag(sqrt(diag(loads.SS)));
                LA = diag(sqrt(diag(matAllow.SS)));
                LM = diag(sqrt(diag(mParam.SS)));
                
                % Compute Eigenvectors and Eigenvalues for Sampling
                [VL,lL] = eig(loads.SS);
                [VA,lA] = eig(matAllow.SS);
                [VM,lM] = eig(mParam.SS);
                
                % Assemble with Sigma Bounds
                L = blkdiag(3*LL,3.5*LA,3*LM);
                                   
                % Store
                obj.transform.V = blkdiag(VL,VA,VM);
                obj.transform.l = [diag(lL);diag(lA);diag(lM)];
                obj.transform.L = L;
                obj.transform.Linv = L\speye(sum(N));
                obj.transform.S = SS;
                                                              
            end
            
            % Compute Feature Vector of Mean
            phi = obj.transform.Linv*(phi - obj.transform.u);
            
            % Augment with Covariance Feature Vector
            S = blkdiag(loads.SS,matAllow.SS,mParam.SS);
            P = (sqrt(diag(S)./diag(obj.transform.S)) - 1);
            phi = 0.9*[phi;P];

            % Get Output Size
            dV = obj.design.region.designVariables;
            No = [size(dV.a,1),size(dV.z,1),size(dV.s,1)*2,size(dV.p,1)];
            
            % Get current design variables  
            ao = [dV.a
                  dV.z
                  dV.so(:)*0
                  dV.p];  
            bo = log(3.5);
            co = bo;   
            
        end                      
        
    end
    
end