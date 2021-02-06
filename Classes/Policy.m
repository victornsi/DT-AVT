classdef Policy < handle
% About: Class definition for Policy
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% /////////////////////////////////////////////////////////////////////////
    properties
        horizonLength         % Horizon length
        coefficients          % Coefficients
        binaryAction          % Binary Action (Fixed)
        mode                  % Policy Mode
        y                     % Discount Factor
        Nb                    % Number of basis functions
        temp                  % Training Temp Variables
    end
    
    methods     
        % Copy Class
        function P = copy(obj)
            % Initialize
            P = Policy(0,0,0,0);
            
            % Copy horizonLength
            P.horizonLength = obj.horizonLength;
            
            % Copy Coefficients
            P.coefficients = obj.coefficients;
            
            % Copy binaryAction
            P.binaryAction = obj.binaryAction;
            
            % Copy Mode
            P.mode = obj.mode;
            
            % Copy y
            P.y = obj.y;
            
            % Copy Nb
            P.Nb = obj.Nb;
            
            % Copy temp
            %P.temp = obj.temp;
                             
        end
        
        % Class constructor
        function obj = Policy(templateTrajectory,binaryAction,discount,modeName) 
            % Specify horizon length
            obj.horizonLength = size(templateTrajectory,2)-1;
            
            % Specify Policy Binary Action
            obj.binaryAction = binaryAction;
                        
            % Specify Discount
            obj.y = discount;
            
            % Specify Mode
            obj.mode = modeName;
            
            % Iniitialize temp
            obj.temp = [];
            
            % Sample Supports and Initialize Coefficients
            if templateTrajectory ~= 0
                obj.sampleSupports(templateTrajectory)
            end
            
        end   
        
        % Sample Supports
        function sampleSupports(obj,templateTrajectory)
            % Specify Number of Basis Functions
            obj.Nb = 2000;
                                       
            % Specify coefficient matrices
            for ii = 1:obj.horizonLength
                % Get feature vector at this stage
                [phi,ao,bo,co,No] = obj.featureVector(templateTrajectory(ii));
                                             
                % Calculate initial coefficients
                obj.coefficients(ii).ao = ao;
                obj.coefficients(ii).bo = bo;
                obj.coefficients(ii).co = co;

                % Initialize Matrices
                N = size(phi,1);
                ar = [zeros(sum(No(1:3)),1);ao(end)*ones(No(4),1)];
                obj.coefficients(ii).A = [zeros(length(ao),obj.Nb),ar];
                obj.coefficients(ii).B = zeros(1,obj.Nb+1);
                obj.coefficients(ii).C = zeros(1,obj.Nb+1);
                obj.coefficients(ii).H = [];
                obj.coefficients(ii).I = [];
                obj.coefficients(ii).n = 0;
                
                % Initialize Information Space and Control Supports
                if ii == 1
                    XY = getappdata(1,'XY');
                    if isempty(XY)
                        % Initialize Information Space Supports
                        Lskip = 1000;
                        hX = haltonset(N);
                        hX = scramble(hX,'RR2');
                        X = 2*(net(hX,obj.Nb+Lskip)-0.5);
                        X = X((Lskip+1):end,:);

                        % Initialize Control Space Supports
                        hY = haltonset(length(ao));
                        hY = scramble(hY,'RR2');
                        Y = 2*(net(hY,obj.Nb+Lskip)-0.5);
                        Y = Y((Lskip+1):end,:);
                        
                        % Store
                        XY.X = X;
                        XY.Y = Y;
                        setappdata(1,'XY',XY)
                        
                    else
                        % Use Previous
                        X = XY.X;
                        Y = XY.Y;
                    end
                      
                end
                             
                % Store
                obj.coefficients(ii).X = X;
                obj.coefficients(ii).Y = Y;
                
                % Create Training Storage Variables
                obj.coefficients(ii).KJ = [];
                obj.coefficients(ii).HJ = [];
                obj.coefficients(ii).VT = [];
                obj.coefficients(ii).VTT = [];
                obj.coefficients(ii).AT = [];
                obj.coefficients(ii).TT = [];
                obj.coefficients(ii).PT = [];
                obj.coefficients(ii).ZP = [];
            end
              
            % Add Regularization
            obj.computeControl(templateTrajectory(1),[],0,0);
            
        end
        
        % Calculate control
        function [U] = computeControl(obj,digitalThread,a,hflg,hessflg)
            % Get feature vector
            [phit,~,~,~,No] = obj.featureVector(digitalThread);
            region = digitalThread.design.region;
            t = digitalThread.stage+1;

            % Apply Discount
            U.y = obj.y;
            if t == obj.horizonLength
                U.y = 0;
            end

            % Sensor Placement and Selection Flag
            mb = obj.mode; 
            ub = obj.binaryAction{t};
            U.us = 1;
            if ~strcmp(mb,'S') || (U.y == 0 || strcmp(ub,'E'))
                U.us = 0;
            end
                        
            % Get Variables
            X = obj.coefficients(t).X;
            Y = obj.coefficients(t).Y;
            H = obj.coefficients(t).H;
            ao = obj.coefficients(t).ao;
            bo = obj.coefficients(t).bo;
            co = obj.coefficients(t).co;
                        
            % Specify scaling parameters
            % For a
            N = digitalThread.transform.N;
            g1 = sqrt([1/N(1) + zeros(N(1),1)
                       1/N(2) + zeros(N(2),1)
                       1/N(3) + zeros(N(3),1)]);
            g1 = [g1;g1]';
            U.g1 = g1;
            gc1 = 1;
                      
            % For Z
            g2 = sqrt([1/No(1) + zeros(No(1),1)
                       1/No(2) + zeros(No(2),1)
                       1/No(3)*strcmp(mb,'P') + zeros(No(3),1)
                       1/No(4)*U.us + zeros(No(4),1)])';
            U.g2 = g2;
            gc2 = 1;
            
            % Calculate magnitude adjustments
            L = [ones(No(1),1)*0.75
                 ones(No(2),1)*0.2
                 ones(No(3),1)
                 ones(No(4),1)];
            Linv = 1./L;
            sd = 2; 
                                                       
            % Compute Projection. This ensures that components that should
            % not be active for a particular policy is forced to be so.
            if strcmp(ub,'E') || (U.y == 0)
                % Adjust for Experiment and Final Stage
                kx = ([ones(No(1) + No(2),1);zeros(No(3),1);zeros(No(4),1)]);
            else
                % Adjust for Sensor Placement and Selection
                kx = ([ones(No(1) + No(2),1);strcmp(mb,'P')*ones(No(3),1);U.us*ones(No(4),1)]);
                ao = ao.*[ones(sum(No(1:3)),1);(~U.us).*ones(No(4),1)];
            end
            
            % Initialize H and I matrices
            if isempty(H)
                obj.applyRegularization(X,Y,g1,g2,gc1,gc2);
            end
                        
            % Get A, B, and C
            A = obj.coefficients(t).A;
            B = obj.coefficients(t).B;
            C = obj.coefficients(t).C;
            
            % Compute Basis Function
            kj = kernel2(transpose(phit).*g1,X.*g1,gc1,0);
            kj = [kj(:);1];
            
            % Calculate Control
            if isempty(a)
                a = A*kj;
            end
            U.a = a;           
            u = kx.*a + ao;
           
            % Calculate Value Function
            V = B*kj + bo;
     
            % Calculate Control Basis
            Ckj = transpose(C).*kj;
            s = Linv.*(u-ao);
            
            % Pass through saturation. Need to do this because future value
            % function approximation is designed on a biunit hypercube.
            % Points outside this hypercube cannot be adequately captured
            % and degrades performance
            [thet,dthet,d2thet] = sigmoid(sd*s);
            thet = 2*thet - 1;
            dthet = 2*dthet*sd;
            d2thet = 2*d2thet*sd.^2;
            
            % Do if future contributes
            if U.y > 0 && hflg == 1
                dX2 = sum(((transpose(phit) - X).*g1).^2,2);
                [hj,dhj1,d2hj] = kernel(transpose(thet).*g2,Y.*g2,gc2,hessflg,transpose(C(1:end-1)),transpose(dX2));
                dhj1 = squeeze(dhj1);
                dhj = dhj1.*dthet.*Linv.*g2'; 
                d2hj = (dthet.*Linv.*g2').*(d2hj.*transpose(dthet.*Linv.*g2')) + diag(dhj1.*d2thet.*Linv.^2.*g2');
            else
                hj = zeros(obj.Nb,1);
                dhj = 0*Linv;      
                d2hj = diag(dhj);
            end
            % Construct basis    
            hj = [hj(:);1];
               
            % Build Functions 
            U.Z = C*hj + co;
            U.dZ = dhj;
            U.dZZ = d2hj;
                        
            % Store values in struct
            U.u = u;
            U.t = t;
            U.ub = ub;
            U.V = V;
            U.ao = ao;
            U.bo = bo;
            U.co = co;
            U.L = L;                                                
            U.kj = kj;
            U.hj = hj;
            U.phit = phit;
            U.thet = thet;
            U.No = No;           
            U.kx = kx;
            U.Ckj = Ckj;
                       
            % Compute Sensor Selection Cost
            p = u(sum(No(1:3)) + (1:No(4)));
            U.p = sum(region.sensorProbabilities(p))/length(p);
            sPw = 0.005/length(p)*10*5/3;
            [sP,dsP,d2sP] = region.sensorProbabilities(p);
                        
            % Add Sensor Selection Cost to immediate cost
            U.Y = sPw*sum(sP);
            U.dY = ([zeros(sum(No(1:3)),1);sPw*dsP]);
            U.dYY = diag([zeros(sum(No(1:3)),1);sPw*d2sP]);
            
            % Add Barrier to Control Components
            du = u-ao;
            sb = 2;
            U.P = sum(((du>=sb).*(du-sb).^2 + (du<=-sb).*(du+sb).^2))/sum(No)*100;
            U.dP = 2*((du>=sb).*(du-sb) + (du<=-sb).*(du+sb))/sum(No)*100;
            U.dPP = 0;%2*((du>=sb) + (du<=-sb))/sum(No)*100; This hessian leads to stability issues so leave as zero
                        
            % Adjust Z and Y based on high level decision and stage
            if strcmp(ub,'E') || t == obj.horizonLength
                U.Z = sum(Ckj) + co;
                U.dZ = U.dZ*0;
                U.dZZ = U.dZZ*0;
                U.hj = kj;

                U.Y = 0;
                U.dY = U.dY*0;
                U.dYY = 0;
            end
                                        
        end
        
        % Compute feature vector
        function [phi,ao,bo,co,No] = featureVector(obj,digitalThread)
            % Extract Statistics Information from Digital Thread
            % Feature vector for now is just a vector of statistics
            [phi,ao,bo,co,No] = digitalThread.getFeatureVector;   
            
        end
        
        % Apply Regularization
        function applyRegularization(obj,X,Y,gx,gy,gc1,gc2)  
            % Specify Matrices
            H = blkdiag(1e2*speye(obj.Nb),1e9);
            I = blkdiag(1e2*speye(obj.Nb),1e9);
                        
            % Store
            for ii = 1:obj.horizonLength
                obj.coefficients(ii).H = H;
                if ii < obj.horizonLength
                    obj.coefficients(ii).I = I;
                end
            end
        end
                
    end
    
end