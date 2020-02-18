classdef Region < handle
% About: Class definition for Region
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% /////////////////////////////////////////////////////////////////////////
    properties
        curves                 % Boundary Curves
        box                    % Bounding Box
        field                  % Displacement Field
        boundaryConditions     % Boundary Conditions
        laminate               % Laminate Definition
        sensor                 % Sensor Locations
        complexity             % Object Complexity
        margins                % Margins
        designVariables        % Design Variables
        basis                  % Basis for FEM
        mesh                   % Mesh Definition
    end
    
    methods
        % Class constructor
        function obj = Region()
                                
        end
        
        % Class Copy
        function region = copy(obj)
           % Copy parameters
           region = Region();
           region.curves = obj.curves;
           region.box = obj.box;
           region.field = obj.field;
           region.boundaryConditions = obj.boundaryConditions;
           region.laminate = obj.laminate;
           region.sensor = obj.sensor;
           region.complexity = obj.complexity;
           region.margins = obj.margins;
           region.designVariables = obj.designVariables;
           region.basis = obj.basis;
           region.mesh = obj.mesh;
            
        end
        
        % Region Update
        function update(obj,G)
            % Find bounding curve
            for ii = size(G.Z,2):-1:1
                % Compute Maximum Distance
                Z = G.Z{ii};
                r(ii) = max(sum(Z.^2,2));
            end
            [~,indr] = max(r);
            
            % Loop through curves
            X = [];
            I = [];
            for ii = 1:1:size(G.Z,2)                
                % Get Points
                B(ii).x = G.Y{ii};
                B(ii).I = G.I{ii};
                
                % Get Design Variable Updates
                X = [X;B(ii).x(1:end-1,:)];
                I = [I;B(ii).x(1:end-1,1)*0+ii];
                
                % Assign Types and Color Plots
                if indr == ii
                    B(ii).type = 'Fill';
                    B(ii).plotColor = 'bo-';
                else
                    B(ii).type = 'Hole';
                    B(ii).plotColor = '.-';
                end
                
            end
            
            % Update Design Variables
            obj.designVariables.xs = X;
            
            % Compute Reference Map
%             if size(B,2) ~= size(obj.curves,2)
                % Build Reference Map Coordinates
                obj.buildReferenceMapCoordinates(B);
                                
                % Build Maps
                obj.buildMaps();
                
                % Reset mesh
                obj.mesh = [];

%             end
                        
        end
        
        % Reference Map Builder
        function C = buildReferenceMapCoordinates(obj,B)                         
            % Loop through curves and determine mapping coordinates
            U = [];
            for ii = 1:1:size(B,2)                              
                % Inherit B
                C(ii).x = B(ii).x;
                C(ii).type = B(ii).type;
                C(ii).plotColor = B(ii).plotColor;
                C(ii).curveType = B(ii).I;
            end
            
            % Update Object
            obj.curves = C;
            
        end
        
        % Build Maps
        function buildMaps(obj)                        
            % Loop through curves
            V = [];
            
            for ii = 1:size(obj.curves,2)
                % Parameterize Curve
                c = obj.curves(ii);
                N = size(c.x,1);
                s = linspace(0,1,N)';
                V = [V;c.x];
                
                % Solve System
                yS = 0.01;
                v = 0.0;
                [K,dK,ddK] = kernel1(s,s,yS,v,1,[]);
                A = [K s s.^2
                     dK(1,:)-dK(end,:) 0 -2 % Match Tangents
                     ddK(1,:)-ddK(end,:) 0 0]; % Match Curvature;
                Kinv = A\eye(size(A,1));
                
                % Store Mapping Parameters
                obj.curves(ii).mapping.param.s = s;
                obj.curves(ii).mapping.param.yS = yS;
                obj.curves(ii).mapping.param.v = v;
                obj.curves(ii).mapping.param.Kinv = Kinv;
                                                
            end
            obj.box =  [min(V,[],1);max(V,[],1)];
                                                                                       
        end
        
        % Helper Functions
        function [x,xT,xN,t,xTdt] = computeTN(obj,t,ii)
            % Compute Tangents and Normals
            c = obj.curves(ii);
            xs = c.x;
            
            % Compute Map
            if isempty(t)
                t = linspace(0,1,size(xs,1))';
            end
           
            [x,T,Tdt] = obj.mappingC(t,xs,c.mapping.param);
            
            % Normalize
            xTn = sqrt(sum(T.^2,2));
            xT = bsxfun(@times,T,1./xTn);
                        
            % Compute Normalized Derivative wrt arc length
            xTdt = bsxfun(@times,1./(xTn).^2,Tdt - bsxfun(@times,sum(xT.*Tdt,2),xT));
                        
            % Compute Normals
            xN = xT*[0 -1;1 0];
            
            % Reverse Normal Direction if Curve is Hole
            if strcmp(c.type,'Hole')
                xN = -xN;
            end
            
        end
                        
        function [f,df,ddf,fk,dfk,ddfk] = mappingC(obj,t,x,param)
            % Get Parameters
            s = param.s;
            yS = param.yS;
            v = param.v;
            Kinv = param.Kinv;

            % Get kernel and derivatives
            [k,dk,ddk] = kernel1(t,s,yS,v,1,[]);
                                    
            % Compute Mapping
            Kinv = Kinv(:,1:end-2);
            l = ones(size(k,1),1);
            o = l*0;
            k = [k, t t.^2];
            dk = [dk, l 2*t];
            ddk = [ddk, o 2*l];
            Kx = Kinv*x;
            f = k*Kx;
            df = dk(:,:,1)*Kx;
            ddf = ddk(:,:,1)*Kx;
            fk = k*Kinv;
            dfk = dk(:,:,1)*Kinv;
            ddfk = ddk(:,:,1)*Kinv;
              
        end
        function [f,df,d2f] = mappingA(obj,G)
            % Get Parameters
            abar = obj.designVariables.a;
            
            % Compute Transformation and Derivatives
            N = G.N;
            dN1 = G.dN1;
            dN2 = G.dN2;
            dN3 = G.dN3;
            
            % Function
            f = N*abar;

            % First Derivative
            df.dA = N;
            df.d1 = dN1*abar;
            df.d2 = dN2*abar;
            df.d3 = dN3*abar;
                        
            % Second Derivative
            d2f.d1.dA = dN1;
            d2f.d2.dA = dN2; 
            d2f.d3.dA = dN3;               
        end
        
        function [f,df,d2f,d3f] = mappingZ(obj,G)
            % Get Parameters
            param = obj.laminate.paramz;
            b1 = param.b(1);
            b2 = param.b(2);
            b3 = param.b(3);
            zbar = obj.designVariables.z;
            
            % Compute Transformation and Derivatives
            N = G.N*b1;
            dN1 = G.dN1*b1;
            dN2 = G.dN2*b1;
            
            % Function
            f = b3*exp(N*zbar);
            
            % First Derivative
            Sp = @(x) spdiags(x,0,size(x,1),size(x,1));
            df.dZ = Sp(f)*(N);
            df.d1 = f.*(dN1*zbar);
            df.d2 = f.*(dN2*zbar);
            
            % Second Derivative
            d2f.d1.dZ = Sp(df.d1)*(N) + Sp(f)*(dN1);
            d2f.d2.dZ = Sp(df.d2)*(N) + Sp(f)*(dN2); 
            d2f.dZ.dZ = @(v) transpose(Sp(v)*df.dZ)*(N);
            
            % Third Derivative
            d3f.d1.dZ.dZ = @(v) transpose(Sp(v)*d2f.d1.dZ)*(N) + transpose(Sp(v)*df.dZ)*(dN1);
            d3f.d2.dZ.dZ = @(v) transpose(Sp(v)*d2f.d2.dZ)*(N) + transpose(Sp(v)*df.dZ)*(dN2);
            
            % Add Offset
            f = f + b2;
               
        end
                
        % Check Region Interior
        function [I,H] = checkInterior(obj,X)
            I = zeros(size(X,1),1);
            H = I;
            for ii = 1:length(obj.curves)
                % Get region points
                x = obj.curves(ii).x;
                
                % Test if Points are in Polygon of Curve
                ind  = inpolygon(X(:,1),X(:,2),x(:,1),x(:,2));
                
                % Screen based on type
                if strcmp(obj.curves(ii).type,'Hole')
                    ind = ~ind;
                end
                
                % Augment Index
                I = I | ~ind;
                if strcmp(obj.curves(ii).type,'Hole')
                    H = H | ~ind;
                end
                
            end

        end
        
        % Generate Quadrature Points and Weights
        function getMesh(obj,flg)
            % Get Mesh Nodes and Elements Using Dist Mesh
            d = @(u) obj.distanceFunction(u);
            s = 0.015;
            b = @(u) obj.densityFunction(u).*s;
            
            % Compute Fixed Points
            O = [];
            for ii = 1:size(obj.curves,2)
                c = obj.curves(ii);
                O = [O;c.x];
            end
            
            % Generate nodes and elements
            if flg%isempty(obj.mesh)
                % Run distmesh
                [nodes,elements] = distmesh2d(d,b,s,obj.box,O);
                
                % Store distmesh nodes and elements
                triangulation.nodes = nodes;
                triangulation.elements = elements;     
            else
                triangulation = obj.mesh.triangulation; 
            end
            
            % Refine Mesh and Nodes
            % FEM Order (Planar/Boundary, Thickness)
            p.plate = [1 13];
            p.parameter = [1 7];
            p.type = {'Planar/Boundary','Thickness'};
            obj.mesh = refineMesh(triangulation.nodes,triangulation.elements,p,d,obj);
            obj.mesh.triangulation = triangulation;
            
            % Initialize Laminate
            obj.getLaminate()
                       
        end
        
        function [N,dN] = getBasis(obj,x,param)
            % Get param values
            Kinv = param.Kinv;
            
            % Compute
            [k,dk] = interpolator(x,param.order,param.P);
            N = k*Kinv;
            for ii = size(dk,3):-1:1
                dN(:,:,ii) = dk(:,:,ii)*Kinv;
            end
            
        end
        
        function [d] = distanceFunction(obj,U)
            % Calculate Distance Function
            d = 0;
                        
            % Get Box
            B = obj.box;
            
            % Compute Rescale
            u1 = 0.5*(B(1)+B(2));
            u2 = 0.5*(B(3)+B(4));
            L1 = 0.5*(B(2)-B(1));
            L2 = 0.5*(B(4)-B(3));
            L = max(L1,L2);
            
            % Rescale
            U = 1./L.*(U-[u1 u2]);
            
            for ii = 1:length(obj.curves)
                % Get region points
                u = obj.curves(ii).x;
                u = 1./L.*(u-[u1 u2]);
                
                if strcmp(obj.curves(ii).type,'Hole')
                    d = ddiff(d,dpoly(U,u));
                else
                    d = dpoly(U,u); 
                end
   
            end
            
        end
        
        function d = densityFunction(obj,U)
            % Calculate Distance Function
            d = 0;
            s = 1;
            
            % Get Box
            B = obj.box;
            
            % Compute Rescale
            u1 = 0.5*(B(1)+B(2));
            u2 = 0.5*(B(3)+B(4));
            L1 = 0.5*(B(2)-B(1));
            L2 = 0.5*(B(4)-B(3));
            L = [L1,L2];
            
            % Rescale
            U = 1./L.*(U-[u1 u2]);
            
            for ii = 1:length(obj.curves)
                % Get region points
                u = obj.curves(ii).x;
                u = 1./L.*(u-[u1,u2]);
                
                if strcmp(obj.curves(ii).type,'Hole')
                    d = dunion(d,s*dpoly(U,u));
                else
                    d = 1-0*dpoly(U,u); 
                end
            end
            d = (1+s*d);
        end
        
        function [G] = buildGlobalBasis(obj,flg)
            % Compute Global Basis for plate
            HP = obj.getGlobalBasis2D(obj.mesh.plate);
            
            % Compute Global Basis for plate Boundary ---------------------
            % Get Values at Quadrature Point Locations
            M = obj.mesh.plate;
            boundary = M.boundary;
            F2 = M.boundaryElem;
            
            % Get Nodal
            T = M.boundary.T;
            C = M.boundary.C;
            X = M.boundary.X(:,:,1);
            Y = M.boundary.X(:,:,2);
            
            % Get Quadrature point values
            t = T*transpose(F2.N);
            x = X*transpose(F2.N);
            y = Y*transpose(F2.N);
            dx = X*transpose(F2.dN);
            dy = Y*transpose(F2.dN);
            
            % Compute Quadrature
            w = repmat(F2.w',size(t,1),1);
            w = sqrt(dx.^2 + dy.^2).*w;
            
            % Unwrap
            C = repmat(C(:,1),1,size(F2.N,1));
            t = t(:);
            C = C(:);
            w = w(:);
            x = x(:);
            y = y(:);
            
            % Sort in Ascending c and then t
            [~,ind] = sortrows([C,t],[1,2]);
            
            % Store
            HB.tq = t(ind);
            HB.cq = C(ind);
            HB.wq = w(ind);
            HB.xq = [x(ind),y(ind)];
              
            % Compute Basis for Boundary Quadrature
            % Determine boundary quadrature locations in element
            % coordinates
            P2 = obj.mesh.plate.planarElem.param.P;
            EB = M.boundary.elements;
            E = M.elements(EB(:,1),:);
            X = repmat(permute(P2,[3 1 2]),size(E,1),1,1);
            Y = zeros(size(X(:,1,:)));
            for ii = (size(EB,2)-1):-1:1
                Y(:,ii,:) = sum((E == EB(:,1+ii)).*X,2);
            end 

            % Compute Coordinates - Note: Basis are interpolatory
            A =  obj.getBasis(F2.q,obj.mesh.plate.boundaryElem.param);
            for ii = size(X,3):-1:1
                Q(:,ii) = reshape(Y(:,:,ii)*transpose(A),[],1);
            end

            % Evaluate on planar basis functions
            E = M.elements(M.boundary.elements(:,1),:);
            Ni = repmat(E,size(Q,1),1);      
            N = obj.getBasis(Q,M.planarElem.param);
            N = repmat(permute(N,[3 1 2]),size(E,1),1,1);
            N = reshape(N,[],size(N,3));
            
            % Loop through quadrature
            rq = repmat(1:size(Ni,1),1,size(Ni,2));
            cq = Ni;
            vq = N; 
            
            % Determine Sparse Matrix
            Nq = sparse(rq(:),cq(:),vq(:));

            % Store - Need to sort from above
            HB.Nq = Nq(ind,:);
                  
            % Get Node, t index, and ci
            I = zeros(size(M.nodes,1),1);
            TI = I;
            CI = I;
            TI(boundary.edges(:,2)) = boundary.edges(:,3);
            CI(boundary.edges(:,2)) = boundary.edges(:,1);
            
            % Store
            HB.Ni = speye(size(I,1)); % Nodal basis is interpolatary
            HB.ti = TI;
            HB.ci = CI;
            
            % Compute Global Matrix for Ply Angle for plate
            HA = obj.getGlobalBasis2D(obj.mesh.parameter);
            HZ = HA;
            [z,dz] = obj.mappingZ(HZ);
            
            % Augment with Thickness for Ply Angle
            F1 = obj.mesh.plate.thicknessElem;
            q = F1.q;
            N = obj.getBasis(q,obj.mesh.parameter.thicknessElem.param);
            t = kron(q,ones(size(HA.N,1),1));
            HV.N = kron(N,HA.N);
            HV.dN1 = 0;
            HV.dN2 = 0;
            HV.dN3 = 0;
            HV.t = t(:);
            HV.w = kron(0.5*F1.w,HA.w.*z);
            
            % Create Basis for Margin Calculation (on nodes and quadrature)  
            t = kron(q,ones(size(HA.N,1)+size(HA.N,2),1));
            M = [HA.N;speye(size(HA.N,2),size(HA.N,2))];
            HM.N = kron(N,M);
            HQ.N = M;
            HQ.dN1 = 0*M;
            HQ.dN2 = 0*M;
            HM.dN1 = 0;
            HM.dN2 = 0;
            HM.dN3 = 0;
            HM.t = t(:);    
          
            % Compute Global Basis for Ply Angle and Thickness ------------           
            % Augment with Thickness for Ply Angle
            M = obj.mesh.parameter;
            F1 = M.thicknessElem;
            N = HA.N;
            HA.N = kron(F1.N,N); 
            HA.dN1 = kron(F1.N,HA.dN1);
            HA.dN2 = kron(F1.N,HA.dN2);
            HA.dN3 = kron(F1.dN,N);
            HA.w = kron(0.5*F1.w,HA.w.*z); % 0.5*|J|*wx*wt*z

            % Update Jacobian for Ply Angle
            J = repmat(HA.J,1,1,size(F1.w,1)); 
            t = F1.q;
            dz1 = dz.d1;
            dz2 = dz.d2;
            zo = permute([0.5*kron(t,[dz1,dz2]),0.5*kron(t*0+1,z)],[3 2 1]);
            o = zeros(size(zo(1,1,:)));
            J = [J, [o;o]
                 zo];
            Jinv = multinv(J);
            
            % Compute Jacobian Derivatives
            t = reshape(kron(t,ones(size(HZ.N,1),1)),1,1,[]);
            dJ.dz1 = [0 0 0
                      0 0 0
                      0.5 0 0].*t;
            dJ.dz2 = [0 0 0
                      0 0 0
                      0 0.5 0].*t;
            dJ.dz = [0 0 0
                     0 0 0
                     0 0 0.5].*(t*0+1);
                 
            % Compute derivative wrt z
            dJinv.dz1 = -multiprod3(Jinv,dJ.dz1,Jinv);
            dJinv.dz2 = -multiprod3(Jinv,dJ.dz2,Jinv);
            dJinv.dz =  -multiprod3(Jinv,dJ.dz,Jinv);
            if flg
                % Compute Second Derivative wrt z 
                d2Jinv.dz1.dz1 = permute(-multiprod3(dJinv.dz1,dJ.dz1,Jinv) - multiprod3(Jinv,dJ.dz1,dJinv.dz1),[3 1 2]);
                d2Jinv.dz1.dz2 = permute(-multiprod3(dJinv.dz1,dJ.dz2,Jinv) - multiprod3(Jinv,dJ.dz2,dJinv.dz1),[3 1 2]);
                d2Jinv.dz1.dz = permute(-multiprod3(dJinv.dz1,dJ.dz,Jinv) - multiprod3(Jinv,dJ.dz,dJinv.dz1),[3 1 2]);

                d2Jinv.dz2.dz1 = d2Jinv.dz1.dz2;
                d2Jinv.dz2.dz2 = permute(-multiprod3(dJinv.dz2,dJ.dz2,Jinv) - multiprod3(Jinv,dJ.dz2,dJinv.dz2),[3 1 2]);
                d2Jinv.dz2.dz = permute(-multiprod3(dJinv.dz2,dJ.dz,Jinv) - multiprod3(Jinv,dJ.dz,dJinv.dz2),[3 1 2]);

                d2Jinv.dz.dz1 = d2Jinv.dz1.dz;
                d2Jinv.dz.dz2 = d2Jinv.dz2.dz;
                d2Jinv.dz.dz = permute(-multiprod3(dJinv.dz,dJ.dz,Jinv) - multiprod3(Jinv,dJ.dz,dJinv.dz),[3 1 2]);  
            else
                d2Jinv = [];
            end
            % Permute
            dJinv.dz1 = permute(dJinv.dz1,[3 1 2]);
            dJinv.dz2 = permute(dJinv.dz2,[3 1 2]);
            dJinv.dz =  permute(dJinv.dz,[3 1 2]);
            
            % Store
            HA.J = J;
            HA.Jinv = Jinv;
            HA.dJinv = dJinv;
            HA.d2Jinv = d2Jinv;
            HA.dJ = dJ;
            HA.t = t(:);        
            
            % Store All
            G.plate.planar = HP;
            G.plate.boundary = HB;
            G.plate.volume = HV;
            G.angle.volume = HA;
            G.thickness.planar = HZ;
            G.margin.volume = HM;
            G.margin.planar = HQ;
                          
        end
                        
        % Generate Laminate Definition
        function getLaminate(obj) 
            % Specify Thickness Function
            L = 0.02;
            H = 0.05;
            pZ = @(u) 0*u(:,1) + H;
            
            % Get Nodal Locations
            E = obj.mesh.parameter.elements;
            E = unique(E);
            u = obj.mesh.parameter.nodes(E,:);
            
            % Compute Thickness
            b2 = 0.019;
            b1 = 0.5*log((H-b2)/(L-b2));
            b3 = exp(0.5*log((L-b2)*(H-b2)));
            z = pZ(u);
            zi = log((z - b2)/b3)/b1;
                        
            % Determine Box Parameters
            B = obj.box(:);
            u1 = 0.5*(B(1)+B(2));
            u2 = 0.5*(B(3)+B(4));
            L1 = 0.5*(B(2)-B(1));
            L2 = 0.5*(B(4)-B(3));
            
            % Compute Ply Angle
            zq = obj.mesh.parameter.thicknessElem.param.P;
            ind = fullfact([size(u,1),length(zq)]);
            U = [u(ind(:,1),:),zq(ind(:,2))]; 
            ai = angleFunction(obj,[1./L1 1./L2 1].*(U - [u1 u2 0]));
            
            % Store Variables
            obj.designVariables.a = ai;
            obj.designVariables.z = zi;
            
            % Store Parameters
            obj.laminate.parama.x = U;     
            obj.laminate.paramz.b = [b1,b2,b3];
            obj.laminate.paramz.x = u;
            
            % Specify Sensor Locations
            % Build sensor location
            u = (linspace(B(1),B(2),11)-u1)*1.03 + u1;
            v = (linspace(B(3),B(4),10)-u2)*1.03 + u2;
            [U,V] = meshgrid(u,v);
            UV = [U(:),V(:)];
            I = obj.checkInterior(UV);
            
            % Filter on sensor locations away from boundary points
            d = obj.distanceFunction(UV) > -2.2e-2;
            UV((sum(I,2)>0) | d,:) = [];
            
            % Add sensor locations
            obj.designVariables.s = zeros(size(UV));
            obj.designVariables.so = UV;
            obj.designVariables.sa = [L1,L2];
            disp(['Number of sensors: ',num2str(size(UV,1))])
            
            % Add sensor selections
            obj.designVariables.p = ones(size(UV,1),1);
            
            % Add interpolator for sensor measurements
            G = obj.buildGlobalBasis(0);
            y = 20;
            v = 1;
            x = 1./[L1 L2].*(G.plate.planar.x-[u1 u2]);
            K = kernel1(x,x,y,v,0,[]);
            Kinv = (K + 1e-9*eye(size(K)))\eye(size(K));
            obj.sensor.param.X = x;
            obj.sensor.param.y = y;
            obj.sensor.param.v = v;
            obj.sensor.param.Kinv = Kinv;
            
        end
        
        % Sensor probabilities
        function [s,ds,d2s] = sensorProbabilities(obj,p)
            % Specify Scaling Parameter
            a = 4;
            [s,ds,d2s] = sigmoid(a*p);
            
            % Rescale Derivatives
            ds = a*ds;
            d2s = a^2*d2s;
            
            
        end
        
        % Sensor Locations
        function [s] = sensorLocations(obj,p)
            % Just use functional saturation form of sensor probabilities for now
            [s] = obj.sensorProbabilities(p(:));
            
            % Shift and scale
            s = 2*s - 1;
            
            % Reshape and apply bounds
            a = obj.designVariables.sa/10;
            s = a.*reshape(s,size(p)) + obj.designVariables.so;
                   
        end
        
        % Sensor Interpolation
        function [T] = sensorInterp(obj,x)
            % Get Parameters
            X = obj.sensor.param.X;
            y = obj.sensor.param.y;
            v = obj.sensor.param.v;
            Kinv = obj.sensor.param.Kinv;
            
            % Determine Box Parameters
            B = obj.box(:);
            u1 = 0.5*(B(1)+B(2));
            u2 = 0.5*(B(3)+B(4));
            L1 = 0.5*(B(2)-B(1));
            L2 = 0.5*(B(4)-B(3));
            
            % Rescale
            x = 1./[L1 L2].*(x-[u1 u2]);
            
            % Compute
            k = kernel1(x,X,y,v,0,[]);
            T = k*Kinv;
            
        end
        
        % Angle initialization
        function [a] = angleFunction(obj,U)
            t = U(:,3);
%             a = (-pi/2)*(1-t)/2 + (pi/2)*(1+t)/2;
            a = (-pi/2)*(1-t) + (3*pi/2)*(t);
            a = a + 0*0.5*cos(pi*U(:,1)).*sin(pi*U(:,2));
            
        end
        
        % Calculate Global Basis
        function G = getGlobalBasis2D(obj,M)
            % Compute Global Basis for Field
            % N, dN, w, Jinv (q x n)           
            % Evaluate N and dN on the quadrature points
            F1 = M.planarElem;
            N = F1.N;  
            dN = F1.dN;
            w = F1.w;
            N = repmat(permute(N,[3 1 2]),size(M.elements,1),1,1);
            dN = repmat(permute(dN,[4 1 2 3]),size(M.elements,1),1,1);
            w = repmat(w',size(M.elements,1),1);
            NT = reshape(N,[],size(N,3));
            dN = reshape(dN,[],size(dN,3),size(dN,4));
            w = w(:);
            
            % Get elements
            To = M.elements;
            T = repmat(To,size(F1.N,1),1);
                       
            % Determine Sparse Matrix
            rq = repmat(1:size(T,1),1,size(T,2));
            cq = T;
            vq = N;          
            dvq = reshape(dN,[],size(dN,3));
  
            N = sparse(rq(:),cq(:),vq(:));
            dN1 = sparse(rq(:),cq(:),dvq(:,1));
            dN2 = sparse(rq(:),cq(:),dvq(:,2));   
            
            % Compute System Jacobian
            Mn = reshape(M.nodes(To(:),:),size(To,1),size(To,2),[]);
            Mn = repmat(Mn,size(F1.N,1),1,1);
            du1 = permute(sum(dN(:,:,1).*Mn,2),[3 2 1]);
            du2 = permute(sum(dN(:,:,2).*Mn,2),[3 2 1]);
            J  = [du1(1,1,:) du2(1,1,:)
                  du1(2,1,:) du2(2,1,:)];
              
            % Compute Quadrature
            detJ = J(1,1,:).*J(2,2,:) - J(2,1,:).*J(1,2,:);
            w = w.*detJ(:);
     
            % Compute Coordinates
            X = M.nodes(M.elements,1);
            Y = M.nodes(M.elements,2);
            X = reshape(X,size(M.elements));
            Y = reshape(Y,size(M.elements));
            x = X*transpose(F1.N);
            y = Y*transpose(F1.N);
                                                
            % Store (these are on quadrature points)
            G.N = N;
            G.dN1 = dN1;
            G.dN2 = dN2;
            G.w = w;
            G.Jinv = multinv(J);
            G.J = J;
            G.x = [x(:),y(:)];
            G.NT = NT;
            G.dNT = dN;
 
        end
        
        % Calculate Boundary Conditions
        function [LOADS,DISPS] = calculateBoundaryConditions(obj,t)
            % Calculate Boundary Condition Values based on t parameter                        
            % Get Normals and Tangents
            iFill = find(strcmp({obj.curves(:).type},'Fill'));
            [~,xt,xn] = obj.computeTN(t,iFill);
            
            % Loop through all BC's
            BC = obj.boundaryConditions;
            o = zeros(size(xt(:,1,:)));
            l = o+1;
            LOADS = 0;
            DISPS = 0;
            for ii = 1:size(BC,2)
                % Get strings
                loadStr = BC(ii).loads;
                filterStr = BC(ii).filter;
                displacementStr = BC(ii).displacements;
                
                % Evaluate to function
                filter = eval(['@(t) ',filterStr]);
                if strcmp(loadStr,'') || isempty(loadStr)
                    loadStr = '0*[xt(:,[1 1 1 1 1],:)]';
                end
                if strcmp(displacementStr,'') || isempty(displacementStr)
                    displacementStr = '0*[xt(:,[1 1 1 1 1 1 1 1 1 1],:)]';
                end
                loads = eval(['@(t) ',loadStr]);
                displacements = eval(['@(t) ',displacementStr]);
                
                % Evaluate Numerically
                L = loads(t);
                F = filter(t);
                D = displacements(t);
                
                % Store Thus Far
                LOADS = LOADS + bsxfun(@times,F,L);
                DISPS= DISPS + bsxfun(@times,F,D);
                
            end
            
            % Parse Arrays
            DISPS = DISPS(:,:,1);
            LOADS = LOADS(:,:,1);                
        end
                  
    end
        
end