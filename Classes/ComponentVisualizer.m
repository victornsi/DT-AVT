classdef ComponentVisualizer < handle
% About: Class definition for Component Visualizer
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% /////////////////////////////////////////////////////////////////////////
    properties
        parentPanel     % Parent panel for visualizer
        component       % Current component
    end
    
    methods
        % Class constructor
        function obj = ComponentVisualizer(parentPanel) 
            obj.parentPanel = parentPanel;
        end       
        
        % Plot Functions
        function plot(obj,component,plotType)
            % Plot based on requested plot type
            if size(plotType,2)>1
                p1 = plotType{1};
                p2 = plotType{2};
            else
                p1 = plotType{1};
                p2 = 'Reference';
            end
            
            switch  p1
                case {'Geometry',''}
                    obj.plotGeometry(component);
                case 'Input Disturbances'
                    obj.plotInputDisturbances(component);
                case 'Mesh'
                    obj.plotMesh(component);
                case 'Laminate'
                    obj.plotLaminate(component,p2);
                case 'Material Length'
                    obj.plotMaterialLength(component)
                case 'Complexity'
                    % Run Compute Complexity
                    obj.plotComplexity(component);
                case 'Displacements'
                    % Displacements
                    obj.plotDisplacements(component);
                case 'Strains'
                    % Run FEM
                    obj.plotStrains(component);
                case 'Margins'
                    % Run Margins
                    obj.plotMargins(component);
           
            end
            
        end
        function plotGeometry(obj,component)
            % Update component
            obj.component = component;
            
            % Plot Geometry
            mP = obj.parentPanel;
            
            % Clear Children
            delete(findobj(mP,'tag','Colorbar'))
            delete(get(mP,'children'))
            
            % Create an Axes of appropriate shape
            positions = obj.getPositions(4);
%             pos = [positions{1}(1:2),positions{2}(1)-positions{1}(1) + positions{1}(3),positions{1}(4)];
            ca = axes('units','normalized','position',positions{1});
                        
            % Reset uicontext menu
            set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'),'parent',mP)
            
            % Plot Thickness
            R = component.design.region;
            P = R.mesh.plate;
            b = R.laminate.paramz.b;
            XY = R.laminate.paramz.x;
            Z = b(2) + b(3)*exp(b(1)*R.designVariables.z);
            trisurf(P.elements(:,1:3),XY(:,1),XY(:,2),Z,'edgecolor','none','parent',ca,'facealpha',0.5)
            colorbar('peer',ca,'location','Eastoutside')
            shading(ca,'interp')
                         
            hold(ca,'on')
            title(ca,'Thickness','fontunits','normalized'), xlabel(ca,'x (m)','fontunits','normalized'), ylabel(ca,'y (m)','fontunits','normalized')
            view(ca,[0 90]), axis(ca,'tight')
            set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'))
            set(ca,'fontunits','normalized')
            
            
            % Plot Sensor Geometry
            for ii = 1:size(R.curves,2)
                [x,xt,xn,~,xtd] = R.computeTN([],ii);
                plot(ca,x(:,1),x(:,2),R.curves(ii).plotColor), hold(ca,'all')
                
                % Get current and reference location
                xz = R.sensorLocations(R.designVariables.s);
                xzo = R.designVariables.so;
                
                % Get selection probabilities
                s = R.sensorProbabilities(R.designVariables.p);        
                col = [0 0 0].*s  + (1-s).*[1 1 1]*.75;
                
                % Plot
                plot3(ca,xz(:,[1 1])',xz(:,[2 2])',[0*s,s]',':','color',0.7*[1 1 1])
                plot(ca,xzo(:,1),xzo(:,2),'o','color',0.7*[1 1 1])
                scatter3(ca,xz(:,1),xz(:,2),s,81,col,'+','linewidth',1)  
%                 plot(ca,xz(:,1),xz(:,2),'k+','linewidth',0.75)   
                 
            end
  
            % Turn on grid
            hold(ca,'off')
            grid(ca,'on')
            title(ca,'Geometry, Sensor Locations, and Thickness','fontunits','normalized')
            xlabel(ca,'x (m)','fontunits','normalized')
            ylabel(ca,'y (m)','fontunits','normalized')
            set(ca,'fontunits','normalized')
            
            % Plot Mesh
            ca = axes('position',positions{2},'parent',mP);
            if ~isempty(R.mesh)
                simpplot(R.mesh.triangulation.nodes,R.mesh.triangulation.elements(:,1:3),ca), hold(ca,'all')
            end
                
            % Plot Curves
            t = linspace(0,1-1/500,500)';
            for jj = 1:size(R.curves,2)
                [x,xt,xn,~,xtd] = R.computeTN(t,jj);
                plot(ca,x(:,1),x(:,2),R.curves(jj).plotColor), hold(ca,'all')
                quiver(ca,x(:,1),x(:,2),xn(:,1),xn(:,2))
                quiver(ca,x(:,1),x(:,2),xt(:,1),xt(:,2))
%                 quiver(ca,x(:,1),x(:,2),xtd(:,1),xtd(:,2))
           
            end
            hold(ca,'off')
            grid(ca,'on')
            axis(ca,'tight')
            %                 axis(ca,'equal')
            title(ca,['Configuration Map']),xlabel(ca,'x (m)'),ylabel(ca,'y (m)')
%                         
            % Plot Boundary Conditions
            pos = [positions{3}(1:2),positions{end}([1 3])*[1;1] - positions{3}(1),positions{4}(4)];
            ca = axes('position',pos,'parent',mP);
            if ~isempty(R.boundaryConditions) && ~isempty(R.curves)
                iFill = find(strcmp({R.curves(:).type},'Fill'));
                [x,xt,xn,t] = R.computeTN([],iFill);
                [L,D] = R.calculateBoundaryConditions(t);
                di = sum(D(:,1:5),2);
                
                L = [bsxfun(@times,L(:,1),xn) + bsxfun(@times,L(:,2),xt),...
                    (bsxfun(@times,L(:,3),xt) - bsxfun(@times,L(:,4),xn)),...
                     L(:,5)];
                
                quiver3(ca,x(:,1),x(:,2),x(:,1)*0,L(:,1),L(:,2),L(:,2)*0), hold(ca,'all')
                quiver3(ca,x(:,1),x(:,2),x(:,1)*0,L(:,3),L(:,4),L(:,2)*0)
                quiver3(ca,x(:,1),x(:,2),x(:,1)*0,L(:,1)*0,L(:,2)*0,L(:,5))
                plot(ca,x(:,1),x(:,2),'k.-')
                plot(ca,x(di==3,1),x(di==3,2),'k^','markersize',20,'linewidth',2)
                plot(ca,x(di==5,1),x(di==5,2),'ks','markersize',20,'linewidth',2)
                plot(ca,x(di<3 & di>0,1),x(di<3 & di>0,2),'ko','markersize',15,'linewidth',2)
                hold(ca,'off')
            end
            title(ca,'Boundary Conditions','fontunits','normalized')
            xlabel(ca,'x (m)','fontunits','normalized')
            ylabel(ca,'y (m)','fontunits','normalized'), grid(ca,'on')
            set(ca,'fontunits','normalized')
            view([0 90])
                        
        end   
        function plotMesh(obj,component)
            % Update component
            obj.component = component;
            R = component.design.region;
            
            % Get panel
            mP = obj.parentPanel;
            
            % Clear Children
            delete(findobj(mP,'tag','Colorbar'))
            delete(get(mP,'children'))
            
            % Plot Nodes and Quadrature Points
            names = {'Nodes','Quadrature','Plate-Domain','Plate-Boundary','Plate-Thickness','Parameter-Domain','Parameter-Thickness'};
            positions = {[0.05 .55 .65 .4],[0.05 0.075 0.65 0.4],...
                                           [0.75 0.75 .2 .2],...
                                           [0.73 0.55 0.1 0.14],...
                                           [0.87 0.55 0.1 0.14],...
                                           [0.75 0.275 0.2 0.2],...
                                           [0.75 0.075 0.2 0.14]};
                                       
            % Get Global Basis
            G = R.buildGlobalBasis(0);
  
            for ii = 1:length(positions)
                % Specify Axes
                ca = axes('units','normalized','position',positions{ii});
                                                
                % Switch based on type
                switch names{ii}
                    case 'Nodes'
                        % Get Mesh data
                        M = R.mesh.plate;
                        
                        % Plot Mesh
                        simpplot(M.nodes,M.elements(:,1:3),ca), hold(ca,'all')
                        
                        % Plot Boundaries
                        for jj = 1:size(R.curves,2)
                            c = R.curves(jj);
                            plot(ca,c.x(:,1),c.x(:,2),'-'), hold(ca,'all')
                        end
                        
                        % Plot Orientation
                        nodes = M.nodes;
                        elements = M.elements;
                        edges = M.boundary.edges;
                        eMid = 0;
                        n = 3;
                        for jj = 1:n
                            eMid = eMid + nodes(elements(:,jj),:);
                        end
                        eMid = eMid/n;
                        hold on
                        plot(ca,eMid(:,1),eMid(:,2),'ro')
                        
                        % Calculate Orientations
                        theta = elements;
                        c = {'b','g','m'};
                        for jj = 1:n
                            dx = nodes(elements(:,jj),:)-eMid;
                            p = eMid + 0.5*dx;
                            theta(:,jj) = atan2(dx(:,2),dx(:,1));
                            plot(ca,p(:,1),p(:,2),['.',c{jj}])
                        end
                        
                        plot(ca,nodes(:,1),nodes(:,2),'.k')
                        
                        % Plot boundary nodes
                        ind = edges(edges(:,3)==0,2);
                        plot(ca,nodes(ind,1),nodes(ind,2),'ok','markersize',15)
                        ind = edges(:,2);
                        plot(ca,nodes(ind,1),nodes(ind,2),'^b')
                        
%                         % Plot Curves
%                         for jj = 1:size(R.curves,2)
%                             t = edges(edges(:,1)==jj,3);
%                             c = R.curves(jj);
%                             u =  R.mappingC(t,c.x,c.mapping.param);
%                             plot(u(:,1),u(:,2),'*')
%                         end
                        xlabel(ca,'x (m)','fontunits','normalized'),ylabel(ca,'y (m)','fontunits','normalized')

                    case 'Quadrature'
                        % Plot Mesh
                        M = R.mesh.plate;
                        simpplot(M.nodes,M.elements(:,1:3),ca), hold(ca,'all')
                        
                        % Plot Boundaries
                        for jj = 1:size(R.curves,2)
                            c = R.curves(jj);
                            plot(ca,c.x(:,1),c.x(:,2),'-'), hold(ca,'all')
                        end
                        
                        % Get Domain Quadrature Points
                        x = G.plate.planar.x(:,1);
                        y = G.plate.planar.x(:,2);
                        
                        % Plot
                        plot(ca,x(:),y(:),'.m')
                        
                        % Get Boundary Quadrature Points
                        x = G.plate.boundary.xq(:,1);
                        y = G.plate.boundary.xq(:,2);
                        t = G.plate.boundary.tq;
                        C = G.plate.boundary.cq;                        
                        
                        % Plot
                        plot(ca,x(:),y(:),'.r')

                        % Plot Curves
                        for jj = 1:size(R.curves,2)
                            c = R.curves(jj);
%                             sort(t(C==jj))
                            u =  R.mappingC(t(C==jj),c.x,c.mapping.param);
                            plot(ca,u(:,1),u(:,2),'o')
%                             [u,t(C==jj)]
                           
%                             sum(t(C==jj)<0)
                        end
                        xlabel(ca,'x (m)','fontunits','normalized'),ylabel(ca,'y (m)','fontunits','normalized')
                        
                    case {'Plate-Domain','Parameter-Domain'}
                        % Compute Points
                        u = linspace(0,1,50);
                        [U,V] = meshgrid(u,u);
                        
                        % Compute Boundary
                        switch names{ii}
                            case 'Plate-Domain'
                                B = R.mesh.plate.planarElem;
                            case 'Parameter-Domain'
                                B = R.mesh.parameter.planarElem;
                        end

                        P = B.param.P;
                        ind = inpolygon(U(:),V(:),P(1:3,1),P(1:3,2));
                        
                        % Compute basis
                        Z = R.getBasis([U(:),V(:)],B.param);
                        Z(~ind,:) = NaN;
                                                
                        for jj = 1:size(Z,2)
                            z = reshape(Z(:,jj),size(U));
                            surf(ca,U,V,z,z*0+jj,'facealpha',0.3),hold(ca,'all')
                        end
                        shading(ca,'interp')
                        
                        % Plot Boundary
                        plot(ca,P([1:3,1],1),P([1:3,1],2),'-ok')
                        plot(ca,P(:,1),P(:,2),'ok')
                        
                        % Plot Quadrature
                        Qp = B.q;
                        plot(ca,Qp(:,1),Qp(:,2),'*')
                        
                        % Rename
                        names{ii} = [names{ii}, ' - Order ',num2str(B.param.order)]; 
                        xlabel(ca,'u_1','fontunits','normalized'),ylabel(ca,'u_2','fontunits','normalized')
                        
                    case {'Plate-Boundary','Parameter-Thickness','Plate-Thickness'}
                        % Compute Boundary
                        switch names{ii}
                            case 'Plate-Boundary'
                                B = R.mesh.plate.boundaryElem;
                            case 'Parameter-Thickness'
                                B = R.mesh.parameter.thicknessElem;
                            case 'Plate-Thickness'
                                B = R.mesh.plate.thicknessElem;
                        end
                        
                        % Get Points
                        u = linspace(B.param.P(1),B.param.P(end),100)';
 
                        % Compute basis
                        Z = R.getBasis(u,B.param);
                        
                        % Plot
                        P = B.param.P;
                        plot(ca,u,Z,'-'), hold(ca,'on')
                        plot(ca,P,0*P,'ok-')
                        Qp = B.q;
                        plot(ca,Qp,0*Qp,'k*')
                        names{ii} = [names{ii}, ' - Order ',num2str(B.param.order)];   
                        xlabel(ca,'u','fontunits','normalized')
                end
                
                % Add Labels
                title(ca,names{ii},'fontunits','normalized'), grid(ca,'on')
                hold(ca,'off')%, axis(ca,[-1 1 -1 1])
                
                % Add UIContext Menu
                set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'),'parent',mP)
                set(ca,'fontunits','normalized')
            end
        end
        function plotInputDisturbances(obj,component)
            % Update component
            obj.component = component;
            I = component.design.inputDisturbances;
            R = component.design.region;
            inputs = component.inputs;
            
            % Get panel
            mP = obj.parentPanel;
            
            % Clear Children
            delete(findobj(mP,'tag','Colorbar'))
            delete(get(mP,'children'))
            
            % Get positions
            N = 15;
            positions = obj.getPositions(N);
            
            % Get Statistics
            t = [I.loads.t;1];
            [loads,matAllow,mParam] = I.getStatistics();
            
            % Get Initial Component information
            DATA = getappdata(1,'DATA');
            G = DATA.GeometryMaker;
            [loadsi,matAllowi,mParami] = G.inputDisturbances.getStatistics();
                        
            % Plot Loads
            L = reshape(inputs.loads,[],5,size(inputs.loads,2));
            M = inputs.matAllow;
            for ii = 1:5
                % Specify Axes
                ca = axes('units','normalized','position',positions{ii});
                V1 = loads.u([1:end,end],ii);
                V2 = loadsi.u([1:end,end],ii);
                V3 = squeeze(L([1:end,end],ii,:));
                
%                 V = [loads.u([1:end,end],ii),loadsi.u([1:end,end],ii),L([1:end,end],ii)]*1e6;
                
                % Plot
                fill([t;flipud(t)],[loadsi.U([1:end,end],ii);flipud(loadsi.L([1:end,end],ii))],[1 0 0],'facealpha',0.1,'EdgeColor',[1 0 0],'parent',ca), grid(ca,'on'),hold(ca,'on')
                plot(ca,t,V1,'r-.','linewidth',1)
                
                
                fill([t;flipud(t)],[loads.U([1:end,end],ii);flipud(loads.L([1:end,end],ii))],[0 0 1],'facealpha',0.15,'EdgeColor',[0 0 1],'parent',ca), grid(ca,'on'),hold(ca,'on')
                plot(ca,t,V2,'b--','linewidth',1)
                plot(ca,t,V3,'m-','linewidth',1.5)
                xlabel(ca,'s'),title(I.loads.names{ii},'fontunits','normalized')
               
%                 % Plot
%                 plot(ca,t,[loads.u(:,ii),loads.U(:,ii),loads.L(:,ii),L(:,ii)],'linewidth',2), grid(ca,'on')
%                 xlabel(ca,'s'),title(I.loads.names{ii},'fontunits','normalized')
% %                 l=legend(ca,'\mu','+2\sigma','-2\sigma','c');
                
                % Add UIContext Menu
                set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'),'parent',mP)
                set(ca,'fontunits','normalized')
            end
            
            % Plot Materials
            UT = matAllow.u;
            DT = matAllow.S;
            Ui = matAllowi.u;
            Di = matAllowi.S;
%             D = [matAllow.u,matAllow.U,matAllow.L,inputs.matAllow];
            for ii = 1:5
                % Specify Axes
                ca = axes('units','normalized','position',positions{ii+5});
                
                % Plot pdf
                % Get Points
                xT = UT(ii) + linspace(-1,1,100)'*DT(ii)*3;
                xi = Ui(ii) + linspace(-1,1,100)'*Di(ii)*3;
                
                % Create pdf
                pT = normpdf(xT,UT(ii),DT(ii));
                pT = [0*xT;flipud(pT)];
                
                pi = normpdf(xi,Ui(ii),Di(ii));
                pi = [0*xi;flipud(pi)];
                
                % Plot
                fill([xT;flipud(xT)],pT,[0 0 1],'facealpha',0.15,'EdgeColor',[0 0 1],'parent',ca), grid(ca,'on'),hold(ca,'on')
                fill([xi;flipud(xi)],pi,[1 0 0],'facealpha',0.15,'EdgeColor',[1 0 0],'parent',ca)  
                plot([M(ii,:);M(ii,:)]',[0 max(pT)*1.1],'m-','linewidth',1.5)
                box(ca,'on')
%                 plot(ca,d,0,'o','markersize',10,'linewidth',2), grid(ca,'on')
                xlabel(ca,I.matAllow.names{ii},'fontunits','normalized'),title(ca,I.matAllow.names{ii},'fontunits','normalized')
                ax = axis;
                axis(ca,[ax(1:3) max(pi)*5])
%                 legend(ca,'\mu','+2\sigma','-2\sigma','c')
                
                % Add UIContext Menu
                set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'),'parent',mP)
                set(ca,'fontunits','normalized')
            end
            
            % Plot load variance
            p = positions{11};
            ca = axes('units','normalized','position',p);
            S = I.loads.S;
            surf(ca,S,'edgecolor','none')
            view(ca,[0,90])
            title(ca,'Loads Covariance')
            colorbar('peer',ca,'location','eastoutside')
            axis(ca,'image')
            
            % Add UIContext Menu
            set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'),'parent',mP)
            set(ca,'fontunits','normalized')
            
            % Plot Manufacturing
            % Specify Axis
            c = R.complexity;
            p = [positions{1+11}(1:2),positions{end-1}(1) - positions{1+11}(1),positions{1+11}(4)];
            p2 = [positions{end-1}(1:2),positions{end}(1:2)+positions{end}(3:4)-positions{end-1}(1:2)];
            ca = axes('units','normalized','position',p);
            h = bar(ca,60*[mParami.U,mParam.U,c.mTime],'w'); hold(ca,'all')
            g=bar(ca,60*[mParami.u,mParam.u,c.mTime]);
            obj.manufacturingTable(p2,mP,c.mTime',c.mLabel)
            box(ca,'on')
            title(ca,'Manufacturing Complexity','fontunits','normalized')
            xlabel(ca,'Manufacturing Step','fontunits','normalized')
            ylabel(ca,'Time (min)','fontunits','normalized'), grid(ca,'on')
            set(g(1),'FaceColor','r');
            set(g(2),'FaceColor','b');
            set(g(3),'FaceColor','c');
            set(h,'edgecolor','k')
            axis(ca,'tight')            
            
            % Add UIContext Menu
            set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'),'parent',mP)
            set(ca,'fontunits','normalized')          
            
        end
        function plotLaminate(obj,component,plotType)
            % Update component
            obj.component = component;
            R = component.design.region;
            DATA = getappdata(1,'DATA');
            G = DATA.GeometryMaker.geometry;
            
            % Get panel
            mP = obj.parentPanel;
            
            % Clear Children
            delete(findobj(mP,'tag','Colorbar'))
            delete(get(mP,'children'))
            
            % Specify Number of Layers to Plot
            N = 8;
            tq = linspace(-1,1,N)';
                                    
            % Get Axis Locations
            positions = obj.getPositions(N+1);
            
            % Get Level Coordinates
            XY = R.laminate.paramz.x;
            
            % Compute Global Basis for centroid
            P = R.mesh.parameter;
            M = R.getBasis([1/3,1/3],P.planarElem.param);
            M = repmat(M,size(P.elements,1),1);
            V = zeros(size(M,1),size(XY,1));
            for ii = 1:size(M,1)
                V(ii,P.elements(ii,:)) = M(ii,:);
            end
            
            % Calculate Coordinates for centroid
            XY = V*XY;
            
            % Compute Thickness
            b = R.laminate.paramz.b;
            Z = b(2) + b(3)*exp(b(1)*R.designVariables.z);
            Zm = b(2) + b(3)*exp(b(1)*V*R.designVariables.z);
            for ii = 1:length(tq)  
                % Compute t coordinate
                % Reference
                tr = tq(ii);
                Nr = R.getBasis(tr,R.mesh.parameter.thicknessElem.param);
                Ar = reshape(permute(Nr,[1 3 2]).*V,size(V,1),[])*R.designVariables.a;
                
                % Max
                ta = tq(ii)*max(Z)./Zm;
                Na = R.getBasis(ta,R.mesh.parameter.thicknessElem.param);
                Aa = reshape(permute(Na,[1 3 2]).*V,size(V,1),[])*R.designVariables.a;
                
                % Specify Axes
                ca = axes('unit','normalized','position',positions{ii});
                
                % Null out things not in domain
                ind = abs(ta) > 1;
                Aa(ind,:) = NaN;
                               
                % Plot
                switch plotType
                    case 'Reference'
                        quiver(ca,XY(:,1),XY(:,2),cos(Ar),sin(Ar),0.7), hold(ca,'all')
                    case 'Actual'
                        quiver(ca,XY(:,1),XY(:,2),cos(Aa),sin(Aa),0.7), hold(ca,'all')
                    case 'Both'
                        quiver(ca,XY(:,1),XY(:,2),cos(Ar),sin(Ar),0.7), hold(ca,'all')
                        quiver(ca,XY(:,1),XY(:,2),cos(Aa),sin(Aa),0.7)
                        
                end
                 
                for jj = 1:size(G.Y,2)
                    y = G.Y{jj};
                    z = G.Z{jj};
                    plot(ca,y(:,1),y(:,2),'-',z(:,1),z(:,2),'ko:')
                end

%                 streamslice(ca,Xt,Yt,cos(Z),sin(Z))
                grid(ca,'on'),axis(ca,'tight')
                
                % Add Labels
                xlabel(ca,'x (m)','fontunits','normalized')
                ylabel(ca,'y (m)','fontunits','normalized')
                title(ca,['Ply Angle at t_n = ',num2str(tq(ii)),' (rad)'],'fontunits','normalized')
                
                % Add UIContext Menu
                set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'),'parent',mP)
                set(ca,'fontunits','normalized')
            end
            
            % Plot Thickness  
            % Create a container object for surface plots
            ca = axes('units','normalized','position',positions{ii+1},'parent',mP);
            
            % Plot
            b = R.laminate.paramz.b;
            XY = R.laminate.paramz.x;
            Z = b(2) + b(3)*exp(b(1)*R.designVariables.z);
            trisurf(P.elements(:,1:3),XY(:,1),XY(:,2),Z,'edgecolor',0.5*[1 1 1])
            
            colorbar('peer',ca,'location','Eastoutside'),hold(ca,'on')
            zm = min(Z);
            for jj = 1:size(G.Y,2)
                y = G.Y{jj};
                z = G.Z{jj};
                plot3(ca,y(:,1),y(:,2),y(:,1)*0+zm,'o-')
                plot3(z(:,1),z(:,2),z(:,1)*0+zm,'k-')
            end
            
            title(ca,'Thickness (m)','fontunits','normalized'), xlabel(ca,'x (m)','fontunits','normalized'), ylabel(ca,'y (m)','fontunits','normalized')
            view(ca,[0 90]), axis(ca,'tight')
            set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'))
            set(ca,'fontunits','normalized')
            
            colormap(linspecer)
            set(ca,'CLim',[0.999*min(Z),max(Z)*1.001])
           
            
        end
                
        function plotDisplacements(obj,component)            
            % Get Reference to Main Panel
            obj.component = component;
            mP = obj.parentPanel;
            
            DATA = getappdata(1,'DATA');
            G = DATA.GeometryMaker.geometry;
            
            % Clear Children
            delete(findobj(mP,'tag','Colorbar'))
            delete(get(mP,'children'))
            
            % Specify labels
            names = {'u_1 (m)','u_2 (m)','p_1 (rad)','p_2 (rad)','w (m)'};
            
            % Get Grid
            positions = obj.getPositions(6);
            
            % Get Fields
            region = component.design.region;
            WS = reshape(region.field.x,[],5);
            
            % Loop through components
            for ii = 1:length(names)
                % Specify Axes
                ca = axes('unit','normalized','position',positions{ii});
                
                % Plot field
                XY = region.mesh.plate.nodes;
                T = region.mesh.plate.elements(:,1:3);
                trisurf(T,XY(:,1),XY(:,2),WS(:,ii),'edgecolor',0.5*[1 1 1],'parent',ca)
                view(ca,[0 90]), hold(ca,'on')
                axis(ca,'tight')
                
                for jj = 1:size(G.Y,2)
                    y = G.Y{jj};
                    z = G.Z{jj};
                    plot(ca,y(:,1),y(:,2),'-',z(:,1),z(:,2),'k-')
                end
                              
                % Add labels
                title(ca,names{ii}),xlabel(ca,'x (m)'),ylabel(ca,'y (m)'), grid(ca,'on')
                hold(ca,'off')
                
                % Add UIContext Menu
                set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'),'parent',mP)
                
                % Add colorbar last
                colorbar('peer',ca,'location','Eastoutside')
            end
            
            % Plot As
            A = region.field.A;
            ca = axes('unit','normalized','position',positions{ii+1});
            spy(A)
            title(ca,'Sparsity Pattern')
            set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'),'parent',mP)
            
        end
        
        function plotStrains(obj,component)            
            % Get Reference to Main Panel
            obj.component = component;
            mP = obj.parentPanel;
            
            % Get geometry
            DATA = getappdata(1,'DATA');
            G = DATA.GeometryMaker.geometry;
            
            % Clear Children
            delete(findobj(mP,'tag','Colorbar'))
            delete(get(mP,'children'))
            
            % Specify labels
            names = {'e_{11}^*','e_{22}^*','2e_{12}^*','e_{11}','e_{22}','2e_{12}','k_{11}','k_{22}','2k_{12}','2y_{13}','2y_{23}'};
            
            % Get Grid
            positions = obj.getPositions(12);
            
            % Get Fields
            region = component.design.region;         

            % Get Strains on Surface
            Gb = region.buildGlobalBasis(0);
            [S,~,~,Sinplane] = computeStrain(region,Gb,[],Gb.plate.planar.x,1);
            WSinplane = reshape(Sinplane,[],3);
            WS = reshape(S,[],8);
            WS1 = [WSinplane,WS];
            
            % Project onto Nodes from Quadrature
            N = Gb.plate.planar.N;
            w = Gb.plate.planar.w;
            WS = ((w.*N)'*N)\((w.*N)'*(WS1));
            
            for ii = 1:length(names)
                % Specify Axes
                ca = axes('unit','normalized','position',positions{ii});
                
                % Plot field
                XY = region.mesh.plate.nodes;
                T = region.mesh.plate.elements(:,1:3);
                trisurf(T,XY(:,1),XY(:,2),WS(:,ii),'edgecolor',0.5*[1 1 1],'parent',ca)
                view(ca,[0 90]), hold(ca,'on')
                axis(ca,'tight')
%                 Gb.plate
%                 plot3(Gb.plate.planar.x(:,1),Gb.plate.planar.x(:,2),WS1(:,ii),'*')
                
                for jj = 1:size(G.Y,2)
                    y = G.Y{jj};
                    z = G.Z{jj};
                    plot(ca,y(:,1),y(:,2),'-',z(:,1),z(:,2),'k-')
                end
                
                % Add labels
                title(ca,names{ii}),xlabel(ca,'x (m)'),ylabel(ca,'y (m)'), grid(ca,'on')
                hold(ca,'off')
                
                % Add UIContext Menu
                set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'),'parent',mP)
                
                % Add colorbar last
                colorbar('peer',ca,'location','Eastoutside')
            end
            
        end
        
        function plotComplexity(obj,component)   
            % Update component
            obj.component = component;
            R = component.design.region;
            DATA = getappdata(1,'DATA');
            G = DATA.GeometryMaker.geometry;
            
            % Get panel
            mP = obj.parentPanel;
            
            % Clear Children
            delete(findobj(mP,'tag','Colorbar'))
            delete(get(mP,'children'))
            
            % Specify Number of Layers to Plot
            N = 8;
            tq = linspace(-1,1,N)';
                                    
            % Get Axis Locations
            positions = obj.getPositions(N+1);
            
            % Get Level Coordinates
            XY = R.laminate.paramz.x;
            
            % Project Complexity Back onto Nodes
            c = R.complexity.Hfd;
            Gb = R.buildGlobalBasis(0);
            N = Gb.angle.volume.N;
            w = Gb.angle.volume.w;
            Cp = ((w.*N)'*N)\((w.*N)'*c);
                        
            % Compute Global Basis at centroids
            P = R.mesh.parameter;
            M = R.getBasis([1/3,1/3],P.planarElem.param);
            M = repmat(M,size(P.elements,1),1);
            V = zeros(size(M,1),size(XY,1));
            for ii = 1:size(M,1)
                V(ii,P.elements(ii,:)) = M(ii,:);
            end
            
            % Compute Global Basis at Nodal Locations
            Vn = speye(size(P.nodes,1),size(P.nodes,1));

            % Calculate Coordinates for centroid
            XY = V*XY;
                        
            % Compute Thickness
            for ii = 1:length(tq)  
                % Compute t coordinate
                % Reference
                t = tq(ii);
                N = R.getBasis(t,R.mesh.parameter.thicknessElem.param);
                C = kron(N,Vn)*Cp;
                A = kron(N,V)*R.designVariables.a;
                                
                % Specify Axes
                ca = axes('unit','normalized','position',positions{ii});
                               
                % Plot
                T = R.mesh.parameter.elements(:,1:3);
                nodes = R.mesh.parameter.nodes;
                trisurf(T,nodes(:,1),nodes(:,2),C*0,C,'edgecolor',0.5*[1 1 1],'parent',ca,'facealpha',0.9)
                view(ca,[0 90]), hold(ca,'all')
                axis(ca,'tight')
                quiver(ca,XY(:,1),XY(:,2),cos(A),sin(A),0.7,'color',0.9*[1 1 1])
                                 
                for jj = 1:size(G.Y,2)
                    y = G.Y{jj};
                    z = G.Z{jj};
                    plot(ca,y(:,1),y(:,2),'-',z(:,1),z(:,2),'ko:')
                end

%                 streamslice(ca,Xt,Yt,cos(Z),sin(Z))
                grid(ca,'on'),axis(ca,'tight')
                
                % Add Labels
                xlabel(ca,'x (m)','fontunits','normalized')
                ylabel(ca,'y (m)','fontunits','normalized')
                title(ca,['Complexity at t_n = ',num2str(tq(ii))],'fontunits','normalized')
                
                % Add UIContext Menu
                set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'),'parent',mP)
                set(ca,'fontunits','normalized')
                
                % Add colorbar last
                colorbar('peer',ca,'location','Eastoutside')
            end
          
        end
        function plotMargins(obj,component)
            % Get Reference to Main Panel
            obj.component = component;
            mP = obj.parentPanel;
            
            % Get geometry
            DATA = getappdata(1,'DATA');
            G = DATA.GeometryMaker.geometry;
            R = component.design.region;
            
            % Clear Children
            delete(findobj(mP,'tag','Colorbar'))
            delete(get(mP,'children'))
            
            % Get Margins Data
            margins = R.margins;
            
            % Specify Number of Layers to Plot
            N = 8;
            tq = linspace(-1,1,N)';
                                    
            % Get Axis Locations
            positions = obj.getPositions(N+1);
                        
            % Project Margins Back onto Nodes
            f = margins.fTL;
            Gb = R.buildGlobalBasis(0);
            N = Gb.plate.volume.N;
            w = Gb.plate.volume.w;
            if size(N,1) ~= size(f,1)
                N = Gb.margin.volume.N;
                w = 1;
            end
            F = ((w.*N)'*N)\((w.*N)'*f);
            max(f)
            
            % Compute Global Basis at Nodal Locations
            Vn = speye(size(Gb.thickness.planar.N,2),size(Gb.thickness.planar.N,2));
            
            for ii = 1:(length(tq)+1)
                % Specify Axes
                ca = axes('unit','normalized','position',positions{ii});
                
                % Calculate field
                if ii < length(tq) + 1
                    % Get Margins on layer
                    % Reference
                    t = tq(ii);
                    N = R.getBasis(t,R.mesh.parameter.thicknessElem.param);
                    FM(:,ii) = kron(N,Vn)*F;   
                else
                    FM(:,ii) = max(FM,[],2);
                end

                % Plot
                T = R.mesh.parameter.elements(:,1:3);
                nodes = R.mesh.parameter.nodes;
                trisurf(T,nodes(:,1),nodes(:,2),FM(:,ii),'edgecolor',0.5*[1 1 1],'parent',ca)
                view(ca,[0 90]), hold(ca,'all')
                axis(ca,'tight')

                for jj = 1:size(G.Y,2)
                    y = G.Y{jj};
                    z = G.Z{jj};
                    plot(ca,y(:,1),y(:,2),'o-',z(:,1),z(:,2),'k-')
                end
                
                % Add labels
                if ii < length(tq) + 1
                    title(ca,['Margin at t_n = ',num2str(tq(ii))],'fontunits','normalized')
                else
                    title(ca,['Max Margin Through Thickness'],'fontunits','normalized')
                end
                xlabel(ca,'x (m)'),ylabel(ca,'y (m)'), grid(ca,'on')
                hold(ca,'off')
                
                % Add UIContext Menu
                set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'),'parent',mP)
                
                % Add colorbar last
                colorbar('peer',ca,'location','Eastoutside')
            end
            
        end
        
        % UI Context menu functions
        function plotmenu = createPlotMenu(obj,pType)
            % Create context menu
            cmenu = uicontextmenu;
            
            switch pType
                case 'mainWindow'
                    % Set tag for cmenu
                    set(cmenu,'tag',['Menu_',num2str(1)])
                    
                    % Create the hierarchy of menu options
                    menuOptions{1} = {'Geometry'
                        'Mesh'
                        'Input Disturbances'
                        'Laminate'
                        'Complexity'
                        'Displacements'
                        'Strains'
                        'Margins'}';
                    menuOptions{2} = {{'View','Sketch','Drag','Delete','Rectangle','Circle','Import'},{'Manufacture and Deploy','Experiment','Wait'},{'Reference','Actual','Both'}};
                    menuOptions{3} = {{'Selection','All'}};
                    
                    % Create the association matrix
                    assocMat = {1,[], [], []
                                2,[], [], []
                                3,[], [], []
                                4,3, [], []
                                5,[], [], []
                                6,[], [], []
                                7, [],[],[]
                                8, [],[],[]};
                                       
            end
            
            % Build menu and submenu
            for i = 1:1:size(assocMat,1)
                % Get starting menu options
                menu1 = uimenu(cmenu, 'label',menuOptions{1}{i});
                
                % Get the 2nd sub menu
                subID2 = assocMat{i,2};
                if ~isempty(subID2)
                    labels2 = menuOptions{2}{subID2};
                    
                    % Loop through second menu
                    for j = 1:1:size(labels2,2)
                        menu2 = uimenu(menu1,'label',labels2{j});
                        
                        % Get the 3rd sub menu
                        subID3 = assocMat{i,3};
                        if ~isempty(subID3)
                            % Determine if there are dependencies to second menu
                            if size(subID3,2) > 1
                                if subID3(j) == 0
                                    labels3 = {};
                                    set(menu2,'callback',@obj.menuButtonCallback)
                                else
                                    labels3 = menuOptions{3}{subID3(j)};
                                end
                            else
                                labels3 = menuOptions{3}{subID3};
                            end
                            
                            % Loop through third sub menu
                            for k = 1:1:size(labels3,2)
                                menu3 = uimenu(menu2,'label',labels3{k});
                                
                                % Get the 4th sub menu
                                subID4 = assocMat{i,4};
                                if ~isempty(subID4)
                                    % Determine if there are dependencies to second menu
                                    if size(subID4,2) > 1
                                        if subID4(j) == 0
                                            labels4 = {};
                                            set(menu3,'callback',@obj.menuButtonCallback)
                                        else
                                            labels4 = menuOptions{4}{subID4(j)};
                                        end
                                    else
                                        labels4 = menuOptions{4}{subID3};
                                    end
                                    
                                    % Loop through 4th sub menu
                                    for l = 1:1:size(labels4,2)
                                        menu4 = uimenu(menu3,'label',labels4{l});
                                        set(menu4,'callback',@obj.menuButtonCallback)
                                    end
                                else
                                    set(menu3,'callback',@obj.menuButtonCallback)
                                    if (i+j+k) == 3
                                        % Initialize menu start and place check mark here
                                        set(menu3,'Checked', 'on')
                                    end
                                end
                            end
                        else
                            set(menu2,'callback',@obj.menuButtonCallback)
                        end
                    end
                else
                    set(menu1,'callback',@obj.menuButtonCallback)
                end
            end
            
            % Return cmenu
            plotmenu = cmenu;
            
        end
        function menuButtonCallback(obj,hObject,eventData)
            % Find chain of handles up to parent
            handleID = [];
            menuLabels = {};
            handleTemp = hObject;
            
            if isfield(eventData,'menuLabels')
                menuLabels = eventData.menuLabels;
            else
                while ~strncmp(get(handleTemp,'tag'),'Menu',4)
                    % Get handles and labels
                    handleID = [handleID,handleTemp];
                    menuLabels = [menuLabels,get(handleTemp,'label')];
                    
                    handleTemp = get(handleTemp,'Parent');
                end
                
                % Get parent menu id
                parentID = str2num(strrep(get(handleTemp,'tag'),'Menu_',''));
                
                % Get Top Level Parent to Child Reference
                menuLabels=fliplr(menuLabels);
            end
                               
            % Get applicaiton data
            DATA = getappdata(1,'DATA');
            
            % Do based on menulabel
            switch menuLabels{1}
                case 'FEM'
                    
                case 'Margins'
                    % Run Margins                              
                case 'Complexity'
                    % Run Complexity
                    
            end
            
            % Update DATA menuID
            DATA.CurrentMenuOption = menuLabels;
            setappdata(1,'DATA',DATA);
            
            % Specify Plot
            plot(obj,obj.component,menuLabels)
                          
        end
        
        % Plotter grid functions
        function positions = getPositions(obj,N)
            % Get Grid
            [M,N] = obj.subplotGrid(N);
            
            % Get Axis Positions
            xTol = 0.04;
            yTol = 0.05;
            xGap = 0.045;
            yGap = 0.07;
            H = 1/M*(1-2*yTol - (M-1)*yGap);
            L = 1/N*(1-2*xTol - (N-1)*xGap);
            counter = 1;
            for ii = 1:M
                for jj = 1:N
                    % Calculate new position
                    x(1) = xTol + (jj-1)*(L+xGap);
                    x(2) = 1 - yTol - H - (ii-1)*(H+yGap);
                    x(3) = L;
                    x(4) = H;
                    positions{counter} = x;
                    counter = counter + 1;
                end
            end
            
        end
        function [rowI,colI] = subplotGrid(obj,num)
            % Function to develop best sized grid for numPoints
            if isprime(num) && num~=2
                num = num+1;
            end
            
            % Factor
            fact1 = [1:1:num];
            fact2 = num./fact1;
            
            % Apply filters
            facFilter = {@(x,y)(y==fix(y))         % Find integer factors
                @(x,y)((x+y)==min(x+y))   % Minimize Perimeter
                @(x,y)(x>=y)};            % Bias towards more rows
            
            % Loop through filters
            for i = 1:1:size(facFilter)
                f = facFilter{i};
                ind = f(fact1,fact2);
                fact1 = fact1(ind);
                fact2 = fact2(ind);
            end
            rowI = fact2;
            colI = fact1;
        end
        
        % Table Functions
        function manufacturingTable(obj,position,parent,t,l)
            % Augment with Total
            total = sum(t)*60;
            l = [l,'Total'];
            
            % Specify Values
            [~,indr] = sort(t,'descend');
            indI = 1:length(t);
            indI(indr) = indI;
            
            values = [l;num2cell([t*60,total]);num2cell([indI,0])]';
            headers = {'Description','Minutes','Rank'};
            
            % Make Table
            mtable = uitable('Units','normalized','Position',position,...
                'Data',values,'tag','ComplexityTable','ColumnName',headers,...
                'ColumnFormat',{'numeric'},'ColumnEditable',false,...
                'Visible','on','RowName',[1:length(t)]','Fontsize',10,'parent',parent);
            set(mtable,'fontunits','normalized')
            
            % Autoresize
            jScroll = findjobj(mtable);
            jTable = jScroll.getViewport.getView;
            jTable.setAutoResizeMode(jTable.AUTO_RESIZE_SUBSEQUENT_COLUMNS);
        end
        
        
    end
    
end