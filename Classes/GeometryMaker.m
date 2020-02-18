classdef GeometryMaker < handle
% About: Class definition for Geometry Maker
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% /////////////////////////////////////////////////////////////////////////
    properties
        parentPanel     % Parent panel for visualizer
        region          % Region
        inputDisturbances % Input Disturbances
        geometry        % Geometry information
        menuOption      % Selected Menu Option
        inputs          % Unknown inputs
        seed            % Geometry Maker seed
    end
    
    methods
        % Class constructor
        function obj = GeometryMaker(parentPanel) 
            obj.parentPanel = parentPanel;
            obj.initializeGeometry;
        end       
        
        % Geometry initialization function
        function initializeGeometry(obj)
            obj.geometry.Y = [];
            obj.geometry.X = [];
            obj.geometry.Z = [];
            obj.geometry.I = [];
            obj.geometry.activeCurve = [];
            obj.inputDisturbances = [];
            obj.region = Region();
            obj.menuOption = 'View';
        end
        
        % Plot Functions
        function plot(obj)
            % Plot based on requested plot type
%             obj.menuOption = 'Sketch'
            switch  obj.menuOption               
                case {'Sketch','Drag','Rectangle','Circle','View',''}
                    obj.plotGeometry();      
                    
                case 'Mesh'
                    obj.plotMesh();
                    
                case 'Input Disturbances'
                    obj.plotInputDisturbances();
                    
                case 'Laminate'
                    obj.plotLaminate();
            end
            
        end
        function plotGeometry(obj)            
            % Get Reference to Main Panel
            mP = obj.parentPanel;
            
            % Get Handles
            ca = findobj(mP,'tag','SketchAxes');
            if isempty(ca)
                % Clear chilren
                delete(findobj(mP,'tag','Colorbar'))
                delete(get(mP,'children'))
                
                % Create an Axes of appropriate shape
                ca = axes('unit','normalized','position',[0.05 0.075 0.45 0.865]);
                
                % Turn on grid
                grid(ca,'on')
                
                % Reset uicontext menu
                set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'),'parent',mP)
                
                % Initialize axis
                ax = [];
                
            else
                % Preserve Axis
                ax = axis(ca);
            end
            % Plot Sketch Data
            G = obj.geometry;
            R = obj.region;
            
            % Plot Control Points
            if ~isempty(G.X)
                x = G.X;
                plot(ca,x(:,1),x(:,2),'k-o'), hold(ca,'on')
                grid(ca,'on')
            end
            % Plot Interpolation Points
            if ~isempty(G.Y)
                for ii = 1:size(G.Y,2)
                    y = G.Y{ii};
                    z = G.Z{ii};
                    plot(ca,y(:,1),y(:,2),'.-',z(:,1),z(:,2),'ko:'), hold(ca,'all')
                end
                
                % Plot Sensor Locations
                if isfield(R.designVariables,'s')
                    xz = R.sensorLocations(R.designVariables.s);
                    plot(ca,xz(:,1),xz(:,2),'m+')
                end

                grid(ca,'on')
            end
            hold(ca,'off')
            
            % Specify Labels
            set(ca,'tag','SketchAxes')
            title(ca,'Sketch Region'), xlabel(ca,'x (m)'),ylabel(ca,'y (m)')
                        
            % Add Button Down Listener
            set(ca,'ButtonDownFcn',@obj.buttonDownFunction)
            
            % Preserve axis limits if present
            if ~isempty(G.I) && strcmp(G.I{1},'import')
                axis(ca,'equal')
            elseif ~isempty(ax)
                axis(ca,ax)
            end
            
            % Plot Maps
            names = 'Configuration';
            positions = [0.54 0.5 0.43 0.44];

            % Find Exisiting Axis
            ca = findobj(mP,'tag',names);
            if isempty(ca)
                ca = axes('position',positions,'parent',mP);
            else
                set(ca,'position',positions,'parent',mP);
            end
            
            % Plot Mesh
            if ~isempty(R.mesh)
                simpplot(R.mesh.triangulation.nodes,R.mesh.triangulation.elements(:,1:3),ca), hold(ca,'all')
            end
                
            % Plot Curves
            t = linspace(0,1,500)';
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
            set(ca,'tag',names)
            title(ca,[names,' Map']),xlabel(ca,'x (m)'),ylabel(ca,'y (m)')
                
            
            % Plot Boundary Conditions
            ca = findobj(mP,'tag','BC');
            if isempty(ca)
                ca = axes('position',[0.54 0.075 0.43 0.33],'parent',mP);
            end
            t = linspace(0,1,500)';
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
                plot(ca,x(di>=4,1),x(di>=4,2),'ks','markersize',20,'linewidth',2)
                plot(ca,x(di<3 & di>0,1),x(di<3 & di>0,2),'ko','markersize',15,'linewidth',2)
                hold(ca,'off')
            end
            title(ca,'Boundary Conditions'),xlabel(ca,'x (m)'),ylabel(ca,'y (m)'), grid(ca,'on')
            set(ca,'tag','BC')
            view(ca,[0 90])
%             axis(ca,'equal')
            
        end
        function plotMesh(obj)
            % Update component
            R = obj.region;
            
            % Get panel
            mP = obj.parentPanel;
            
            % Clear Children
            delete(findobj(mP,'tag','Colorbar'))
            delete(get(mP,'children'))
            fIndex = get(obj.parentPanel,'Parent');
            
            % Plot Nodes and Quadrature Points
            figure(fIndex)
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
                        
%                         % Plot boundary nodes
%                         nds = elements(M.boundary.elements(:,1),1);
%                         plot(ca,nodes(nds,1),nodes(nds,2),'*')
                        
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
                            plot(u(:,1),u(:,2),'o')
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
        function plotInputDisturbances(obj)
            % Get panel
            mP = obj.parentPanel;

            % Clear Children
            delete(findobj(mP,'tag','Colorbar'))
            delete(get(mP,'children'))
            
            % Get positions
            N = 15;
            positions = obj.getPositions(N);
            
            % Get Statistics
            I = obj.inputDisturbances;
            t = [I.loads.t;1];
            [loads,matAllow,mParam] = I.getStatistics();
                        
            % Plot Loads
            L = reshape(obj.inputs.loads,[],5);
            M = obj.inputs.matAllow;
            for ii = 1:5
                % Specify Axes
                ca = axes('units','normalized','position',positions{ii});
                V = [loads.u([1:end,end],ii),L([1:end,end],ii)];            
                
                fill([t;flipud(t)],[loads.U([1:end,end],ii);flipud(loads.L([1:end,end],ii))],[0 0 1],'facealpha',0.15,'EdgeColor',[0 0 1],'parent',ca), grid(ca,'on'),hold(ca,'on')
                plot(ca,t,V(:,1),'b--','linewidth',1)
                
                plot(ca,t,V(:,2),'m-','linewidth',1.5)
                xlabel(ca,'s'),title(I.loads.names{ii},'fontunits','normalized')
                               
                % Add UIContext Menu
                set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'),'parent',mP)
                set(ca,'fontunits','normalized')
            end
            
            % Plot Materials
            UT = matAllow.u;
            DT = matAllow.S;
            for ii = 1:5
                % Specify Axes
                ca = axes('units','normalized','position',positions{ii+5});
                
                % Plot pdf
                % Get Points
                xT = UT(ii) + linspace(-1,1,100)'*DT(ii)*3;
                
                % Create pdf
                pT = normpdf(xT,UT(ii),DT(ii));
                pT = [0*xT;flipud(pT)];
                                
                % Plot
                fill([xT;flipud(xT)],pT,[0 0 1],'facealpha',0.15,'EdgeColor',[0 0 1],'parent',ca), grid(ca,'on'),hold(ca,'on')
                plot([M(ii) M(ii)],[0 max(pT)*1.1],'m-','linewidth',1.5)
                box(ca,'on')
%                 plot(ca,d,0,'o','markersize',10,'linewidth',2), grid(ca,'on')
                xlabel(ca,I.matAllow.names{ii},'fontunits','normalized'),title(ca,I.matAllow.names{ii},'fontunits','normalized')
                ax = axis;
                axis(ca,[ax(1:3) max(pT)*5])
%                 legend(ca,'\mu','+2\sigma','-2\sigma','c')
                
                % Add UIContext Menu
                set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'),'parent',mP)
                set(ca,'fontunits','normalized')
            end
            
            % Plot Manufacturing
            % Specify Axis
            c = obj.region.complexity;
            p = [positions{1+10}(1:2),positions{end-1}(1) - positions{1+10}(1),positions{1+10}(4)];
            p2 = [positions{end-1}(1:2),positions{end}(1:2)+positions{end}(3:4)-positions{end-1}(1:2)];
            ca = axes('units','normalized','position',p);
            h = bar(ca,60*[mParam.U,c.mTime],'w'); hold(ca,'all')
            g=bar(ca,60*[mParam.u,c.mTime]);
            obj.manufacturingTable(p2,mP,c.mTime',c.mLabel)
            box(ca,'on')
            title(ca,'Manufacturing Complexity','fontunits','normalized')
            xlabel(ca,'Manufacturing Step','fontunits','normalized')
            ylabel(ca,'Time (min)','fontunits','normalized'), grid(ca,'on')
            set(g(1),'FaceColor','r');
            set(g(2),'FaceColor','c');
            set(h,'edgecolor','k')
            axis(ca,'tight')            
            
            % Add UIContext Menu
            set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'),'parent',mP)
            set(ca,'fontunits','normalized')
            
        end
        function plotLaminate(obj)
            % Update component
            R = obj.region;
            G = obj.geometry;
            
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
            N = R.getBasis(tq,P.thicknessElem.param);

            % Calculate Coordinates for centroid
            XY = V*XY;
                         
            for ii = 1:length(tq)
                % Compute Angle
                A = kron(N(ii,:),V)*R.designVariables.a;
                 
                % Specify Axes
                ca = axes('unit','normalized','position',positions{ii});
                
                % Plot                
                quiver(ca,XY(:,1),XY(:,2),cos(A),sin(A)),hold(ca','on')
%                 streamslice(ca,XY(:,1),XY(:,2),cos(A),sin(A)),hold(ca','on')
                grid(ca,'on'),axis(ca,'tight')
                
                for jj = 1:size(G.Y,2)
                    y = G.Y{jj};
                    z = G.Z{jj};
                    plot(ca,y(:,1),y(:,2),'-',z(:,1),z(:,2),'ko-')
                end
                
                % Add Labels
                xlabel(ca,'x (m)','fontunits','normalized')
                ylabel(ca,'y (m)','fontunits','normalized')
                title(ca,['Ply Angle at t_n = ',num2str(tq(ii)), ' (rad)'],'fontunits','normalized')
                
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
        
        % BC table functions
        function launchBCTable(obj)
            % Get Application Data
            R = obj.region;
            
            % Launch UITable
            Nr = 20;
            headers = {'Filter','Load [N_n (N/m), N_t (N/m), M_n (N), M_t (N), Q (N/m)]','Displacement [u_n (m), u_t (m), p_n (rad), p_t (rad), w (m)]'};
            colFormat = {'char','char','char'};
            colEdit = [true true true];
            values = repmat({'','',''},Nr,1);
            if ~isempty(R.boundaryConditions)
                % Use existing data
                for ii = 1:size(R.boundaryConditions,2)
                    b = R.boundaryConditions(ii);
                    values(ii,:) = {b.filter,b.loads,b.displacements};
                end
            end
            figure(2)
            set(2,'name','Boundary Conditions')
            set(2,'position',[161         342        1221         420])
            mtable = uitable('Parent',2,'Units','normalized','Position',[0.2 0.02 0.79 0.97],...
                'Data',values,'tag','BCTable','ColumnName',headers,...
                'ColumnFormat',colFormat,'ColumnEditable',colEdit,...
                'Visible','on','RowName',[1:Nr]','Fontsize',12);
            
            % Make Table to show available points
            iFill = find(strcmp({R.curves(:).type},'Fill'));
            [x,~,~,t] = R.computeTN([],iFill);
            values = [t,x];
            mtable2 = uitable('Parent',2,'Units','normalized','Position',[0.01 0.02 0.19 0.97],...
                'Data',values,'tag','BCTable2','ColumnName',{'t','x (m)','y (m)'},...
                'ColumnFormat',{'numeric'},'ColumnEditable',false,...
                'Visible','on','RowName',[1:size(x,1)]','Fontsize',12);
            
            % Autoresize
            jScroll = findjobj(mtable);
            jTable = jScroll.getViewport.getView;
            jTable.setAutoResizeMode(jTable.AUTO_RESIZE_SUBSEQUENT_COLUMNS);
            
            jScroll = findjobj(mtable2);
            jTable = jScroll.getViewport.getView;
            jTable.setAutoResizeMode(jTable.AUTO_RESIZE_SUBSEQUENT_COLUMNS);
            
            % Set Callback
            set(mtable,'CellEditCallback',@obj.tableEditCallback);
            set(mtable,'CellSelectionCallback',@obj.tableEditCallback)
            set(mtable,'UIContextMenu',obj.createPlotMenu('BCWindow'))
        end
        function tableEditCallback(obj,hObject,eventData)
            % Get Table Data
            values = get(hObject,'Data');
            
            % Store Indices of Current EventData
            set(hObject,'UserData',eventData.Indices)
            
            if min(sum(strcmp(values,''),2))>1
                return
            end
            
            % Get Functions
            counter = 1;
            for ii = 1:size(values,1)
                % Store Functions
                if sum(strcmp(values(ii,:),'')) > 1
                    break
                end
                
                % Add BC's
                d = values{ii,3};
                l = values{ii,2};
                f = values{ii,1};
                if strcmp(f,'')
                    f = '0 <= t & t <=1';
                end
                BC(counter).loads = l;
                BC(counter).filter = f;
                BC(counter).displacements = d;
                counter = counter + 1;
                
            end
            
            % Store Info
            obj.region.boundaryConditions = BC;
            
        end
          
        % UI Context menu functions
        function plotmenu = createPlotMenu(obj,pType)
            % Create context menu
            cmenu = uicontextmenu;
            
            switch pType
                case 'mainWindow'
                    % Set tag for cmenu
                    set(cmenu,'tag',['Menu_',num2str(3)])
                    
                    % Create the hierarchy of menu options
                    menuOptions{1} = {'Geometry'
                                      'Boundary Conditions'
                                      'Mesh'
                                      'Input Disturbances'
                                      'Laminate'}';
                    menuOptions{2} = {{'View','Sketch','Drag','Delete','Rectangle','Circle','Import'},{'View','Generate'}};
                    menuOptions{3} = {{'Selection','All'}};
                    
                    % Create the association matrix
                    assocMat = {1, 1, [0 0 0 1 0 0 0], []
                                2,[], [], []
                                3,2, [], []
                                4,[], [], []
                                5,[], [], []};
                            
                case 'BCWindow'
                    % Set tag for cmenu
                    set(cmenu,'tag',['Menu_',pType])
                    
                    % Create the hierarchy of menu options
                    menuOptions{1} = {'Displacement'
                        'Loading'}';
                    menuOptions{2} = {{'Pin-Pin-Pin  [l,l,o,o,l,o,o,o,o,o]'
                        'Free-Pin-Pin [o,l,o,o,l,o,o,o,o,o]'
                        'Pin-Free-Pin [l,o,o,o,l,o,o,o,o,o]'
                        'Free-Free-Pin [o,o,o,o,l,o,o,o,o,o]'
                        'Hard-Clamp [l,l,l,l,l,o,o,o,o,o]'}',...
                        {'Normal Compression [-l,o,o,o,o]'
                        'Normal Tension [l,o,o,o,o]'}'};
                    % Create the association matrix
                    assocMat = {1, 1, [], []
                        2, 2, [], []};
                    
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
            
            % Get Axis Handle
            ca = findobj('tag','SketchAxes');
            
            % Update based on menu
            switch menuLabels{1}
                case 'Geometry'
                    % Do Based on Geometry
                    obj.menuOption = menuLabels{end};
                    
                    % Delete Based on Case
                    switch menuLabels{end}
                        case 'Import'
                            % Load Data
                            [FileName,PathName,~] = uigetfile({'*.mat'},'Choose a MAT File',pwd);
                            FilePath = [PathName,FileName];
                            display(['Opened: ',FilePath])
                            load(FilePath)
                            
                            % Inherit geometry properties
                            Rbc = obj.region.boundaryConditions;
                            obj.region = Region();
                            obj.geometry.Z = G.Y;
                            obj.geometry.Y = G.Y;
                            obj.geometry.I = G.I;
                            obj.region.boundaryConditions = Rbc;
                            
                            obj.region.update(obj.geometry)
                            obj.menuOption = 'View';
                              
                        case 'Selection'
                            % Turn on brush command
                            brush on
                            brushObj = brush(gcf);
                            set(brushObj,'Color',[0 215 0]/255)
                            
                            % Wait until user makes selection and then turn off brush command
                            choices = uigetpref(...
                                'mygraphics',...
                                'stuff',...
                                'Brushing Options',...
                                {'Select Objects for Brushing/Filtering'
                                ''},...
                                {'Remove Selected','Cancel'},...
                                'DefaultButton','Cancel');
                            brush off
                            
                            % Remove Selected (if any)
                            if strcmpi(choices,'Remove Selected')
                                % Get Selected Data Points
                                % Loop through each axis child
                                c = get(ca,'children');
                                X = [];
                                for ii = 1:length(c)
                                    % Get Points
                                    b = logical(get(c(ii),'BrushData'));
                                    x = [get(c(ii),'XData');get(c(ii),'YDATA')]';
                                    X = [X;x(b,:)];
                                end
                                
                                % Now Compare with every curve
                                for ii = 1:size(obj.geometry.Y,2)
                                    % Compute Distance
                                    Y = obj.geometry.Y{ii};
                                    r = bsxfun(@plus,Y(:,1),-X(:,1)').^2 + bsxfun(@plus,Y(:,2),-X(:,2)').^2;
                                    d(ii) = min(r(:));
                                end
                                indR = find(d<0.0001^2);
                                
                                % Remove Entities
                                obj.geometry.Z(indR) = [];
                                obj.geometry.Y(indR) = [];
                                obj.geometry.I(indR) = [];
                                obj.geometry.Region = [];
                            end
                            obj.region.update(obj.geometry)
                            obj.menuOption = 'View';
                        case 'All'
                            % Clear All Geometry but preserve BC
                            Rbc = obj.region.boundaryConditions;
                            obj.initializeGeometry();
                            obj.region.boundaryConditions = Rbc;
                            
                            % Remove All Children from axis
                            delete(get(ca,'children'))
                    end
                                
            case 'Boundary Conditions'
                obj.launchBCTable();
                obj.menuOption = 'View';
                
                case {'Displacement','Loading'}
                    % Update Table Data
                    mtable = findobj('tag','BCTable');
                    ind = get(mtable,'UserData');
                    
                    if strcmpi(menuLabels{1},'Displacement')
                        ind = [ind(1),3];
                    else
                        ind = [ind(1),2];
                    end
                    
                    % Update Values
                    values = get(mtable,'Data');
                    
                    % Parse MenuLabel
                    str = regexprep(menuLabels{end},'.+(?=[)','');
                    
                    % Update Values
                    values{ind(1),ind(2)} = str;
                    set(mtable,'Data',values)
                    
                    % Update table Values
                    eData.Indices = ind;
                    obj.tableEditCallback(mtable,eData)

                case 'Mesh'  
                    % Set whether to rerun mesher
                    switch menuLabels{end}
                        case 'View'
                            flg = 0;
                        case 'Generate'
                            flg = 1;
                    end
                    
                    tic
                    obj.region.getMesh(flg)
                    obj.seed = obj.region.designVariables;
                    t=toc;
                    display(['Time to Compute Mesh: ',num2str(t),' sec'])
                    obj.menuOption = menuLabels{end-1};
                    
                case 'Input Disturbances'
                    % Initialize if non-existent
                    if 1%isempty(obj.inputDisturbances)
                        % Get Materials
                        load Materials
                        M = Materials;
                        m = [M.F11t;
                            M.F11c;
                            M.F22t;
                            M.F22c;
                            M.F12];
                        
                        G = obj.region.buildGlobalBasis(1);
                        
                        % Get Loads
                        tbq = G.plate.boundary.tq;
                        ibq = G.plate.boundary.cq;
                        iFill = find(strcmp({obj.region.curves(:).type},'Fill'));
                        L = obj.region.calculateBoundaryConditions(tbq(ibq==iFill));
                        LC = L;
                        
                        % Get Manufacturing Parameters
                        computeComplexity(obj.region,G,0,'');
                        A = obj.region.complexity.A;
                        P = obj.region.complexity.mParam;
                        
                        % Initialize Input Disturbances
                        obj.inputDisturbances = InputDisturbances(tbq(ibq==iFill),1.5*LC(:),0.9*m,P*1.2,A);
                        
                        % Initialize true inputs
                        obj.inputs.loads = L(:);
                        obj.inputs.matAllow = m;
                        obj.inputs.manufacturingParam = P; 
                            
                    end
                    computeComplexity(obj.region,G,obj.inputs.manufacturingParam,'');
                    obj.menuOption = menuLabels{end};
                    
                case 'Laminate'
                    obj.menuOption = menuLabels{end};
                              
            end
            
            % Update plot
            if ~isempty(hObject)
                obj.plot()
            end
        end
        
        % Sketch primitives functions
        function xs = generateSpline(obj,x,m)
            % Close Loop
            x = x([1:end,1],:);
            
            % Specify parameters
            d = 4;
            tb = linspace(0,1,(d+1)*size(x,1)+d);
            t = linspace(0,1,size(x,1));
            
            % Generate Spline Boundary
            switch m
                case 'spline'
                    s = 0.5*([-1 1]*x(1:2,:) + [-1 1]*x(end-1:end,:));
                    xs = spline(t,[s;x;s]',tb);
                    xs = xs';
                    
                case 'Rectangle'
                    d = 4;
                    tb = linspace(0,1,(d+1)*size(x,1)+d);
                    xs = interp1(t,x,tb,'linear');
                    
                case {'Circle'}
                    d = 2;
                    tb = linspace(0,1,(d+1)*size(x,1)+d);
                    xs = interp1(t,x,tb,'pchip');
                    
                case 'import'
                    xs = x(1:end-1,:);
            end
            
        end
        function [xs,xc] = generateRectangle(obj,x,xp)
            % Generate Rectangle from x and xp
            % Compute mean
            xm = mean([x;xp],1);
            
            % Specify Parameterization
            d = 4;
            t = linspace(0,1,(d+1)*9+d)';
            
            % Get Major and Minor Radius
            L = abs(xp-xm);
            
            % Create Ellipse
            r = [cos(2*pi*t),sin(2*pi*t)];
            
            % Saturate to Edges
            r = bsxfun(@times,r,1./max(abs(r),[],2));
            
            % Scale and Shift
            r = bsxfun(@plus,r*diag(L),xm);
            xs = r;
            
            % Return Corner Points
            xc = [1 0
                1 1
                0 1
                -1 1
                -1 0
                -1 -1
                0 -1
                1 -1];
            %   xc = [1 -1;1 1;-1 1;-1 -1];
            xc = bsxfun(@plus,xc*diag(L),xm);
            
        end
        function [xs,xc] = generateCircle(obj,x,xp)
            % Generate Circle from x to xp
            % Compute mean
            xm = x;
            
            % Specify Parameterization
            d = 2;
            t = linspace(0,1,(d+1)*9+d)';
            tc = linspace(0,1,9)';
            
            % Get Radius
            R = sqrt(sum((x-xp).^2));
            
            % Create Ellipse
            r = [cos(2*pi*t),sin(2*pi*t)];
            
            % Scale and Shift
            r = bsxfun(@plus,r*R,xm);
            xs = r;
            
            % Return Corner Points
            rc = [cos(2*pi*tc),sin(2*pi*tc)];
            r = bsxfun(@plus,rc*R,xm);
            xc = r(1:end-1,:);
        end

        % Mouse Click Events
        % Cursor Callbacks
        function buttonDownFunction(obj,hobj,event_obj)
            % Set Up Function For Axes (via gcf)
            set(gcf,'WindowButtonUpFcn',@obj.WindowButtonUpFunction)
            
            % Do based on menu selection
            if strcmp(obj.menuOption,'Drag') || ...
                    strcmp(obj.menuOption,'Rectangle') || ...
                    strcmp(obj.menuOption,'Circle')
                % Add Motion Button
                set(gcf,'WindowButtonMotionFcn',@obj.WindowButtonMotionFunction)
                
                % Do Additional Steps for Rectangle and Circle
                if strcmp(obj.menuOption,'Rectangle') || strcmp(obj.menuOption,'Circle')
                    % Initialize Data
                    ca = findobj(hobj,'tag','SketchAxes');
                    cpx = get(ca,'CurrentPoint');
                    cpx = cpx(1,1:2);
                    obj.geometry.X = [];
                    obj.geometry.X(1,:) = cpx;
                    obj.geometry.activeCurve = size(obj.geometry.Y,2)+1;
                end
                
            end

        end
        function WindowButtonUpFunction(obj,hobj,event_obj)
            % Remove Motion and Self
            set(hobj,'WindowButtonUpFcn','')
            set(hobj,'WindowButtonMotionFcn','')
            
            % Get Handle to Plotter
            ca = findobj(hobj,'tag','SketchAxes');
            
            % Get Selection Type for Figure
            selectionType = get(hobj,'SelectionType');
            
            % Get Current Point
            cpx = get(ca,'CurrentPoint');
            cpx = cpx(1,1:2);
            
            % Add Point to Data Struct
            % Only add point if left click
            obj.geometry.activeCurve = [];
            menuOption = obj.menuOption;
            if strcmp(selectionType,'normal') && strcmp(obj.menuOption,'Sketch')
                % Add to X
                obj.geometry.X = [obj.geometry.X;cpx];
                
            elseif (strcmp(selectionType,'extend') || strcmp(selectionType,'alt')) && (strcmp(obj.menuOption,'Sketch'))
                % Clear menu option
                menuOption = '';
                
                % Interpolate Curve Data
                Z = obj.geometry.X;
                Y = obj.generateSpline(Z,'spline');
                
                % Add New Curve (Interpolated and Control Points)
                obj.geometry.Y = [obj.geometry.Y,{Y}];
                obj.geometry.Z = [obj.geometry.Z,{Z}];
                obj.geometry.I = [obj.geometry.I,{'spline'}];
                
                % Clear Point Storage
                obj.geometry.X = [];
                
                % Update Region
                obj.region.update(obj.geometry)
                
            elseif (strcmp(selectionType,'extend') || strcmp(selectionType,'alt')) && (strcmp(obj.menuOption,'Drag'))
                % Clear menu option
                menuOption = '';
                
            elseif strcmp(obj.menuOption,'Rectangle') || strcmp(obj.menuOption,'Circle')
                % Deactivate Rectangle or Circle
                menuOption = '';
                
                % Clear Point Storage
                obj.geometry.X = [];
                
                % Update Region
                obj.region.update(obj.geometry)
                
            elseif strcmp(selectionType,'normal') && strcmp(obj.menuOption,'Drag')
                % Clear menu option
                menuOption = '';
                obj.region.update(obj.geometry)
                
            end
                        
            % Replot
            obj.plot()
            
            % Update Menu Options
            obj.menuOption = menuOption;
            
        end
        function WindowButtonMotionFunction(obj,hobj,event_obj)
            % Get Handle to Plotter
            ca = findobj(hobj,'tag','SketchAxes');
            
            % Get Current Point
            cpx = get(ca,'CurrentPoint');
            cpx = cpx(1,1:2);
            
            % Get Active Curve
            ac = obj.geometry.activeCurve;
            
            if strcmp(obj.menuOption,'Drag')
                if ~isempty(ac)
                    % Get Control Points of Current Curve
                    Z = obj.geometry.Z{ac};
                    
                    % Update Control Point
                    [~,ind] = min(sum(bsxfun(@plus,Z,-cpx).^2,2));
                    obj.geometry.Z{ac}(ind,:) = cpx;
                    
                    % Update Current Set of Points
                    obj.geometry.X = obj.geometry.Z{ac};
                    
                    % Interpolate Curve Data
                    Z = obj.geometry.X;
                    Y = obj.generateSpline(Z,obj.geometry.I{ac});
                    
                    % Update Interpolant
                    obj.geometry.Y{ac} = Y;
                    
                    % Update Region and Plot
                    obj.region.update(obj.geometry)
                    obj.plot()
                    
                    % Clear DATA.X
                    obj.geometry.X = [];
                    
                else
                    % Find Closest
                    r = [];
                    
                    for ii = size(obj.geometry.Z,2):-1:1
                        % Compute Minimum Distance
                        Z = obj.geometry.Z{ii};
                        r(ii) = min(sum(bsxfun(@plus,Z,-cpx).^2,2));
                    end
                    
                    % Store Current Active Curve
                    [~,obj.geometry.activeCurve] = min(r);
                    
                end
            elseif strcmp(obj.menuOption,'Rectangle')
                % Draw Appropriate Size Rectangle
                xb = obj.geometry.X(1,:);
                xp = cpx;
                
                % Get Points for Rectangle
                [xs,xc] = obj.generateRectangle(xb,xp);
                
                % Store in Z and Y
                obj.geometry.Z{ac} = xc;
                obj.geometry.Y{ac} = xs;
                obj.geometry.I{ac} = 'Rectangle';
                
                % Update Plot
                obj.plot()
                
            elseif strcmp(obj.menuOption,'Circle')
                % Draw Appropriate Size Circle
                xb = obj.geometry.X(1,:);
                xp = cpx;
                
                % Get Points for Circle
                [xs,xc] = obj.generateCircle(xb,xp);
                
                % Store in Z and Y
                obj.geometry.Z{ac} = xc;
                obj.geometry.Y{ac} = xs;
                obj.geometry.I{ac} = 'Circle';
                
                % Update Plot
                obj.plot()
                
            end

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
               
        
    end
    
end