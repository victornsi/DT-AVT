classdef PolicyVisualizer < handle
% About: Class definition for Policy Plot Visualizer
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% /////////////////////////////////////////////////////////////////////////
    properties
        parentPanel     % Parent panel for visualizer
        policies        % Current poicies
    end
    
    methods
        % Class constructor
        function obj = PolicyVisualizer(parentPanel) 
            obj.parentPanel = parentPanel;
        end       
        
        % Plot Functions
        function plot(obj,policies,plotType)
            % Plot based on requested plot type
            if ~iscell(plotType)
                h = plotType;
            else
                h = plotType{1};
            end
            switch h
                case {'Costs and Constraints'}
                    obj.plotCostsConstraints(policies);
                case 'Basis Functions'
                    obj.plotBasisFunctions(policies,plotType{2});           
            end
            
        end
        function plotCostsConstraints(obj,policies)
            % Update policies
            obj.policies = policies;
            
            % Plot Geometry
            mP = obj.parentPanel;
            
            % Clear Children
            delete(findobj(mP,'tag','Colorbar'))
            delete(get(mP,'children'))
                                                
            % Get position of axis
            Np = size(policies,2);
            M = 4;
            positions = obj.getPositions(M,size(policies,2),[]);
            
            % Loop through each policy
            for ii = 1:Np
                % Plot Costs
                ca = axes('units','normalized','position',positions{ii});
                R = policies(ii).Costs;
                Rs = cumsum(R,2);
                t = 0:size(R,2)-1;
                plot(ca,t,Rs,'.-'), grid(ca,'on')
                xlabel(ca,'t'),title(ca,['Costs | Policy ',num2str(ii),' [',[policies(ii).rules{:},' - ',policies(ii).modes],']'],'fontunits','normalized')
                set(ca,'parent',mP)
                
                % Set Axis Menu
                set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'),'parent',mP)
                
                % Plot Margin Constraints
                ca = axes('units','normalized','position',positions{ii+Np});
                Gs = policies(ii).MConstraints;
                plot(ca,t,Gs,'.-'), grid(ca,'on')
                xlabel(ca,'t'),title(ca,['Margin Constraints | Policy ',num2str(ii),' [',[policies(ii).rules{:},' - ',policies(ii).modes],']'],'fontunits','normalized')
                set(ca,'parent',mP)
                axis(ca,[0 t(end) -1 1])
                
                % Set Axis Menu
                set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'),'parent',mP)
                
                % Plot Margin Constraints
                ca = axes('units','normalized','position',positions{ii+2*Np});
                Gs = policies(ii).SConstraints;
                plot(ca,t,Gs,'.-'), grid(ca,'on')
                xlabel(ca,'t'),title(ca,['Senser Constraints | Policy ',num2str(ii),' [',[policies(ii).rules{:},' - ',policies(ii).modes],']'],'fontunits','normalized')
                set(ca,'parent',mP)
                axis(ca,[0 t(end) -1 1])
                
                % Set Axis Menu
                set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'),'parent',mP)

                % Store Means of R
                Rm(ii,:) = mean(R,1);
                
                % Build Legend Substring
                l{ii} = [[policies(ii).rules{:}],' - ',policies(ii).modes];
                
            end
            
            % Add Sum
            Rm = [Rm,sum(Rm,2)];
            Rm = Rm./min(Rm(1,end));
%             Rm = Rm/min(Rm(:,end));
            
            % Create xlabel
            xL = cellfun(@(x) ['Stage ',num2str(x)],num2cell(t),'uniformoutput',false);
            
            % Plot Mean Comparisons
            p = [positions{(M-1)*Np+1}(1:2),positions{(M-1)*Np+1}(3:4) + positions{(M)*Np}(1:2)-positions{(M-1)*Np+1}(1:2)];
            ca = axes('units','normalized','position',p,'parent',mP);
            bar(ca,Rm'), grid(ca,'on')
            set(ca,'xticklabel',[xL,'Total'])
            title(ca,'Comparison of Policies | Normalized Mean Costs')
            legend(ca,l,'Location','NorthWest')
            colormap(linspecer)
            
            % Set Axis Menu
            set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'),'parent',mP)
           
                        
        end   
        function plotBasisFunctions(obj,policies,pType)
            % Update policies
            obj.policies = policies;
            
            % Plot Geometry
            mP = obj.parentPanel;
            
            % Get Policy Number
            pn = str2num(regexprep(pType,'Policy ',''));
            
            % Clear Children
            delete(findobj(mP,'tag','Colorbar'))
            delete(get(mP,'children'))
            
            % Get position of axis
            T = length(policies(1).rules);
                        
            % Get Test Data
            TEMP = getappdata(1,'TEMP');
            policiesT = TEMP.Policies;            
            
            % Loop through Stages
            l = {'KJ','HJ','TT','PT'};
            M = length(l);
            positions = obj.getPositions(M,T,[]);
            
            for m = 1:M
            for t = 1:T
                % Plot Axis
                ca = axes('units','normalized','position',positions{t+T*(m-1)});
                
                % Plot Training Data for each policy
                for p = pn%1:size(policies,2)
                    % Get Training Basis
                    if ~isempty(policies(p).policy)
                        IJ = policies(p).policy.coefficients(t).(l{m});
                    else
                        IJ = 0;
                    end
                    ind = 1:size(IJ,1);

                    if size(IJ,2) > 0
%                         IJ = IJ(:,:);%1:min(size(IJ,2),policies(p).policy.coefficients(t).n));
                        plot(ca,ind,IJ,'-','color',0.8*[1 1 1]) 
%                         display(sprintf('p = %d | t = %d | N = %d',p,t,size(IJ,2)))
                        hold(ca,'all')
                        
%                         if strcmp(l{m},'TT')
%                             plot(ca,ind,policies(p).policy.coefficients(t).ao,'r-')
%                         end
                        
                    end
%                     
%                     switch l{m}
%                         case 'TT'
%                             [~,indTT] = min(IJ(2700,:))
%                     end
                    
                    % Get Test Basis
                    IJ2 = policiesT(p).data(t).(l{m});
                    ind = 1:size(IJ2,1);
                    plot(ca,ind,IJ2,'-'), grid(ca,'on')   
                    
                    % Plot Closest training point
                    if size(IJ,2) > 0
                        for k = 1:size(IJ2,2)
                            [~,ii] = min(sum((IJ2(:,k)-IJ).^2,1),[],2);
                           plot(ca,ind,IJ(:,ii),'--','color',0.4*[1 1 1])
                        end
                    end
                           
                end
                hold(ca,'off')
                
                % Add titles
                title(ca,['Policy ',num2str(p),' | ', l{m},' | t = ',num2str(t-1)])
                
                % Set Axis Menu
                set(ca,'uicontextmenu',obj.createPlotMenu('mainWindow'),'parent',mP)
                axis(ca,'tight')
                             
            end
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
                    menuOptions{1} = {'Costs and Constraints'
                                      'Basis Functions'}';
                    pStr = cellfun(@(x) ['Policy ',num2str(x)],num2cell(1:size(obj.policies,2)),'uniformoutput',false);
                    menuOptions{2} = {pStr, {'View','Sketch','Drag','Delete','Rectangle','Circle','Import'},{'Manufacture and Deploy','Experiment','Wait'}};
                    menuOptions{3} = {{'Selection','All'}};
                    
                    % Create the association matrix
                    assocMat = {1,[], [], []
                                2,1, [], []};
%                                 3,[], [], []
%                                 4,[], [], []
%                                 5,[], [], []
%                                 6,[], [], []
%                                 7,[], [], []};
                                       
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
                                           
            % Specify Plot
            plot(obj,obj.policies,menuLabels)
                          
        end
        
        % Plotter grid functions
        function positions = getPositions(obj,M,N,g)
            % Get Axis Positions
            if isempty(g)
                xTol = 0.04;
                yTol = 0.05;
                xGap = 0.03;
                yGap = 0.07;
            else
                xTol = g(1);
                yTol = g(2);
                xGap = g(3);
                yGap = g(4);
            end
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