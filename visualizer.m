function visualizer()
% About: Visualizer for Digital Thread
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% ///////////////////////////////////////////////////////////////////////// 
addpath(genpath('Classes'))
addpath(genpath('Functions'))
addpath(genpath('Materials'))
addpath(genpath('Examples'))
addpath(genpath('Toolbox'))
addpath(genpath('Debug'))

global figObj
figObj = figure(1);

% Setup main panel
set(0, 'FixedWidthFontName', 'Courier New');
set(0,'DefaultAxesFontSize',8)
delete(get(figObj,'children'))
set(figObj,'name','Digital Thread Visualization Tool v3')
mP = uipanel(figObj,'units','normalized','position',[0.13 0.01 0.865 0.98],...
                  'tag','mainPanel','uicontextmenu',createOptionsMenu);
             
% Launch Tree
buildTree
              
% Initialize Data Struct
[DATA] = initializeData(1);

% Store Application Data
setappdata(1,'DATA',DATA);

% Make Tree
updateTree()

end

% -------------------------------------------------------------------------
% Initialize
% -------------------------------------------------------------------------
function [DATA,TEMP] = initializeData(flg)
% Initialize a Digital Thread Trajectory
for jj = 1:2
for ii = 1:1:1
    TEMP.Policies(1).Trajectories(jj).states(ii) = DigitalThread(ii-1,ii-1,Design(Region(),InputDisturbances(1,1,1,1,1)));
end
end

% Specify Policies
p = {{'E','D','D'}
     {'E','D','D'}
     {'E','D','D'}
     {'D','E','D'}
     {'D','E','D'}
     {'D','E','D'}};
m = {'G','P','S','G','P','S'};
d = [0 1.0 1.0 0 1.0 1.0];
ind = [1 2 3 4 5 6];

p = p(ind);
m = m(ind);
d = d(ind);

for ii = 1:size(p,1)
    Policies(ii).rules = p{ii};
    Policies(ii).discounts = d(ii);
    Policies(ii).modes = m{ii};
    Policies(ii).Costs = [];
    Policies(ii).MConstraints = [];
    Policies(ii).SConstraints = [];
    Policies(ii).ESU = [];
    TEMP.Policies(ii) = TEMP.Policies(1);
end
DATA.Policies = Policies;

% Get Reference to Main Panel
if flg
    mP = findobj('tag','mainPanel');
    DATA.PanelID = mP;
    DATA.CurrentMenuOption = {'Geometry'};
end
setappdata(1,'TEMP',TEMP);
end

function D0 = initializeDigitalThread()
DATA = getappdata(1,'DATA');
G = DATA.GeometryMaker;
design = Design(G.region,G.inputDisturbances);
D0 = DigitalThread(0,0,design);
D0.getFeatureVector();
end

function DATA = initializePolicies()
% Get application data
DATA = getappdata(1,'DATA');
G = DATA.GeometryMaker;
% inputs = G.inputs;

% Initialize Inputs
eventData.menuLabels = {'Input Disturbances'};
G.menuButtonCallback([],eventData)
setappdata(1,'XY',[]);

% Initialize DATA
PoliciesOld = DATA.Policies;
Dp = initializeData(0);
TEMP = getappdata(1,'TEMP');
DATA.Policies = Dp.Policies;
TEMP.inputs = G.inputs;

% Initialize Design and Digital Thread
D0 = initializeDigitalThread();

% Loop through each policy
P = DATA.Policies;

for p = 1:size(P,2)
    % Get Action
    A = P(p).rules;
    T = TEMP.Policies(p).Trajectories;
    Y = P(p).discounts;
    M = P(p).modes;

    % Loop through number of trajectories
    for t = 1:size(T,2)
        % Initialize beginning of Trajectory
        D(1) = D0.copy;
        
        % Modify Control and Get Design Variables
        dV = D(1).design.region.designVariables;
        u = [dV.a;dV.z;dV.s(:);dV.p];
         
        % Loop through each action
        for a = 1:size(A,2)
            % Get Action
            ub = A{a};
            
            % Get Input
            iD = D(1).design.inputDisturbances;
            inputs.loads = iD.loads.u;
            inputs.matAllow = iD.matAllow.u;
            inputs.manufacturingParam = iD.manufacturingParam.u;
       
            % Transition Digital Thread
            D(a+1) = transitionModel(D(a),ub,u,inputs);

        end
        
        % Update Trajectory
        TEMP.Policies(p).Trajectories(t).states = D;
                
        % Create Policy
        if t == 1
           DATA.Policies(p).policy = Policy(D,A,Y,M);     
        end
        
    end
end
setappdata(1,'TEMP',TEMP)
% DATA.Policies([1 2 3 4 6]) = PoliciesOld([1 2 3 4 6]);% Uncomment to reuse previous
disp('All Policies Successfully Initialized')

end
% -------------------------------------------------------------------------
% Operations
% -------------------------------------------------------------------------
function DATA = optimizePolicies()
% Get application data
DATA = getappdata(1,'DATA');

% Initialize Design and Digital Thread
D0 = initializeDigitalThread();

% Do for each policy
tic
for ii = 1:size(DATA.Policies,2)
    P = DATA.Policies(ii);
    bellmanBackup(D0.copy,P.policy,ii); 
end
toc
disp('Policies Optimized')
% saveFunction(DATA,'DATA.mat');

end

function DATA = evaluatePolicies()
% Evaluate Policies
% Get application data
DATA = getappdata(1,'DATA');
G = DATA.GeometryMaker;
inputs = G.inputs;

% Initialize Design and Digital Thread
D0 = initializeDigitalThread();
N = D0.design.inputDisturbances.samples.N;
D0.design.FEMInfo.N = 30;

% Loop through each policy
P = DATA.Policies;
T = 2;
for p = 1:size(P,2)
    % Get Action      
    A = P(p).rules;
    Py = P(p).policy;
%     ind = 10;
%     inputs.loads = DATA.GeometryMaker.inputDisturbances.samples.L(:,ind);
%     inputs.matAllow = DATA.GeometryMaker.inputDisturbances.samples.A(:,ind);
%     inputs.manufacturingParam = DATA.GeometryMaker.inputDisturbances.samples.P(:,ind);

    p
    % Initialize reward and constraint matrix
    R = [];
    G = [];
    H = [];
    E = [];
    bData = [];
    % Loop through number of trajectory samples
    for t = 1:T
        % Initialize beginning of Trajectory
        D(1) = D0.copy();
                
        % Loop through each action
        r = [];
        g = [];
        h = [];
        flg = 1;
        for a = 1:size(A,2)            
            % Get Action
            U = Py.computeControl(D(a),[],1,0);
            
            % Collect Stage Reward
            r(a) = stageCosts(D(a),U,1:N,'',inputs);
            r(a) = r(a) + U.Y;
                     
            % Collect Stage Constraint
            M = Py.coefficients(a).n;
            [c,~,margins] = stageConstraints(D(a),U);
            g(a) = c(1);
            h(a) = c(2);
            s(a) = U.p;

            F = r(a);

            % Display Costs
            disp(['Cost: ',sprintf('(%d,%d): ',p,a),sprintf('R = %6.4f   V = %6.4f   Z = %6.4f   Y = %6.4f   N = %6d   p = %6.4f   y = %6.4f',F,exp(U.V),exp(U.Z),U.Y,M,U.p,U.y)])
      
            % Transition Digital Thread   
            D(a).design.region.margins = margins;
            D(a+1) = transitionModel(D(a),U.ub,U.u,inputs);
            
            % Store TEMP Data
            bData(a).AT(:,t) = U.a;
            bData(a).VT(:,t) = exp(U.V);
            bData(a).KJ(:,t) = U.kj;
            bData(a).HJ(:,t) = U.hj;
            bData(a).PT(:,t) = U.phit;
            bData(a).TT(:,t) = U.thet;
            
        end
        % Update Trajectory
        if t < 3
            Trajectory(t).states = D; 
        end
        R(t,:) = r;
        G(t,:) = g;
        H(t,:) = h;
        E(t,:) = s;
    end
    
    % Update DATA Struct
    DATA.Policies(p).Costs = R;
    DATA.Policies(p).ESU = E;
    DATA.Policies(p).MConstraints = G;
    DATA.Policies(p).SConstraints = H;
    DATA.Policies(p).rules = A;
    
    % Update TEMP Struct
    TEMP.Policies(p).data = bData;
    TEMP.Policies(p).Trajectories = Trajectory;
end
TEMP.inputs = inputs;
disp('All Policies Successfully Evaluated')
setappdata(1,'TEMP',TEMP)
figure(1)
end

% -------------------------------------------------------------------------
% Tree Functions
% -------------------------------------------------------------------------
% Build tree structure
function buildTree
global figObj
% Remove exisiting tree
try
     delete(findobj('tag','dataTree'))
catch 
end

% Create UITREE
% uitreenode(~,value,name,iconpath,isleaf)
Root = uitreenode('v0','Policies','Policies',fullfile(matlabroot,'toolbox','matlab','icons','webicon.gif'),false);

% Place uitree in window
[mtree, container] = uitree('v0', 'Root',Root, 'Parent',figObj); % parent is ignored
set(container, 'units','normalized','Parent',figObj,'position',[0.005 0.01 .12 .98]);
set(container,'tag','dataTree')
set(mtree, 'NodeSelectedCallback',@treeCallback); 
set(mtree, 'multipleSelectionEnabled',true);
set(container,'UserData',mtree)
end

% Update Tree
function updateTree()
% Get Application Data
DATA = getappdata(1,'DATA');
TEMP = getappdata(1,'TEMP');

% Remove existing elements from tree
jTree = get(findobj('tag','dataTree'),'UserData');
Root = get(jTree,'Root');
Root.removeAllChildren()

% Create nodes for each Policy, Trajectory, and Digital Thread
P = DATA.Policies;
colorVec = @(x,y)(([0 0 1]-[1 0 0])*x/y+[1 0 0])*255;

% Specify icon path
iconPath1 = fullfile(matlabroot,'/toolbox/matlab/icons/bookicon.gif');
iconPath2 = fullfile(matlabroot,'/toolbox/matlab/icons/tool_data_cursor.gif');

% For each policy
for p = 1:1:size(TEMP.Policies,2)
    % Create policy node
    ind = p;
    name = [' Policy ',num2str(p)];
    policyNode = uitreenode('v0',{ind,name},['<html><font color="',sprintf('rgb(%g,%g,%g)',colorVec(p,size(P,2))),'">',name],iconPath1,false);
    
    % For each Trajectory
    T = TEMP.Policies(p).Trajectories;
    for tt = 1:1:size(T,2)
        % Create trajectory node
        ind = 1;
        name = [' Trajectory ',num2str(tt)];
        trajectoryNode = uitreenode('v0',{ind,name},['<html><font color="',sprintf('rgb(%g,%g,%g)',colorVec(tt,size(T,2))),'">',name],iconPath2,false);
        
        % For each Digital Thread
        DT = T(tt).states;
        for d = 1:1:size(DT,2)
            colInd = d;
            t = DT(d).stage;
            dT = DT(d);
            trajectoryNode = addNode([' Digital Thread at t =  ',num2str(t)],d,sprintf('rgb(%g,%g,%g)',colorVec(colInd,size(DT,2))),dT,trajectoryNode);
        end
        policyNode.add(trajectoryNode)
    end
    Root.add(policyNode)
end


% % Update Tree
% for ii = 1:1:size(T,2)
%     colInd = ii;
%     t = T(ii).stage;
%     DT = T(ii);
%     addNode([' Digital Thread at t =  ',num2str(t)],ii,sprintf('rgb(%g,%g,%g)',colorVec(colInd)),DT);
% end
jTree.reloadNode(Root)

end

% Add Nodes
function TJ = addNode(name,ind,col,DT,TJ)
% Specify icon path
iconPath = fullfile(matlabroot,'/toolbox/matlab/icons/greencircleicon.gif');
iconPath2 = fullfile(matlabroot,'/toolbox/matlab/icons/pageicon.gif');

% Add nodes for new data Set
mainNode = uitreenode('v0',{ind,name},['<html><font color="',col,'">',name],iconPath,false);

% Add leafs for Components in Operation and Design
No = size(DT.operation,2);
for ii = 1:No
    str = ['<html><font color="',col,'"> Operation ',num2str(ii)];
    leafNode = uitreenode('v0',ii,str,iconPath2,true);
    mainNode.add(leafNode)
end
str = ['<html><font color="',col,'"> Design'];
leafNode = uitreenode('v0',ii,str,iconPath2,true);
mainNode.add(leafNode)

% Add Node and reload
TJ.add(mainNode)

end

% Callback when tree structure is clicked **
function treeCallback(tree, value)
% Get handles of all selected nodes
N = tree.getSelectedNodes();
selected = {};
for ii = 1:1:size(N,1)
    if N(ii).getLevel == 4
        % Build vector information of selected nodes
        % Get Parent Name
        DTParent = N(ii).getParent;
        DTName = regexprep(char(DTParent.getName),'.+>.','');
        trajParent = DTParent.getParent;
        trajName = regexprep(char(trajParent.getName),'.+>.','');
        policyParent = trajParent.getParent;
        policyName = regexprep(char(policyParent.getName),'.+>.','');
        
        % Get Current Name
        currentName = regexprep(char(N(ii).getName),'.+>.','');
        selected = [selected,{policyName;trajName;DTName;currentName}];
        
    elseif N(ii).getLevel == 0
        selected = {N(ii).getName;''};
    else
        return
    end
end

% Store Selected
DATA = getappdata(1,'DATA');
DATA.Selected = selected;
setappdata(1,'DATA',DATA);

% Update Plot 
updatePanel()
end

% Menu when tree is right clicked 
function optionsMenu = createOptionsMenu()
% Create context menu
cmenu = uicontextmenu;

% Set tag for cmenu
set(cmenu,'tag',['Menu_',num2str(2)])

% Create the hierarchy of menu options
menuOptions{1} = {'Generate Geometry'
                  'Policies'
                  'Load'
                  'Save'}';
menuOptions{2} = {{'Initialize','Optimize','Evaluate'}};
menuOptions{3} = {{'Selection','All'}};

% Create the association matrix
assocMat = {1,[], [], []
            2,1, [], []
            3,[], [], []
            4,[], [], []};
        
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
                        set(menu2,'callback',@optionsMenuCallback)
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
                                set(menu3,'callback',@optionsMenuCallback)
                            else
                                labels4 = menuOptions{4}{subID4(j)};
                            end
                        else
                            labels4 = menuOptions{4}{subID3};
                        end
                        
                        % Loop through 4th sub menu
                        for l = 1:1:size(labels4,2)
                            menu4 = uimenu(menu3,'label',labels4{l});
                            set(menu4,'callback',@optionsMenuCallback)
                        end
                    else
                        set(menu3,'callback',@optionsMenuCallback)
                        if (i+j+k) == 3
                            % Initialize menu start and place check mark here
                            set(menu3,'Checked', 'on')
                        end
                    end
                end
            else
                set(menu2,'callback',@optionsMenuCallback)
            end
        end
    else
        set(menu1,'callback',@optionsMenuCallback)
    end
end

% Return cmenu
optionsMenu = cmenu;
end

% Callback when tree is right clicked
function optionsMenuCallback(hObject,eventData)
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

% Do based on menuselection
DATA = getappdata(1,'DATA');
switch menuLabels{1}
    case 'Generate Geometry'
        delete(get(DATA.PanelID,'children'))
        % Launch Geometry Maker
        if ~isfield(DATA,'GeometryMaker')
            DATA.GeometryMaker = GeometryMaker(DATA.PanelID); 
        end
        DATA.GeometryMaker.plot()
        
    case 'Save'
        % Prompt user to specify path
        [FileName,PathName,~] = uiputfile({'*.mat'},'Save Mat file as',[pwd]);
        FilePath = [PathName,FileName];

        % Save Data
        saveFunction(DATA,FilePath);
                
    case 'Load'
        % Load Data
        [FileName,PathName,~] = uigetfile({'*.mat'},'Choose a MAT File',pwd);
        FilePath = [PathName,FileName];
        display(['Opened: ',FilePath])
        load(FilePath)
        
        % Get Reference to Main Panel again (doesn't exist anymore)
        mP = findobj('tag','mainPanel');
        DATA.PanelID = mP;
        DATA.GeometryMaker.parentPanel = mP;
        disp('DATA loaded')
        
    case 'Policies'
        % Do based on submenu
        switch menuLabels{end}
            case 'Initialize'
                tic
                DATA = initializePolicies();
                toc
            case 'Optimize'
                DATA = optimizePolicies();
            case 'Evaluate'
                DATA = evaluatePolicies();
                
        end
end
setappdata(1,'DATA',DATA);
updateTree()
end
function saveFunction(DATA,FilePath)
% Save Data
if isfield(DATA,'policies')
    DATA = rmfield(DATA,'policies');
end
if ~isempty(FilePath)
    save(FilePath,'DATA','-v7.3')
    disp('DATA saved')
end

end

% -------------------------------------------------------------------------
% Panel Update Functions
% -------------------------------------------------------------------------
function updatePanel()
% Get application data
DATA = getappdata(1,'DATA');
TEMP = getappdata(1,'TEMP');
selected = DATA.Selected;

% Determine number of selected
Ns = size(selected,2);

% Launch different functions depending on selected
if Ns == 1 && strcmp(selected{1},'Policies')
    % Create a Panel within the main panel
    p = getPositions2(1,[0.01 0.01 0.01 0.01]);
    mP = uipanel(DATA.PanelID,'units','normalized','position',p{1},...
                  'tag','subPanel','uicontextmenu',createOptionsMenu);
              
    % Plot Policy
    PV = PolicyVisualizer(mP);
    PV.plot(DATA.Policies,'Costs and Constraints')
    
    return
end

% Create appropriate number of panels
p = getPositions2(Ns,[0.01 0.01 0.01 0.01]);

% Delete Current Children
delete(get(DATA.PanelID,'children'))

% Create a Component Visualizer for each component
for ii = 1:Ns
    % Get selection
    s = selected(:,ii);

    % Get Component 
    indP = str2num(strrep(s{1},'Policy ',''));
    indT = str2num(strrep(s{2},'Trajectory ',''));
    indDT = str2num(strrep(s{3},'Digital Thread at t = ',''));
    indComp = strrep(s{4},'Operation','');
    if strcmp(indComp,'Design')
        % Create inputs
        inputs = TEMP.inputs;
        
        % Create a dummy component to house design
        component = Component(TEMP.Policies(indP).Trajectories(indT).states(indDT+1).design,inputs);
                
    else
        component = TEMP.Policies(indP).Trajectories(indT).states(indDT+1).operation(str2num(indComp));
    end
        
    % Create a panel within the main panel
    mP = uipanel(DATA.PanelID,'units','normalized','position',p{ii},...
                  'tag','subPanel','uicontextmenu',createOptionsMenu);
              
    % Plot Component
    CV = ComponentVisualizer(mP);
    CV.plot(component,DATA.CurrentMenuOption)

end

end

% -------------------------------------------------------------------------
% Plotter Grid Functions
% -------------------------------------------------------------------------
function positions = getPositions(M,N,g)
% Get Grid
% [M,N] = subplotGrid(N);

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
function positions = getPositions2(N,g)
% Get Grid
[N,M] = subplotGrid(N);

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
function [rowI,colI] = subplotGrid(num)
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
    @(x,y)(x>=y)};            % Bias towards more columns

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