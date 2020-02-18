function genGeometry()
% Specify file names
fileName = 'Examples/RibWeb/airfoil.dat';

% Specify chord length
co = 1;%45.4545*0.0254; % in

% Specify Hole Centers and Radius
rh = [0.02 0.035 0.04 0.035 0.02];
dh = linspace(0.2,0.65,5);
lh = 0.0175;

% Import Airfoil data
Y = importData(fileName);

% Truncate Rear and front
yTrun1 = Y(:,1)>0.69;
yTrun2 = Y(:,1)<0.15;
Y(yTrun1,1) = 0.7;
Y(yTrun2,1) = 0.15;



% Extract Corner Points
a = pi/180*20;
ind = false(size(Y,1),1);
for ii = 2:(size(Y,1)-1)
    % Get current point
    p = Y(ii,:);
    
    % Get next and last
    n = Y(ii+1,:) - p;
    l = -Y(ii-1,:) + p;
    
    % Compute angles
    nA = atan2(n(2),n(1));
    lA = atan2(l(2),l(1));
    
    % Check for removal
    if wrapToPi(abs(nA-lA)) > a && ind(ii-1)==false
        ind(ii) = true;
        abs(nA-lA)
    else
        ind(ii) = false;
    end

end

% Compute Distances
ind = find(ind)-1;
t = cumsum(sqrt(sum((Y(2:end,:)-Y(1:end-1,:)).^2,2)));
tc = [t(ind)];
ti = t/max(t);

% Create a lower interpolation
q = 0.02;
T = linspace(0,max(t),floor(max(t)/q));

% Determine Sections between Each Corner Point
TC = [[0;tc(1:end-1)],[tc]];
L = [];
s = [3 0 3 0 3];
for ii = 1:size(TC,1)
    % Determine Count
    n = sum(T > TC(ii,1) & T < TC(ii,2));
    n = n+s(ii);
    l = linspace(TC(ii,1),TC(ii,2),n);
    L = [L,l(1:n)];

end
T = L;
% % Add corners
% for ii = 1:length(tc)
%     % Check if there is a similiar point
%     ind = find(abs(T-tc(ii)) < q);
%     ind = ind(1);
%     if ~isempty(ind)
%         % Replace
%         T(ind) = tc(ii);
%     else
%         % Add
%         T = [T,tc(ii)];
%     end
%      
% end
T = unique(T);
T = T/max(t);
Y = interp1([0;ti],Y,T);
Y(end,:) = Y(1,:);


% Store
G.Y{1} = Y*co;
G.I{1} = 'import';

% Add Holes
phi = linspace(0,1,15)'*2*pi;
u = [cos(phi),sin(phi)];

for ii = 1:length(rh)
    Y = bsxfun(@plus,[dh(ii),lh],rh(ii)*u);
    G.Y{ii+1} = Y*co;    
    G.I{ii+1} = 'import';
end

for ii = 1:size(G.Y,2)
   plot(G.Y{ii}(:,1),G.Y{ii}(:,2),'o-'), hold all 
end

hold off
axis equal
save(['Examples/RibWeb/G.mat'],'G')

end

function Y = importData(fileName)
% Open File
fileID = fopen(fileName,'r');
data = textscan(fileID,'%f %f','HeaderLines',1);

% Store
Y = [data{1},data{2}];

end