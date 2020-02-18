function X = multiprod3(varargin)
% About: Compute multiarray matrix multiplication A1 * A2 ... * AN of N 3D-Matrix-Arrays
% Higher dimensions are collapsed into 3rd dimension and re-expanded after
% inverse computation
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% ///////////////////////////////////////////////////////////////////////// 
% Get First Matrix
A = varargin{1};

% Get Dimensions
NA = size(A);
indA = fullfact([NA(1:2),prod(NA(3:end))]);

% Compute System Indices
mA = indA(:,[1 2]) + NA([1 2]).*(indA(:,3)-1);

% Build System
AA = sparse(mA(:,1),mA(:,2),A(:));

% Loop through other matrices
X = AA;
N = length(varargin);
for ii = 2:N
    % Get Matrices
    B = varargin{ii};
    
    % Get Dimensions
    NB = size(B);
    indB = fullfact([NB(1:2),prod(NA(3:end))]);
    
    % Compute System Indices
    mB = indB(:,[1 2]) + NB([1 2]).*(indB(:,3)-1);
    
    % Build System
    BB = sparse(mB(:,1),mB(:,2),B(:));
    
    % Compute Operation
    X = X * BB;
    
end

% Return Results
NC = [NA(1) NB(2) NA(3:end)];
Z = repmat(eye(NC(2),NB(2)),prod(NC(3:end)),1);
X = X*Z;
X = permute(reshape(transpose(X),NC([2 1 3:length(NA)])),[2 1 3:length(NA)]);
end