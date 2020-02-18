function derivativeScript()
% About: Script to determine which derivatives need to be evaluated based
% on children and parent variables
% Calculates one stage derivative and hessian
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% /////////////////////////////////////////////////////////////////////////
clc
func = 'Hf';
p = {'Hfd','wfq'};
c = {'A','Z'};

% Form Derivative
s = '';
h = '';

% Specify Type
% typeString = ' %s d%s.d%s.*(d2%s.d%s.d%s).*d%s.d%s';
typeString = ' %s transpose(d%s.d%s)*bsxfun(@times,d2%s.d%s.d%s,d%s.d%s) ... \n           ';

% Compute Fiber Complexity - Calculate New Derivatives 
% Parents ('Hf') Children ('Hfd','wfq')
dHf.dwfq = 1;
dHf.dHfd = 1;

% Compute preliminary quantities
% Parents ('Hfd','wfq') children('A','Z')
dHfd.dA = 1;
dHfd.dZ = 1;
      
dwfq.dA = 1;
dwfq.dZ = 1;

% Compute Fiber Complexity - Calculate New Hessians 
% Parents ('Hf') Children ('Hfd','wfq')
d2Hf.dwfq.dwfq = 0;
d2Hf.dwfq.dHfd = 1;
d2Hf.dHfd.dwfq = d2Hf.dwfq.dHfd;
d2Hf.dHfd.dHfd = 0;

% Compute preliminary quantities
% Parents ('Hfd','wfq') children('A','Z')
d2Hfd.dA.dA = 1;
d2Hfd.dA.dZ = 1;
d2Hfd.dZ.dA = 1;
d2Hfd.dZ.dZ = 1;

d2wfq.dA.dA = 1;
d2wfq.dA.dZ = 1;
d2wfq.dZ.dA = 1;
d2wfq.dZ.dZ = 1;

% 
% % Compute Fiber Complexity - Calculate New Derivatives
% dHf.dyq = 1;
% dHf.dHfd = 1;
% 
% dHfd.ddAs = 1;
% dHfd.dz = 1;
% dHfd.da = 1;
% dHfd.ddZs = 1;
% dyq.ddAs = 0;
% dyq.dz = 0;
% dyq.ddZs = 0;
% 
% ddZs.ddz1 = 1;
% ddZs.ddz2 = 1;
% ddAs.ddA1 = 1;
% ddAs.ddA2 = 1;
% ddAs.ddA3 = 1;
% dyq.dyq = 1;
% dz.dz = 1;
% da.da = 1;
% 
% ddA1.dA = 1;
% ddA2.dA = 1;
% ddA3.dA = 1;
% ddz1.dA = 0;
% ddz2.dA = 0;
% ddz1.dZ = 1;
% ddz2.dZ = 1;
% dz.dZ = 1;
% da.dA = 1;
% dyq.dZ = 1;
% 
% % Compute Fiber Complexity - Calculate New Hessians
% d2Hf.dyq.dyq = 0;
% d2Hf.dyq.dHfd = 1;
% d2Hf.dHfd.dyq = d2Hf.dyq.dHfd;
% d2Hf.dHfd.dHfd = 0;
% 
% d2Hfd.ddAs.ddAs = 0;
% d2Hfd.ddAs.dz = 1;
% d2Hfd.dz.ddAs = d2Hfd.ddAs.dz;
% d2Hfd.ddAs.ddZs = 0;
% d2Hfd.ddZs.ddAs = d2Hfd.ddAs.ddZs;
% d2Hfd.dz.dz = 1;
% d2Hfd.da.da = 1;
% d2Hfd.dz.ddZs = 0;
% d2Hfd.ddZs.ddZs = 0;
% d2Hfd.ddZs.dz = d2Hfd.dz.ddZs;
% 
% d2dAs.ddA1.ddA1 = 2;
% d2dAs.ddA1.ddA2 = 0;
% d2dAs.ddA2.ddA2 = 2;
% d2dAs.ddA2.ddA1 = d2dAs.ddA1.ddA2;
% d2dAs.ddA1.ddA3 = 0;
% d2dAs.ddA3.ddA1 = d2dAs.ddA1.ddA3;
% d2dAs.ddA2.ddA3 = 0;
% d2dAs.ddA3.ddA2 = d2dAs.ddA2.ddA3;
% d2dAs.ddA3.ddA3 = 2;
% 
% d2dZs.ddz1.ddz1 = 2;
% d2dZs.ddz1.ddz2 = 0;
% d2dZs.ddz2.ddz2 = 2;
% d2dZs.ddz2.ddz1 = d2dZs.ddz1.ddz2;
% 
% d2z.dZ.dZ = 1;
% d2dz1.dZ.dZ = 1;
% d2dz2.dZ.dZ = 1;
% d2yq.dZ.dZ = 1;


for kk = 1:length(c)
    % Build string for derivative
    s = [s,sprintf('d%s.d%s =',func,c{kk})];
    for ii = 1:length(p)  
        if ii == 1
            addStr = '';
        else
            addStr = '+';
        end
        s = [s,sprintf(' %s d%s.d%s.*d%s.d%s',addStr,func,p{ii},p{ii},c{kk})];
         
    end
    s = [s,sprintf(';\n')];
    
    % Build String for Hessian
    for ii = 1:length(c)
        h = [h,sprintf('d2%s.d%s.d%s =',func,c{kk},c{ii})];
                
        for jj = 1:length(p)
            for mm = 1:length(p)
                
                if jj == 1 && mm == 1
                    addStr = '';
                else
                    addStr = '+';
                end
                % Do First Term
                v1 = sprintf('d%s.d%s',p{jj},c{kk});
                v2 = sprintf('d%s.d%s',p{mm},c{ii});

                if ~isfield(eval(['d',p{jj}]),['d',c{kk}])
                    eval([v1,'=0;'])
                    
                end
                if ~isfield(eval(['d',p{mm}]),['d',c{ii}])
                    eval([v2,'=0;'])
                end
                
                if any(strcmp(p{jj},c)) && ~strcmp(p{jj},c{kk})
                    eval([v1,'=0;'])
                end
                
                if any(strcmp(p{mm},c)) && ~strcmp(p{mm},c{ii})
                    eval([v2,'=0;'])
                end
                
                if (eval(v1) ~=0) && (eval(v2) ~= 0)
                    h = [h,sprintf(typeString,addStr,p{jj},c{kk},func,p{jj},p{mm},p{mm},c{ii})];
                end
            end
            
        end
        
        % Add Second term
        for jj = 1:length(p)
            v1 = sprintf('d2%s.d%s.d%s',p{jj},c{kk},c{ii});

            % Do First Term
            if exist(['d2',p{jj}],'var') ~= 0
                if isfield(eval(['d2',p{jj}]),['d',c{kk}])
                    if isfield(eval(['d2',p{jj},'.d',c{kk}]),['d',c{ii}])
                        % Check parent child dependence
                        if any(strcmp(p{jj},c)) && (~strcmp(p{jj},c{kk}) || ~strcmp(p{jj},c{ii}))
                            eval([v1,'=0;']);
                        end
                    else
                        eval([v1,'=0;']);
                    end
                else
                    eval([v1,'=0;']);
                end
                
            else
                eval([v1,'=0;']);
            end
            
            if eval(v1) ~= 0
                h = [h,sprintf(' + (d%s.d%s).*d2%s.d%s.d%s',func,p{jj},p{jj},c{kk},c{ii})];
            end
        end
        h = [h,sprintf(';\n')];
        
    end
end
s
strrep(strrep(h,'= +','='),'=;',' = 0;')

end