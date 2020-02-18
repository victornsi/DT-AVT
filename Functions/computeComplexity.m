function [r,dr,d2r] = computeComplexity(region,G,Pinput,derivativeFlg)
% About: Cost Model
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% /////////////////////////////////////////////////////////////////////////
% Get Parameter Values
[a,da,d2a] = region.mappingA(G.angle.volume);
[z,dz,d2z,d3z] = region.mappingZ(G.thickness.planar);

% Replicate z and dz for thickness
Ntq = size(region.mesh.parameter.thicknessElem.q,1);
Mq = size(z,1);
z = repmat(z,Ntq,1);
dz.d1 = repmat(dz.d1,Ntq,1);
dz.d2 = repmat(dz.d2,Ntq,1);
d2z.d1.dZ = repmat(d2z.d1.dZ,Ntq,1);
d2z.d2.dZ = repmat(d2z.d2.dZ,Ntq,1);
dz.dZ = repmat(dz.dZ,Ntq,1);

% Get Values
GA = G.angle.volume;
w = GA.w;

% Get Jacobian operations
dJ = GA.dJ;
dJ.dz1 = permute(dJ.dz1,[3 1 2]);
dJ.dz2 = permute(dJ.dz2,[3 1 2]);
dJ.dz = permute(dJ.dz,[3 1 2]);
Jinv = permute(GA.Jinv,[3 1 2]);
Ji = permute(Jinv,[1 3 2]);
dJinv = GA.dJinv;
d2Jinv = GA.d2Jinv;

% Compute Derivatives in Spatial Coordinates
dA1 = Jinv(:,1,1).*(da.d1) + Jinv(:,2,1).*(da.d2) + Jinv(:,3,1).*(da.d3);
dA2 = Jinv(:,1,2).*(da.d1) + Jinv(:,2,2).*(da.d2) + Jinv(:,3,2).*(da.d3);

dz1 = Jinv(:,1,1).*(dz.d1) + Jinv(:,2,1).*(dz.d2);
dz2 = Jinv(:,1,2).*(dz.d1) + Jinv(:,2,2).*(dz.d2); 

% Compute Gradients
dAs = dA1.^2 + dA2.^2;
dZs = dz1.^2 + dz2.^2;

% Compute Fiber Complexity ------------------------------------
s = 50;
q = 0.002;
k = 0;
l = 500;
m = -1;
Hfd = q*dAs.*z.^m + k*dAs + s*dZs./z + l;

% Compute Complexity
Hf = transpose(w)*Hfd;
% num2cell(transpose(w)*[q*dAs.*z.^m k*dAs s*dZs./z l+0*z])

% Compute Manufacturing Complexity
dS = 1e-7;
Sa = sum(G.thickness.planar.w);
Sb = sum(G.plate.boundary.wq);
Sv = sum(w);
[t,l,Mf,P]=computeManufacturingComplexity(Sa,Sb,Sv,Pinput);
mT = sum(mean(t,2));

% Store Results
region.complexity.Hf = Hf;
region.complexity.Hfd = Hfd;
region.complexity.x = GA.x;

sM = 1; % Factor for Manufacturing Complexity
r = region.complexity.Hf + sM*mT;
dr = [zeros(size(region.designVariables.a))
      zeros(size(region.designVariables.z))
      zeros(size(region.designVariables.s(:)))
      zeros(size(region.designVariables.p))];
d2r = zeros(size(dr,1));

% Update region
region.complexity.r = r;
region.complexity.mTime = mean(t,2);
region.complexity.mLabel = l;
region.complexity.A = Mf.A;
region.complexity.mParam = P;

% Compute Gradient and Hessian
if strcmpi(derivativeFlg,'derivative')
    % Jacobian nightmare derivatives
    Sp = @(x) spdiags(x,0,size(x,1),size(x,1));
    S = @(x) sum(reshape(x,Mq,[]),2);
    R = @(x) x(1:Mq,:);
    dJinv.dZ = @(ii,jj) Sp(dJinv.dz1(:,ii,jj))*d2z.d1.dZ + Sp(dJinv.dz2(:,ii,jj))*d2z.d2.dZ + Sp(dJinv.dz(:,ii,jj))*dz.dZ;
    
    d2Jinv.dZ.dZ = @(ii,jj,v) transpose(Sp(S(d2Jinv.dz1.dz1(:,ii,jj).*v))*R(d2z.d1.dZ) + Sp(S(d2Jinv.dz2.dz1(:,ii,jj).*v))*R(d2z.d2.dZ) + Sp(S(d2Jinv.dz.dz1(:,ii,jj).*v))*R(dz.dZ))*R(d2z.d1.dZ) + ...
                              transpose(Sp(S(d2Jinv.dz1.dz2(:,ii,jj).*v))*R(d2z.d1.dZ) + Sp(S(d2Jinv.dz2.dz2(:,ii,jj).*v))*R(d2z.d2.dZ) + Sp(S(d2Jinv.dz.dz2(:,ii,jj).*v))*R(dz.dZ))*R(d2z.d2.dZ) + ...
                              transpose(Sp(S(d2Jinv.dz1.dz(:,ii,jj).*v))*R(d2z.d1.dZ) + Sp(S(d2Jinv.dz2.dz(:,ii,jj).*v))*R(d2z.d2.dZ) + Sp(S(d2Jinv.dz.dz(:,ii,jj).*v))*R(dz.dZ))*R(dz.dZ) + ...
                              d3z.d1.dZ.dZ(S(dJinv.dz1(:,ii,jj).*v)) + ...
                              d3z.d2.dZ.dZ(S(dJinv.dz2(:,ii,jj).*v)) + ...
                              d2z.dZ.dZ(S(dJinv.dz(:,ii,jj).*v));
    
    % Derivatives wrt A
    ddA1.dA = Sp(Jinv(:,1,1))*(d2a.d1.dA) + Sp(Jinv(:,2,1))*(d2a.d2.dA) +  Sp(Jinv(:,3,1))*(d2a.d3.dA);
    ddA2.dA = Sp(Jinv(:,1,2))*(d2a.d1.dA) + Sp(Jinv(:,2,2))*(d2a.d2.dA) +  Sp(Jinv(:,3,2))*(d2a.d3.dA);
    
    % Derivatives wrt Z
    ddA1.dZ = Sp(da.d1)*dJinv.dZ(1,1) + Sp(da.d2)*dJinv.dZ(2,1) + Sp(da.d3)*dJinv.dZ(3,1);
    ddA2.dZ = Sp(da.d1)*dJinv.dZ(1,2) + Sp(da.d2)*dJinv.dZ(2,2) + Sp(da.d3)*dJinv.dZ(3,2);
    
    ddz1.dZ = Sp(Jinv(:,1,1))*(d2z.d1.dZ) + Sp(Jinv(:,2,1))*(d2z.d2.dZ) + Sp(dJinv.dZ(1,1))*(d2z.d1.dZ) + Sp(dJinv.dZ(2,1))*(d2z.d2.dZ);
    ddz2.dZ = Sp(Jinv(:,1,2))*(d2z.d1.dZ) + Sp(Jinv(:,2,2))*(d2z.d2.dZ) + Sp(dJinv.dZ(1,2))*(d2z.d1.dZ) + Sp(dJinv.dZ(2,2))*(d2z.d2.dZ);
    
    % Second Derivatives
    d2dA1.dA.dZ = @(v) transpose(Sp(v)*dJinv.dZ(1,1))*d2a.d1.dA + transpose(Sp(v)*dJinv.dZ(2,1))*d2a.d2.dA + transpose(Sp(v)*dJinv.dZ(3,1))*d2a.d3.dA;
    d2dA2.dA.dZ = @(v) transpose(Sp(v)*dJinv.dZ(1,2))*d2a.d1.dA + transpose(Sp(v)*dJinv.dZ(2,2))*d2a.d2.dA + transpose(Sp(v)*dJinv.dZ(3,2))*d2a.d3.dA;
    d2dA1.dZ.dZ = @(v) d2Jinv.dZ.dZ(1,1,(v.*da.d1)) + d2Jinv.dZ.dZ(2,1,(v.*da.d2)) + d2Jinv.dZ.dZ(3,1,(v.*da.d3));
    d2dA2.dZ.dZ = @(v) d2Jinv.dZ.dZ(1,2,(v.*da.d1)) + d2Jinv.dZ.dZ(2,2,(v.*da.d2)) + d2Jinv.dZ.dZ(3,2,(v.*da.d3));
    
    d2dz1.dZ.dZ = @(v) d3z.d1.dZ.dZ(S(v.*Jinv(:,1,1))) +  d3z.d2.dZ.dZ(S((v.*Jinv(:,2,1)))) + ...
                       transpose(Sp(v)*dJinv.dZ(1,1))*d2z.d1.dZ + transpose(Sp(v)*dJinv.dZ(2,1))*d2z.d2.dZ + ...
                       d2Jinv.dZ.dZ(1,1,v.*dz.d1) + d2Jinv.dZ.dZ(2,1,v.*dz.d2);
    
    d2dz2.dZ.dZ = @(v) d3z.d1.dZ.dZ(S(v.*Jinv(:,1,2))) +  d3z.d2.dZ.dZ(S((v.*Jinv(:,2,2)))) + ...
                       transpose(Sp(v)*dJinv.dZ(1,2))*d2z.d1.dZ + transpose(Sp(v)*dJinv.dZ(2,2))*d2z.d2.dZ + ...
                       d2Jinv.dZ.dZ(1,2,v.*dz.d1) + d2Jinv.dZ.dZ(2,2,v.*dz.d2);
    
    % Parents ('Hf') Children ('Hfd','w')
    dHf.dw = Hfd;
    dHf.dHfd = w;
                   
    % Parents ('Hfd','w') children('A','Z')
    dHfd.dA = Sp((q.*z.^m + k).*(2*dA1))*ddA1.dA + Sp((q.*z.^m + k).*(2*dA2))*ddA2.dA;
    dHfd.dZ =  Sp(q.*m.*z.^(m-1).*(dAs))*dz.dZ + Sp((q.*z.^m + k).*(2*dA1))*ddA1.dZ + Sp((q.*z.^m + k).*(2*dA2))*ddA2.dZ + ...
               Sp(s./z.*(2*dz1))*ddz1.dZ + Sp(s./z.*(2*dz2))*ddz2.dZ - Sp(s./z.^2.*(dZs))*dz.dZ;
    dw.dA = sparse(size(da.dA,1),size(da.dA,2));
    dw.dZ = Sp(w.*sum(Ji.*dJ.dz1,[2 3]))*d2z.d1.dZ + Sp(w.*sum(Ji.*dJ.dz2,[2 3]))*d2z.d1.dZ + Sp(w.*sum(Ji.*dJ.dz,[2 3]))*dz.dZ;

    % Reconstruct Derivative
    dHf.dA =  transpose(dHf.dHfd)*dHfd.dA + transpose(dHf.dw)*dw.dA;
    dHf.dZ =  transpose(dHf.dHfd)*dHfd.dZ + transpose(dHf.dw)*dw.dZ;
    
    % Manufacturing Complexity Derivative
    dSa.dA = 0;
    dSb.dA = 0;
    dSv.dA = 0;
    dSa.dZ = 0;
    dSb.dZ = 0;
    dSv.dZ = transpose(sum(dw.dZ,1));
    
    mTSa = 0;%sum(mean(computeManufacturingComplexity(Sa+dS,Sp,Sv,Pinput),2));
    mTSb = 0;%sum(mean(computeManufacturingComplexity(Sa,Sp+dS,Sv,Pinput),2));
    mTSvp = sum(mean(computeManufacturingComplexity(Sa,Sb,Sv+dS,Pinput),2));
    mTSvm = sum(mean(computeManufacturingComplexity(Sa,Sb,Sv-dS,Pinput),2));
    
    dmT.dSa = 0;%imag(mTSa-mT)./imag(dS);
    dmT.dSb = 0;%imag(mTSp-mT)./imag(dS);
    dmT.dSv = (mTSvp-mT)./(dS);
    d2mT.dSv.dSv = (mTSvp - 2*mT + mTSvm)./dS.^2;
    
    dmT.dA = dmT.dSa.*dSa.dA + dmT.dSb.*dSb.dA + dmT.dSv.*dSv.dA;
    dmT.dZ = dmT.dSa.*dSa.dZ + dmT.dSb.*dSb.dZ + dmT.dSv.*dSv.dZ;
    
    % Compute Gradient
    dr = [transpose(dHf.dA) + sM*dmT.dA;
          transpose(dHf.dZ) + sM*dmT.dZ;
          zeros(size(region.designVariables.s(:)))
          zeros(size(region.designVariables.p))];
    
    % Compute Fiber Complexity - Calculate New Hessians 
    % Parents ('Hf') Children ('Hfd','w')
    d2Hf.dw.dw = 0;
    d2Hf.dw.dHfd = 1;
    d2Hf.dHfd.dw = d2Hf.dw.dHfd;
    d2Hf.dHfd.dHfd = 0;
    
    % Parents ('Hfd','w') children('A','Z')
    d2Hfd.dA.dA = transpose(ddA1.dA)*(Sp(2*(q.*z.^m + k).*dHf.dHfd)*ddA1.dA) + transpose(ddA2.dA)*(Sp(2*(q.*z.^m + k).*dHf.dHfd)*ddA2.dA);
    d2Hfd.dA.dZ = + transpose(Sp(2*dA1)*ddA1.dA + Sp(2*dA2)*ddA2.dA)*(Sp((m.*q.*z.^(m-1).*dHf.dHfd))*dz.dZ) ...
                  + transpose(Sp(2.*dHf.dHfd)*ddA1.dA)*(Sp(q.*z.^m + k)*ddA1.dZ) + transpose(Sp(2*dHf.dHfd)*ddA2.dA)*(Sp(q.*z.^m + k)*ddA2.dZ) ...
                  + transpose(d2dA1.dA.dZ(2*dA1.*(q.*z.^m + k).*dHf.dHfd)) + transpose(d2dA2.dA.dZ(2*dA2.*(q.*z.^m + k).*dHf.dHfd));
    d2Hfd.dZ.dA = d2Hfd.dA.dZ;
    d2Hfd.dZ.dZ = + transpose(dz.dZ)*(Sp(m*(m-1)*q.*z.^(m-2).*(dAs.*dHf.dHfd))*dz.dZ) ...
                  + transpose(Sp(2*dA1)*ddA1.dZ + Sp(2*dA2)*ddA2.dZ)*(Sp(m.*q.*z.^(m-1).*dHf.dHfd)*dz.dZ) ...
                  + d2z.dZ.dZ(S(dHf.dHfd.*m.*q.*z.^(m-1).*(dAs)))...
                  + transpose(dz.dZ)*(Sp(m.*q.*z.^(m-1).*dHf.dHfd)*(Sp(2*dA1)*ddA1.dZ + Sp(2*dA2)*ddA2.dZ))...
                  + transpose(Sp(2*dHf.dHfd.*(q.*z.^m + k))*ddA1.dZ)*(ddA1.dZ) + transpose(Sp(2*dHf.dHfd.*(q.*z.^m + k))*ddA2.dZ)*(ddA2.dZ)...
                  + d2dA1.dZ.dZ(dHf.dHfd.*(q.*z.^m + k).*2.*dA1) + d2dA2.dZ.dZ(dHf.dHfd.*(q.*z.^m + k).*2.*dA2)...
                  + transpose(dz.dZ)*(Sp(2*s./z.^3.*(dZs.*dHf.dHfd))*dz.dZ) ...
                  - transpose(Sp(2*dz1)*ddz1.dZ + Sp(2*dz2)*ddz2.dZ)*(Sp(s./z.^2.*dHf.dHfd)*dz.dZ) ...
                  - d2z.dZ.dZ(S(dHf.dHfd.*s./z.^2.*(dZs)))...
                  - transpose(dz.dZ)*(Sp(s./z.^2.*dHf.dHfd)*(Sp(2*dz1)*ddz1.dZ + Sp(2*dz2)*ddz2.dZ))...             
                  + transpose(Sp(2*s./z.*dHf.dHfd)*ddz1.dZ)*ddz1.dZ + transpose(Sp(2*s./z.*dHf.dHfd)*ddz2.dZ)*ddz2.dZ...
                  + d2dz1.dZ.dZ(2*s./z.*dHf.dHfd.*dz1) + d2dz2.dZ.dZ(2*s./z.*dHf.dHfd.*dz2); 
      
    d2w.dA.dA = 0;
    d2w.dA.dZ = 0;
    d2w.dZ.dA = d2w.dA.dZ;
    
    d2w.dZ.dZ = transpose(Sp(dHf.dw)*dw.dZ)*(Sp(sum(Ji.*dJ.dz1,[2 3]))*d2z.d1.dZ + Sp(sum(Ji.*dJ.dz2,[2 3]))*d2z.d1.dZ + Sp(sum(Ji.*dJ.dz,[2 3]))*dz.dZ)...
              + transpose(Sp(dJ.dz1(:,3,1))*dJinv.dZ(1,3))*(Sp(dHf.dw.*w)*d2z.d1.dZ) + transpose(Sp(dJ.dz2(:,3,2))*dJinv.dZ(2,3))*(Sp(dHf.dw.*w)*d2z.d2.dZ) + transpose(Sp(dJ.dz(:,3,3))*dJinv.dZ(3,3))*(Sp(dHf.dw.*w)*dz.dZ)...
              + d3z.d1.dZ.dZ(S(dHf.dw.*w.*sum(Ji.*dJ.dz1,[2 3]))) + d3z.d2.dZ.dZ(S(dHf.dw.*w.*sum(Ji.*dJ.dz2,[2 3]))) + d2z.dZ.dZ(S(dHf.dw.*w.*sum(Ji.*dJ.dz,[2 3])));
                        
    d2Sv.dZ.dZ = transpose(dw.dZ)*(Sp(sum(Ji.*dJ.dz1,[2,3]))*d2z.d1.dZ + Sp(sum(Ji.*dJ.dz2,[2,3]))*d2z.d2.dZ + Sp(sum(Ji.*dJ.dz,[2,3]))*dz.dZ)...
               + transpose(Sp(dJ.dz1(:,3,1))*dJinv.dZ(1,3))*(Sp(w)*d2z.d1.dZ) + transpose(Sp(dJ.dz2(:,3,2))*dJinv.dZ(2,3))*(Sp(w)*d2z.d2.dZ) + transpose(Sp(dJ.dz(:,3,3))*dJinv.dZ(3,3))*(Sp(w)*dz.dZ)...
                + d3z.d1.dZ.dZ(S(w.*sum(Ji.*dJ.dz1,[2,3]))) + d3z.d2.dZ.dZ(S(w.*sum(Ji.*dJ.dz2,[2,3]))) + d2z.dZ.dZ(S(w.*sum(Ji.*dJ.dz,[2,3])));
                      
    d2mT.dZ.dZ = transpose(dSv.dZ)*(d2mT.dSv.dSv)*dSv.dZ + dmT.dSv*d2Sv.dZ.dZ;
                 
    % Compute Hessians
%     d2Hf.dA.dA =  transpose(dHfd.dA)*bsxfun(@times,d2Hf.dHfd.dHfd,dHfd.dA) ... 
%                  + transpose(dHfd.dA)*bsxfun(@times,d2Hf.dHfd.dwfq,dwfq.dA) ... 
%                  + transpose(dwfq.dA)*bsxfun(@times,d2Hf.dwfq.dHfd,dHfd.dA) ... 
%                  + transpose(dwfq.dA)*bsxfun(@times,d2Hf.dwfq.dwfq,dwfq.dA) ... 
%                  + (dHf.dHfd).*d2Hfd.dA.dA + (dHf.dwfq).*d2wfq.dA.dA
    d2Hf.dA.dA =  transpose(dHfd.dA)*bsxfun(@times,d2Hf.dHfd.dHfd,dHfd.dA) ...  
                 + d2Hfd.dA.dA;
             
%     d2Hf.dA.dZ =  transpose(dHfd.dA)*bsxfun(@times,d2Hf.dHfd.dHfd,dHfd.dZ) ... 
%                  + transpose(dHfd.dA)*bsxfun(@times,d2Hf.dHfd.dwfq,dwfq.dZ) ... 
%                  + transpose(dwfq.dA)*bsxfun(@times,d2Hf.dwfq.dHfd,dHfd.dZ) ... 
%                  + transpose(dwfq.dA)*bsxfun(@times,d2Hf.dwfq.dwfq,dwfq.dZ) ... 
%                  + (dHf.dHfd).*d2Hfd.dA.dZ + (dHf.dwfq).*d2wfq.dA.dZ;
    d2Hf.dA.dZ =  transpose(dHfd.dA)*bsxfun(@times,d2Hf.dHfd.dHfd,dHfd.dZ) ... 
                 + transpose(dHfd.dA)*bsxfun(@times,d2Hf.dHfd.dw,dw.dZ) ... 
                 + d2Hfd.dA.dZ;% + (dHf.dwfq).*d2wfq.dA.dZ;<-- 0
             
%     d2Hf.dZ.dZ =  transpose(dHfd.dZ)*bsxfun(@times,d2Hf.dHfd.dHfd,dHfd.dZ) ... 
%                  + transpose(dHfd.dZ)*bsxfun(@times,d2Hf.dHfd.dwfq,dwfq.dZ) ... 
%                  + transpose(dwfq.dZ)*bsxfun(@times,d2Hf.dwfq.dHfd,dHfd.dZ) ... 
%                  + transpose(dwfq.dZ)*bsxfun(@times,d2Hf.dwfq.dwfq,dwfq.dZ) ... 
%                  + (dHf.dHfd).*d2Hfd.dZ.dZ + (dHf.dwfq).*d2wfq.dZ.dZ;
    d2Hf.dZ.dZ =  transpose(dHfd.dZ)*bsxfun(@times,d2Hf.dHfd.dHfd,dHfd.dZ) ... 
                 + transpose(dHfd.dZ)*bsxfun(@times,d2Hf.dHfd.dw,dw.dZ) ... 
                 + transpose(dw.dZ)*bsxfun(@times,d2Hf.dw.dHfd,dHfd.dZ) ... 
                 + transpose(dw.dZ)*bsxfun(@times,d2Hf.dw.dw,dw.dZ) ... 
                 + d2Hfd.dZ.dZ + d2w.dZ.dZ;

    % Store Results   
    d2r = [d2Hf.dA.dA   d2Hf.dA.dZ
           transpose(d2Hf.dA.dZ) d2Hf.dZ.dZ + sM*d2mT.dZ.dZ];
    d2r = blkdiag(d2r,diag(region.designVariables.s(:)*0),diag(region.designVariables.p*0));
    
    
end

end
