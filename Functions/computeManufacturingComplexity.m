function [T,l,Mf,P]=computeManufacturingComplexity(Sa,Sp,Sv,Pinput)
% About: Manufacturing Cost Model
% Note: Time is in hours
% Convert Input Area, Perimeters, and Volume to SI units
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% /////////////////////////////////////////////////////////////////////////

% Specify Options
hasTape = 1;
depType = 'Automatic';

% Compute Total Length
stripWidth = 0.01;
Slength = Sa/stripWidth;
Splies = Sv/Sa/0.0125;

% Specify Male and Female Bends
Smale = 0;
Sfemale = 0;

% Specify Radius
Smradius = 0.01;
Sfradius = 0.01;

% Specify Basis Functions

% Clean Lay-up tool surface
ii = 1;
t(ii) = 0.0000006*Sa;
b{ii} = Sa;
p{ii} = [0.0000006];
l{ii} = 'Clean lay-up tool surface';

% Apply release agent to surface
ii = ii + 1;
t(ii) = 0.0000009*Sa;
b{ii} = Sa;
p{ii} = [0.0000009];
l{ii} = 'Apply release agent to surface';

% Position template and tape down
ii = ii + 1;
t(ii) = 0.000107*Sa.^0.77006;
b{ii} = Sa.^0.77006;
p{ii} = [0.000107];
l{ii} = 'Position template and tape down';

% Ply Deposition
ii = ii + 1;
switch depType
    case 'Manual'
        t(ii) = 0.05 + Splies*(0.001454*Slength.^0.8245);
        b{ii} = [1,Splies.*Slength.^0.8245];
        p{ii} = [0.05,0.001454];
        l{ii} = 'Ply deposition - Manual';
    case 'Hand-assist'
        t(ii) = 0.1  + Splies*(0.001585*Slength.^0.583); % Refer back to original ACCEM Cost model for power coefficient
        b{ii} = [1,Splies.*Slength.^0.583];
        p{ii} = [0.1,0.001585];
        l{ii} = 'Ply deposition - Hand-assist';
    case 'Automatic'
        t(ii) = 0.15  + Splies*(0.00063*Slength.^0.4942); % 720 IPM
        b{ii} = [1,Splies.*Slength.^0.4942];
        p{ii} = [0.15,0.00063];
        l{ii} = 'Ply deposition - Automatic';
    otherwise
        error('Check Ply Deposition Flag for Syntax')
end

% Tape Layer
% ii = ii + 1;
% if hasTape
%     t(ii) = 0.15 + Splies*(0.00063*Slength.^0.4942);
%     b{ii} = [1,Splies.*Slength.^0.4942];
%     p{ii} = [0.15,0.00063];
%     l{ii} = 'Tape layer';
% else
%     t(ii) = 0;
%     b{ii} = 1;
%     p{ii} = 0;
% end
% l{ii} = 'Tape layer';

% Transfer from Plate to Stack
ii = ii + 1;
t(ii) = Splies*(0.000145*Sa.^0.6711);
b{ii} = Splies.*Sa.^0.6711;
p{ii} = [0.000145];
l{ii} = 'Transfer from plate to stack';

% Transfer from stack to tool
ii = ii + 1;
t(ii) = 0.000145*Sa.^0.6711;
b{ii} = Sa.^0.6711;
p{ii} = [0.000145];
l{ii} = 'Transfer from stack to tool';

% Clean Curing Tool
ii = ii + 1;
t(ii) = 0.000006*Sa;
b{ii} = Sa;
p{ii} = [0.000006];
l{ii} = 'Clean curing tool';

% Apply release agent to curing tool
ii = ii + 1;
t(ii) = 0.000009*Sa;
b{ii} = Sa;
p{ii} = [0.000009];
l{ii} = 'Apply release agent to curing tool';

% Transfer lay-up to curing tool
ii = ii + 1;
t(ii) = 0.000145*Sa.^0.6711;
b{ii} = Sa.^0.6711;
p{ii} = [0.000145];
l{ii} = 'Transfer lay-up to curing tool';

% Debulking (disposable bag)
ii = ii + 1;
t(ii) = 0.02 + 0.00175*Sa.^0.6911;
b{ii} = [1,Sa.^0.6911];
p{ii} = [0.02 0.00175];
l{ii} = 'Debulking (disposable bag)';

% Sharp male bend
ii = ii + 1;
t(ii) = Splies*(0.00007*Slength)*Smale;
b{ii} = [Splies.*Slength.*Smale];
p{ii} = [0.00007];
l{ii} = 'Sharp male bend';

% Sharp female bend
ii = ii + 1;
t(ii) = Splies*(0.00016*Slength)*Sfemale;
b{ii} = [Splies.*Slength.*Sfemale];
p{ii} = [0.00016];
l{ii} = 'Sharp female bend';

% Male radial
ii = ii + 1;
if Smradius > 2
    t(ii) = 0;
    b{ii} = 1;
    p{ii} = 0;
elseif Smradius == 0
    t(ii) = 0;
    b{ii} = 1;
    p{ii} = 0;
else
    t(ii) = Splies*(0.00007*Slength)*Smale;
    b{ii} = Splies.*Slength.*Smale;
    p{ii} = 0;
end
l{ii} = 'Male radial';

% Female radial
ii = ii + 1;
if Sfradius == 0
    t(ii) = 0;
    b{ii} = 1;
    p{ii} = 0;
elseif Sfradius > 2
    t(ii) = Splies*(Slength*0.00047*Sfradius.^(-1.3585))*Sfemales;
    b{ii} = Splies.*Slength.*Sfradius.^(-1.3585).*Sfemales;
    p{ii} = 0.00047;
else
    t(ii) = Splies*(0.00016*Slength)*Smale;
    b{ii} = Splies.*Slength.*Smale;
    p{ii} = 0.00016;
end
l{ii} = 'Female radial';

% Set up
ii = ii + 1;
t(ii) = 0.07;
b{ii} = 1;
p{ii} = 0.07;
l{ii} = 'Set up';

% Gather details, prefit, disassemble, clean
ii = ii + 1;
t(ii) = 0.001326*Sa.^0.5252;
b{ii} = Sa.^0.5252;
p{ii} = 0.001326;
l{ii} = 'Details, prefit, disassemble, clean';

% Apply adhesive
ii = ii + 1;
t(ii) = 0.000055*Slength;
b{ii} = Slength;
p{ii} = 0.000055;
l{ii} = 'Apply adhesive';

% Assemble detail parts
ii = ii + 1;
t(ii) = 0.000145*Sa.^0.6711;
b{ii} = Sa.^0.6711;
p{ii} = 0.000145;
l{ii} = 'Assemble detail parts';

% Trim part
ii = ii + 1;
t(ii) = 0.00011*Sp;
b{ii} = Sp;
p{ii} = 0.00011;
l{ii} = 'Trim part';

% Apply porous separator film
ii = ii + 1;
t(ii) = 0.000009*(1.5*Sa);
b{ii} = Sa;
p{ii} = 0.000009*(1.5);
l{ii} = 'Apply porous separator film';

% Apply bleeder plies
ii = ii + 1;
t(ii) = Splies*0.00002*Sa;
b{ii} = Splies.*Sa;
p{ii} = 0.00002;
l{ii} = 'Apply bleeder plies';

% Apply non-porous separator film
ii = ii + 1;
t(ii) = 0.000009*(1.5*Sa);
b{ii} = Sa;
p{ii} = 0.000009*(1.5);
l{ii} = 'Apply non-porous separator film';

% Apply vent cloth
ii = ii + 1;
t(ii) = 0.00002*(1.5*Sa);
b{ii} = Sa;
p{ii} = 0.00002*(1.5);
l{ii} = 'Apply vent cloth';

% Install vacuum fittings
ii = ii + 1;
if Sa > 288
    t(ii) = 0.0062*2;
    b{ii} = 1;
    p{ii} = 0.0062*2;
else
    t(ii) = 0.0062;
    b{ii} = 1;
    p{ii} = 0.0062;
end
l{ii} = 'Install vacuum fittings';

% Install thermocouples
ii = ii + 1;
if Sa > 288
    t(ii) = 0.0162*2;
    b{ii} = 1;
    p{ii} = 0.0162*2;
else
    t(ii) = 0.0162;
    b{ii} = 1;
    p{ii} = 0.0162;
end
l{ii} = 'Install thermocouples';

% Apply seal strips
ii=ii + 1;
t(ii) = 0.00016*Sp;
b{ii} = Sp;
p{ii} = 0.00016;
l{ii} = 'Apply seal strips';

% Apply disposable bag
ii = ii + 1;
t(ii) = 0.000006*(1.5*Sa);
b{ii} = Sa;
p{ii} = 0.000006*(1.5);
l{ii} = 'Apply disposable bag';

% Seal edges
ii = ii + 1;
t(ii) = 0.00054*(Sp);
b{ii} = Sp;
p{ii} = 0.00054;
l{ii} = 'Seal edges';

% Connect vacuum lines, apply vacuum
ii = ii + 1;
t(ii) = 0.0061;
b{ii} = 1;
p{ii} = 0.0061;
l{ii} = 'Connect vacuum lines, apply vacuum';

% Smooth down
ii = ii + 1;
t(ii) = 0.000006*(1.5*Sa);
b{ii} = Sa;
p{ii} = 0.000006*(1.5);
l{ii} = 'Smooth down';

% Check seals
ii = ii + 1;
t(ii) = 0.000017*(Sp);
b{ii} = Sp;
p{ii} = 0.000017;
l{ii} = 'Check seals';

% Disconnect vacuum lines
ii = ii + 1;
t(ii) = 0.0031;
b{ii} = 1;
p{ii} = 0.0031;
l{ii} = 'Disconnect vacuum lines';

% Check autoclave interior
ii = ii + 1;
t(ii) = 0.03;
b{ii} = 1;
p{ii} = 0.03;
l{ii} = 'Check autoclave interior';

% Load layup-tray
ii = ii + 1;
t(ii) = 0.000145*Sa.^0.6711;
b{ii} = Sa.^0.6711;
p{ii} = 0.000145;
l{ii} = 'Load lay-up tray';

% Roll tray in
ii = ii + 1;
t(ii) = 0.025;
b{ii} = 1;
p{ii} = 0.025;
l{ii} = 'Roll tray in';

% Connect thermocouple
ii = ii + 1;
if Sa > 288
    t(ii) = 0.0092*2;
    b{ii} = 1;
    p{ii} = 0.0092*2;
else
    t(ii) = 0.0092;
    b{ii} = 1;
    p{ii} = 0.0092;
end
l{ii} = 'Connect thermocouple';

% Connect vacuum lines, apply vacuum
ii = ii + 1;
if Sa > 288
    t(ii) = 0.0061*2;
    b{ii} = 1;
    p{ii} = 0.0061*2;
else
    t(ii) = 0.0061;
    b{ii} = 1;
    p{ii} = 0.0061;
end
l{ii} = 'Connect vacuum lines, apply vacuum';

% Check bag, seal and fittings
ii = ii + 1;
t(ii) = 0.025;
b{ii} = 1;
p{ii} = 0.025;
l{ii} = 'Check bag, seal, and fittings';

% Close autoclave
ii = ii + 1;
t(ii) = 0.0192;
b{ii} = 1;
p{ii} = 0.0192;
l{ii} = 'Close autoclave';

% Set recorders
ii = ii + 1;
t(ii) = 0.056;
b{ii} = 1;
p{ii} = 0.056;
l{ii} = 'Set recorders';

% Start cure cycle and check
ii = ii + 1;
t(ii) = 0.08; % Time to cure?
b{ii} = 1;
p{ii} = 0.08;
l{ii} = 'Start cure cycle and check';

% Shutdown remove charts and open autoclave door
ii = ii + 1;
t(ii) = 0.00332 + 0.00332 + 0.0192;
b{ii} = 1;
p{ii} = 0.00332 + 0.00332 + 0.0192;
l{ii} = 'Remove charts, open autoclave door';

% Disconnect thermocouple leads
ii = ii + 1;
if Sa > 288
    t(ii) = 0.0035*2;
    b{ii} = 1;
    p{ii} = 0.0035*2;
else
    t(ii) = 0.0035;
    b{ii} = 1;
    p{ii} = 0.0035;
end
l{ii} = 'Disconnect thermocouple leads';

% Disconnect vacuum lines
ii = ii + 1;
if Sa > 288
    t(ii) = 0.0031*2;
    b{ii} = 1;
    p{ii} = 0.0031*2;
else
    t(ii) = 0.0031;
    b{ii} = 1;
    p{ii} = 0.0031;
end
l{ii} = 'Disconnect vacuum lines';

% Roll tray out of autoclave
ii = ii + 1;
t(ii) = 0.012;
b{ii} = 1;
p{ii} = 0.012;
l{ii} = 'Roll tray out of autoclave';

% Remove lay-up from tray
ii = ii + 1;
t(ii) = 0.000145*Sa.^0.6711;
b{ii} = Sa.^0.6711;
p{ii} = 0.000145;
l{ii} = 'Remove lay-up from tray';

% Remove disposable bags
ii = ii + 1;
t(ii) = 0.000008*(1.5*Sa);
b{ii} = Sa;
p{ii} = 0.000008*(1.5);
l{ii} = 'Remove disposable bags';

% Remove thermocouples
ii = ii + 1;
if Sa > 288
    t(ii) = 0.0095*2;
    b{ii} = 1;
    p{ii} = 0.0095*2;
else
    t(ii) = 0.0095;
    b{ii} = 1;
    p{ii} = 0.0095;
end
l{ii} = 'Remove thermocouples';

% Remove vacuum fittings
ii = ii + 1;
if Sa > 288
    t(ii) = 0.0029*2;
    b{ii} = 1;
    p{ii} = 0.0029*2;
else
    t(ii) = 0.0029;
    b{ii} = 1;
    p{ii} = 0.0029;
end
l{ii} = 'Remove vacuum fittings';

% Remove vent cloth
ii = ii + 1;
t(ii) = 0.000007*(1.5*Sa);
b{ii} = Sa;
p{ii} = 0.000007*(1.5);
l{ii} = 'Remove vent cloth';

% Remove non-porous separator film
ii = ii + 1;
t(ii) = 0.000007*(1.5*Sa);
b{ii} = Sa;
p{ii} = 0.000007*(1.5);
l{ii} = 'Remove non-porous separator film';

% Remove bleeder plies
ii = ii + 1;
t(ii) = 0.000007*(1.5*Sa);
b{ii} = Sa;
p{ii} = 0.000007*(1.5);
l{ii} = 'Remove bleeder plies';

% Remove porous separator film
ii = ii + 1;
t(ii) = 0.000007*(1.5*Sa);
b{ii} = Sa;
p{ii} = 0.000007*(1.5);
l{ii} = 'Remove porous separator film';

% Put Used material aside
ii = ii + 1;
t(ii) = 0.000005*(1.5*Sa);
b{ii} = Sa;
p{ii} = 0.000005*(1.5);
l{ii} = 'Put used material aside';

% Remove layup
ii = ii + 1;
t(ii) = 0.000006*(1.5*Sa);
b{ii} = Sa;
p{ii} = 0.000006*(1.5);
l{ii} = 'Remove layup';

% Clean tool
ii = ii + 1;
t(ii) = 0.000006*(1.5*Sa);
b{ii} = Sa;
p{ii} = 0.000006*(1.5);
l{ii} = 'Clean Tool';

% Form A Matrix
A = zeros(length(b),length([b{:}]));
A = blkdiag(b{:});

% Specify Time to Manufacture using P input
P = transpose([p{:}]);
T = A*Pinput;

% Store Operators
Mf.A = A;

end