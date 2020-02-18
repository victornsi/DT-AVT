classdef Component < handle
% About: Class definition for Component
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% /////////////////////////////////////////////////////////////////////////
    properties
        design       % Inherited design information
        inputs       % Actual inputs
    end
    
    methods
        % Copy
        function component = copy(obj)
            % Make a copy
            component = Component(0,0);
            
            % Copy Design
            component.design = obj.design.copy();
            
            % Copy inputs (true)
            component.inputs = obj.inputs;
                       
        end
        
        % Class constructor
        function obj = Component(design,inputs)           
           % Initialize Design
           obj.design = design;
           
           % Initialize Inputs
           obj.inputs = inputs;
                      
        end      
        
        % Take Measurements
        function [d] =  measure(obj,tag,G)
            region = obj.design.region;
            % Get Strains
            if strcmp(tag,'operate')
                % Calculate Strains
                Nr = 1; % Number of effective repeat measurements
                
                % Complete filter matrix for region
                if size(region.field.Lf.A,1) == size(region.field.A,1)
                    s = speye(size(region.field.x,1),size(region.field.A,1));
                    region.field.Lf.A = s*(region.field.A\region.field.Lf.A);
                    region.field.Lf.b = s*(region.field.A\region.field.Lf.b);
                end
                
                % Compute Strains
                [~,~,Lf,strains] = computeStrain(region,G,[],region.sensorLocations(region.designVariables.s),1);
                d.strains = strains + 1/sqrt(Nr)*sqrt(obj.design.sensorInfo.loads.variance)*randn(size(strains));
                d.Lf = Lf;
                                
                % Calculate Sensor Selection Probabilities
                sP = region.sensorProbabilities(region.designVariables.p);
                sPi = rand(size(sP)) < sP;
                
                % Remove outputs from sensors that are removed from domain
                nonactive = region.checkInterior(region.sensorLocations(region.designVariables.s));
                nonactive = repmat(nonactive,3,1) | repmat(~sPi,3,1);
                                
                d.strains(nonactive,:) = [];               
                d.Qstrain = obj.design.sensorInfo.loads.variance*eye(size(d.strains,1));
                d.Lf.A(nonactive,:) = [];
                d.Lf.b(nonactive,:) = [];

            end
            
            % Get Failure Stresses
            if strcmp(tag,'experiment')
                d.stresses = obj.inputs.matAllow + sqrt(obj.design.sensorInfo.matAllow.covariance/4)*randn(5,1);
                d.Qstress = obj.design.sensorInfo.matAllow.covariance/4;
            end
            
            % Get Manufacturing Parameters
            if strcmp(tag,'manufacture')
                % Access Region
                time = region.complexity.mTime;
                Mf.A = region.complexity.A;
                d.times = time + sqrt(obj.design.sensorInfo.manufacturingParam.variance)*randn(size(time));
                d.Qtime = obj.design.sensorInfo.manufacturingParam.variance*eye(size(time,1));
                d.Mf = Mf;
            end
            
        end
                
    end
    
end