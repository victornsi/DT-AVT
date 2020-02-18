classdef Design < handle
% About: Class definition for Design
%
% Author: Victor Singh
% Lab: MIT ACDL
% Contact: victorsi@mit.edu
% /////////////////////////////////////////////////////////////////////////
    properties
        region            % Geometry Definition for Design
        inputDisturbances % Input disturbances to design
        sensorInfo        % Sensor Information
        FEMInfo           % FEM solver parameters
    end
    
    methods
        % Class constructor
        function obj = Design(region,inputDisturbances)
            % Copy Region
            r = region.copy();
            obj.region = r;
            
            % Copy Input Disturbances
            id = inputDisturbances.copy();
            obj.inputDisturbances = id;

            % Sensor Information
            obj.sensorInfo.matAllow.covariance = diag([1000 1000 2 30 2])*1e12;
            obj.sensorInfo.loads.variance = (1e-5).^2; % 10 microstrains, don't go below
            obj.sensorInfo.manufacturingParam.variance = (1/1800).^2;
            
            % FEM Information
            obj.FEMInfo.r = 1e-8;
            obj.FEMInfo.N = 30;
            
        end 
        
        % Copy class
        function design = copy(obj)
            % Create new design by copying previous
            design = Design(obj.region.copy,obj.inputDisturbances.copy);
            
            % Copy Other Parameters
            design.sensorInfo = obj.sensorInfo;
            design.FEMInfo = obj.FEMInfo;

        end
        
           
    end
    
end