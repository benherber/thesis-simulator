classdef test < handle
    %UNTITLED7 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Property1
    end

    methods
        function obj = test()
            %UNTITLED7 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = 0;
        end

        function res = add(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.Property1 = obj.Property1 + 1;
            res = obj.Property1;
        end
    end
end