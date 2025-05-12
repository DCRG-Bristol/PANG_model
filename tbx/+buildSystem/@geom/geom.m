classdef geom
    properties
        L double
        b function_handle
        a function_handle
        WTWidth = 1000000000
        sweep = 0 
        twist = @(x)(zeros(size(x)));
    end
end