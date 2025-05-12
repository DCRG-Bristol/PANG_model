classdef basis%<buildSystem.buildBase

    properties
        Nw int64
        Nv int64
        Nthet int64
        xi double
    end

    properties (Dependent)
        Ngam int64
        Ntot int64
        Nstr int64
        xColloc
    end

    methods
        function xColloc = get.xColloc(obj)
            xColloc = 0.5*(obj.xi(2:end)+obj.xi(1:end-1));
        end

        function Ntot = get.Ntot(obj)
            Ntot = obj.Nw+obj.Nv+obj.Nthet;
        end
        function Ngam = get.Ngam(obj)
            Ngam = length(obj.xi)-1;
        end
    end
end