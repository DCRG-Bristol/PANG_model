classdef extAeroBase<analysis.analysisBase

    properties (Dependent)
        q0_aeroStruct
    end

    methods

        function q0_aeroStruct = get.q0_aeroStruct(obj)
            q0_aeroStruct = zeros(2*obj.Nstr,1);
        end

        function [L0, D0, M0] = initLoads(obj)
            L0 = ones(obj.basis.Ngam,1);
            D0 = ones(obj.basis.Ngam,1);
            M0 = ones(obj.basis.Ngam,1);
        end

        function qdt = aero_structDeriv(obj, q_all, L, D, M, ang, varargin)

            %reference parametert values
            p = obj.p;
            p0 = obj.p0;

            %update with user provided values..
            if isempty(varargin)
            else
                for p_inpt=1:length(varargin)/2
                    idx0 = find(strcmp(varargin{2*p_inpt-1}, obj.par0));
                    if isempty(idx0)
                        idx = find(strcmp(varargin{2*p_inpt-1}, obj.par));
                        if isempty(idx)
                            error(['Undefined variable ',varargin{2*p_inpt-1}]);
                        else
                            p(idx) = varargin{2*p_inpt};
                        end
                    else
                        p0(idx0) = varargin{2*p_inpt};
                    end
                end
            end

            %assign memory
            M_glob = zeros(2*obj.Nstr, 2*obj.Nstr);
            F_glob = zeros(2*obj.Nstr, 1);
            [M_glob,F_glob] = obj.buildStruct(M_glob,F_glob, q_all,p0,p);
            [M_glob,F_glob] = obj.buildAero(M_glob,F_glob, L, D, M, ang, q_all,p0,p);

            qdt = M_glob\F_glob;
        end

        %% function to add mesh...
        function [x, y, z] = getMesh(obj,q_all,frame)
            [x, y, z] = analysis.plotFcns.wingDefl(obj, q_all, obj.basis.xi);
            
            availTypes = {'aircraft', 'beamModel'};
            if isempty(find(strcmp(frame, availTypes)))
                error(['unexpected frame type ',frame,', Expected: aircraft or beamModel']);
            else
                switch frame
                    case 'beamModel'
                    case 'aircraft'
                        R(1,:) = [y(1,:), y(2,:)];
                        R(2,:) = [x(1,:), x(2,:)];
                        R(3,:) = [-z(1,:), -z(2,:)];

                        beta = obj.geom.sweep;
                        R = [cos(beta), -sin(beta), 0;...
                            sin(beta), cos(beta), 0;...
                            0, 0, 1]*R;

                        x = [R(1,1:end/2); R(1,end/2+1:end)];
                        y = [R(2,1:end/2); R(2,end/2+1:end)];
                        z = [R(3,1:end/2); R(3,end/2+1:end)];
                end
            end
        end


    end

end