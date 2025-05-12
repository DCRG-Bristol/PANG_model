classdef analysisBase

    properties
        transF
        par
        par0
        p = [];
        p0 = [20, 0, 0, 9.81, 1.225];
        dampMatr
    end

    properties (SetAccess = {?buildSystem.buildBase})
        basis
        geom
        file
    end

    properties (Dependent)
        Nstr
        structDisp
        structVel
        q0_struct
    end

    methods

        function obj = analysisBase(obj_in)
            arguments
                obj_in buildSystem.buildBase
            end
            obj.basis = obj_in.basis;
            obj.geom = obj_in.geom;
            obj.transF = eye(obj_in.basis.Ntot,obj_in.basis.Ntot);
            obj.par = obj_in.par;
            obj.par0 = obj_in.par0;
            obj.dampMatr = zeros(obj_in.basis.Ntot);
            obj.file = obj_in.file;
        end

        %%
%                 function obj = set.dampMatr(obj,val)
%                     sz = size(val);
%                     if sz(1)~=obj.basis.Ntot
%                         error(['Incorrect size for damping matrix, must be: ', num2str(obj.basis.Ntot), ' by ', num2str(obj.basis.Ntot)]);
%                     else
%                         if sz(2)~=obj.basis.Ntot
%                             error(['Incorrect size for damping matrix, must be: ', num2str(obj.basis.Ntot), ' by ', num2str(obj.basis.Ntot)]);
%                         end
%                     end
%                     obj.dampMatr = val;
%                 end

        function q0_struct = get.q0_struct(obj)
            q0_struct = zeros(2*obj.Nstr,1);
        end

        function Nstr = get.Nstr(obj)
            Nstr = length(obj.transF(1,:));
        end

        function structDisp = get.structDisp(obj)
            Nstr = length(obj.transF(1,:));
            structDisp = [1:Nstr];
        end
        function structVel = get.structVel(obj)
            Nstr = length(obj.transF(1,:));
            structVel = [Nstr+1:2*Nstr];
        end

        %%
        function obj = setTransform(obj, type, DoF)

            if nargin<3 %if not specified, use the complete set of DoFs
                DoF = obj.basis.Ntot;
            end

            if DoF>obj.basis.Ntot
                error('Requested degrees of freedom exceeds model size');
            end

            if strcmp(type, 'default')
                if DoF>obj.basis.Ntot
                    DoF = obj.basis.Ntot;
                    warning(['default transformation requires the entire model size...re setting to ',num2str(obj.basis.Ntot)])
                end
            end

            %first, re-set to default..
            obj.transF = eye(obj.basis.Ntot);

            %now change..
            switch type
                case 'default'
                    obj.transF = [eye(DoF); zeros(obj.basis.Ntot - DoF, DoF)];
                case 'modal'
                    obj_temp = obj;
                    obj_temp.dampMatr = zeros(obj.basis.Ntot);
                    [shp, evals] = getStructModes(obj_temp,zeros(2*obj.basis.Ntot,1));
                    for j=1:DoF
                        mc_r=real(shp(:,j));
                        mc_i=imag(shp(:,j));
                        [uu,ss,vv]=svd([mc_r,mc_i]');
                        dq=uu(:,1).'*[mc_r.';mc_i.'];   % real mode
                        shpMat(:,j) = dq';
                    end
                    obj.transF = shpMat;
            end
        end

        %% function for eigenvalue analysis - structural..
        function [shp, evals, K, C, M] = getStructModes(obj,q0,varargin)

            [Mglob_0,Fglob_0] = obj.structModel(q0, varargin{:});
            epsilon = 1e-6;

            K = zeros(obj.Nstr, obj.Nstr);
            M = Mglob_0(obj.Nstr+1:end, obj.Nstr+1:end);
            C = obj.transF'*obj.dampMatr*obj.transF;

            for j=1:obj.Nstr
                q=q0; q(j) = q(j)+epsilon;
                [~,Fglob] = obj.structModel(q,varargin{:});
                K(:,j) = -(Fglob(obj.Nstr+1:end) - Fglob_0(obj.Nstr+1:end))/epsilon;
            end

            [V,D] = polyeig(K, C, M);

            [~,idx] = sort(abs(D));
            for vecNum = 1:length(V(1,:))
                %normalise
                normFac = sqrt(V(:,vecNum)'*M*V(:,vecNum));
                V(:,vecNum) = V(:,vecNum)./normFac; %mass normalised vectors
            end
            evals = D(idx); shp = V(:,idx);

            evals = evals(2:2:end); shp = shp(:,2:2:end);
        end


        %%
        function obj = setPars(obj,varargin)

            %reference parametert values
            p = zeros(1,length(obj.par));
            p0 =  obj.p0;

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

            obj.p = p;
            obj.p0 = p0;

        end

        %%

        function qdt = structDeriv(obj,q_all,varargin)
            [M_glob,F_glob] = obj.structModel(q_all,varargin{:});
            qdt = M_glob\F_glob;
        end

        function [M_glob,F_glob] = structModel(obj,q_all,varargin)

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
            [M_glob,F_glob] = obj.buildStruct(M_glob,F_glob,q_all,p0,p);

        end

        %% function to get displacements..
        function [x, y, z] = getDisplField(obj,q_all,xStat,frame)
            
            [x, y, z] = analysis.plotFcns.wingDefl(obj, q_all, xStat);

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
                        R = [cos(beta), sin(beta), 0;...
                            -sin(beta), cos(beta), 0;...
                            0, 0, 1]*R;

                        x = [R(1,1:end/2); R(1,end/2+1:end)];
                        y = [R(2,1:end/2); R(2,end/2+1:end)];
                        z = [R(3,1:end/2); R(3,end/2+1:end)];
                end
            end
        end
    end

    methods (Static)
        function obj = loadobj(obj)
            addpath([obj.file], '-begin');
            clc;
            fprintf(['...Accessing matrices found in %s'], obj.file);
        end
    end
end

