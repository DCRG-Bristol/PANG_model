classdef oneraBase<analysis.analysisBase
    properties (Dependent)
        unsCircIdx = []
        unsAeroIdx
        q0_aeroStruct
    end

    properties
        Cl_grad = @(alp, U)(2*pi*ones(size(alp)));
        Cl = @(alp, U)(2*pi*alp);
        Cm = @(alp, U)(zeros(size(alp)));
        Cm_grad = @(alp, U)(zeros(size(alp)))
        Cd = @(alp, U)(zeros(size(alp)));

        %%onera stuff....
        lam = 0.17;
        ML = 0.53;

        %pressure loss/scaling for aero laods (e.g. fairing gaps)
        lossFactor = 0;
    end

    %this is temporarily hidden.. until the option for propegation eqn
    %without uLLM is included
    properties (Hidden)
        isUnstCirc = true;
    end

    methods
        function unsAeroIdx = get.unsAeroIdx(obj)
            Nstr = length(obj.transF(1,:));
            Ngam = obj.basis.Ngam;
            unsAeroIdx = [2*obj.Nstr+1:2*Nstr+Ngam];
        end

        function unsCircIdx = get.unsCircIdx(obj)
            Nstr = length(obj.transF(1,:));
            Ngam = obj.basis.Ngam;
            if obj.isUnstCirc
                unsCircIdx = [2*Nstr+Ngam+1:2*Nstr+2*Ngam];
            else
                unsCircIdx = [];
            end
        end
        function q0_aeroStruct = get.q0_aeroStruct(obj)
            if obj.isUnstCirc
                q0_aeroStruct = zeros(2*obj.Nstr+2*obj.basis.Ngam,1);
            else
                q0_aeroStruct = zeros(2*obj.Nstr+obj.basis.Ngam,1);
            end
        end

        %% function for derivaties...

        function qdt = aero_structDeriv(obj,q_all,varargin)
            [M_glob,F_glob] = obj.aero_structModel(q_all,varargin{:});
            qdt = M_glob\F_glob;
        end


        %%
        function [M_glob,F_glob] = aero_structModel(obj,q_all,varargin)

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

            if isempty(obj.unsCircIdx)
                M_glob = zeros(2*obj.Nstr+obj.basis.Ngam, 2*obj.Nstr+obj.basis.Ngam);
                F_glob = zeros(2*obj.Nstr+obj.basis.Ngam, 1);
            else
                M_glob = zeros(2*obj.Nstr+2*obj.basis.Ngam, 2*obj.Nstr+2*obj.basis.Ngam);
                F_glob = zeros(2*obj.Nstr+2*obj.basis.Ngam, 1);
            end
            [M_glob,F_glob] = obj.buildStruct(M_glob,F_glob,q_all,p0,p);
            [M_glob,F_glob] = obj.buildAero(M_glob,F_glob,q_all,p0,p);

        end

        %%

        function [eval, evec, isFlut] = AE_modes(obj, q0, varargin)

            epsilon = 1e-6;

            %initial call...
            [M0,F0] = obj.aero_structModel(q0,varargin{:});
            statSize = size(q0, 1);
            jac = zeros(statSize, statSize);

            for j=1:statSize
                q=q0; q(j) = q(j)+epsilon;
                [~,F] = obj.aero_structModel(q,varargin{:});
                jac(:,j) = (F - F0)/epsilon;
            end

            %run eigenvalue check...
            [V,D] = eig(jac, M0); D = diag(D);
            %find oscillatory modes..
            osc_idx = find(abs(imag(D))>0.01*abs(real(D)));
            D = D(osc_idx);
            V = V(:,osc_idx);

            %sort by frequency..
            [~,odr] = sort(abs(D));
            eval=D(odr(2:2:end));
            evec=V(:,odr(2:2:end));

            %test for flutter...
            isFlut = ~isempty(find(real(eval)>0, 1));

        end

        %% aero coeffs...

        function [Cf, Cf_alp] = Cl_fcn(obj, vz, U_tot)
            Cf = obj.Cl(vz./U_tot + obj.geom.twist(obj.basis.xColloc(:)),...
                U_tot);
            Cf_alp = obj.Cl_grad(vz./U_tot + obj.geom.twist(obj.basis.xColloc(:)),...
                U_tot);
        end

        function [Cf, Cf_alp] = Cd_fcn(obj, vz, U_tot)
            Cf = obj.Cd(vz./U_tot + obj.geom.twist(obj.basis.xColloc(:)),...
                U_tot);
            Cf_alp = 0*ones(size(vz));
        end

        function [Cf, Cf_alp] = Cm_fcn(obj, vz, U_tot)
            Cf = obj.Cm(vz./U_tot + obj.geom.twist(obj.basis.xColloc(:)),...
                U_tot);
            Cf_alp = obj.Cm_grad(vz./U_tot + obj.geom.twist(obj.basis.xColloc(:)),...
                U_tot);
        end

    end

end