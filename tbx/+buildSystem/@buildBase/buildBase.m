classdef buildBase

    properties
        fileLocation = []
        name = 'wingModel'
        basis buildSystem.basis
        elas buildSystem.structure.elasBase
        inertia buildSystem.structure.inertiaBase
        grav buildSystem.structure.gravBase
        geom buildSystem.geom
        par0 = {'U', 'alpha0', 'alpha', 'g', 'rho'};
        par cell
        extForce %not completed yet
    end

    properties(Dependent)
        file
    end

    properties (Hidden)
        aer struct
    end

    methods
        function file = get.file(obj)
            if isempty(obj.fileLocation)
                file = [cd,'\',obj.name];
            else
                file = [obj.fileLocation,'\',obj.name];
            end   
        end

        function grav = inertia2grav(obj, inertia)
            %function to create gravity property frominertia...
            for itm=1:length(inertia)
                grav(itm) = buildSystem.structure.gravBase;
                grav(itm).m = inertia(itm).m;
                grav(itm).e = inertia(itm).e;
                grav(itm).fctrId = inertia(itm).fctrId;
                grav(itm).name = inertia(itm).name;

                for elem_j = 1:length(inertia(itm).elem)
                    elem(elem_j) = buildSystem.structure.descrWeight;
                    elem(elem_j).m = inertia(itm).elem(elem_j).m;
                    elem(elem_j).e = inertia(itm).elem(elem_j).e;
                    elem(elem_j).xp = inertia(itm).elem(elem_j).xp;
                end
                grav(itm).elem = elem;
            end
        end

        function obj = writeModel(obj)
            origCd = cd;
            cd(obj.file);
            obj = obj.genBasis;
            obj = obj.stiffMatr;
            obj = obj.MassMatr;
            obj = obj.gravEqn;
            obj = obj.AerDisps;
            obj = obj.AerMatr;
            cd(origCd);
        end

        function prepFolder(obj)
            
            mkdir([obj.file]);
            fName = [obj.file];
            
            projName = [fName, '\+project'];

            %make project folder
            mkdir(projName);

            %basis
            mkdir([projName,'\+basis']);

            %sturctural parts...
            items = {'inertia', 'grav', 'elas'};
            for i=1:length(items)
                item = [projName,'\+',items{i}];
                mkdir(item);

                for odr=0:3
                    mkdir([item, '\+odr_', num2str(odr)])
                end
            end

            %aero parts...
            items = {'disps', 'flow', 'geom'};
            for i=1:length(items)
                item = [projName,'\+aero\+',items{i}];
                mkdir(item);
                for odr=1:3
                    mkdir([item, '\+odr_', num2str(odr)])
                end
            end
        end

    end

end