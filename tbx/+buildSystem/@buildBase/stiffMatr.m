function obj = stiffMatr(obj)

%clc;
fPath = ['+project\+elas'];

%open saved basis info..
N_tot = obj.basis.Ntot;
varAlloc = {@(x)(project.basis.W1(x)), @(x)(project.basis.W2(x)), @(x)(project.basis.V1(x)),@(x)(project.basis.V2(x)),...
    @(x)(project.basis.T(x)), @(x)(project.basis.T1(x))};

%extract saved symbolic expressions....

exprs = load('stiffExpr.mat');
trms = exprs.U_keep;
varDef = exprs.varDef; parDef = exprs.parDef; 

%extract attachment positions from EoMs...
L = obj.geom.L;

for itm=1:length(obj.elas)


    %item specific properties...
    EI1 = obj.elas(itm).EI1;
    EI2 = obj.elas(itm).EI2;
    EI12 = obj.elas(itm).EI12;
    GJ = obj.elas(itm).GJ;

    %item object name for stiffness function...
    fcnName = obj.elas(itm).name;

    if isempty(obj.elas(itm).fctrId)
        parIdx = [];
    else
        for nameid=1:length(obj.par)
            if strcmp(obj.elas(itm).fctrId, obj.par{nameid})
                parIdx = nameid;
                break;
            end
        end
    end

    disp(['Running stiffness model, item: ', fcnName, '\n'])

    %%
    K{1} = zeros(N_tot,N_tot);
    K{2} = zeros(1,N_tot,N_tot,N_tot);
    K{3} = zeros(N_tot,N_tot,N_tot,N_tot);
    K{3} = zeros(N_tot,N_tot,N_tot,N_tot);
    K{4} = zeros(1,N_tot,N_tot,N_tot,N_tot,N_tot);
    K{5} = zeros(N_tot,N_tot,N_tot,N_tot,N_tot,N_tot);

    for odr=1:3
        for numTrms = 1:length(trms{odr})
            [scalFac_sym, Var_sym] = coeffs(trms{odr}(numTrms), varDef);
            scalFcn = matlabFunction(scalFac_sym, 'vars', parDef);
            scalF = @(x)(scalFcn(EI1(x),EI2(x),EI12(x),GJ(x)));

            %ident the comtributing shape functionals....
            allocPos = 1;
            for varNum=1:length(varDef)
                pow = polynomialDegree(Var_sym,varDef(varNum));
                if pow>0
                    for powCount=1:pow
                        var{allocPos} = varAlloc{varNum};
                        allocPos = allocPos+1;
                    end
                end
            end

            argin.odr = odr; argin.scalF = scalF;
            for j=1:length(var)
                argin.vctr{j} = @(x,J)(var{j}(x));
            end
            Matr0 = matCalc2(argin,N_tot,L); %call integrator...
            clear argin; clear var;

            idx = 1:1:odr+1;
            for vrt = 1:odr+1
                exclIdx = find(idx~=vrt);
                if floor(odr/2)~=odr/2
                    dimOdr = [idx(exclIdx(1:end-1)), vrt, idx(exclIdx(end))];
                else
                    dimOdr = [1,1+idx(exclIdx(1:end-1)), 1+vrt, 1+idx(exclIdx(end))];
                end
                K{odr} = K{odr}  + permute(Matr0, dimOdr);
            end
        end
        disp([' .....', num2str(odr,1), ' order stiffness terms complete\n']);
        writeFunction(K{odr}, odr, [fPath,'\+odr_',num2str(odr),'\'], fcnName, parIdx)
    end

    inpt{itm}.fcnName = fcnName;

end

%% function to combine stiffness items...

writeOperator('matrix', fPath, 'combElasMatr', inpt)

end