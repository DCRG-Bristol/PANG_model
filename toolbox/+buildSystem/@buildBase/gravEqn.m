function [obj] =  gravEqn(obj)

%clc;
fPath = ['+project\+grav'];

%open saved basis info..

%open saved basis info..
N_tot = obj.basis.Ntot;
varAlloc = {@(x)(-project.basis.fw(x)-project.basis.fv(x)), @(x)(project.basis.V(x)), @(x)(project.basis.W(x)), @(x)(project.basis.T(x)),...
    @(x)(project.basis.V1(x)), @(x)(project.basis.W1(x)),...
    @(x)(project.basis.W(x)), @(x)(project.basis.V(x)), @(x)(project.basis.V1(x)), @(x)(project.basis.W1(x)),...
    @(x)(project.basis.T(x)), @(x)(project.basis.T1(x))};

%open saved algebrai expressions..
exprs = open('AlgExprs\gravExpr.mat'); 
vars = exprs.varSet; disps = exprs.disps;
varDef = [disps, vars];
trms = exprs.gravTrm;  par = exprs.pars; 

%extract attachment positions from EoMs...
L = obj.geom.L;
b = obj.geom.b;
a = obj.geom.a;

%% constant vector....

for itm = 1:length(obj.grav)

    %item specific properties...
    m = obj.grav(itm).m;
    e = obj.grav(itm).e;
    elem = obj.inertia(itm).elem;

    %item object name for stiffness function...
    fcnName = obj.grav(itm).name;
    if isempty(obj.grav(itm).fctrId)
        parIdx = [];
    else
        for nameid=1:length(obj.par)
            if strcmp(obj.grav(itm).fctrId, obj.par{nameid})
                parIdx = nameid;
                break;
            end
        end
    end

    disp(['Running garvity model, item: ', fcnName, '\n'])

    for dir=1:2
        %assign null matrices...
        gravCons{dir} = zeros(N_tot,1);
        gravMat{dir}{1} = zeros(N_tot,N_tot);
        gravMat{dir}{2} = zeros(1,N_tot,N_tot,N_tot);
        gravMat{dir}{3} = zeros(N_tot,N_tot,N_tot,N_tot);

        %zero order terms.....................................................
        for numTrms=1:length(trms{dir}{1})
            %ident scaling fcators and contributing terms...
            [scalFac_sym, Var_sym] = coeffs(trms{dir}{1}(numTrms), varDef);
            scalFcn = matlabFunction(scalFac_sym, 'vars', par);
            scalF = @(x)(m(x)*scalFcn(b(x)*(a(x)-e(x))));
            for varNum=1:length(varDef)
                pow = polynomialDegree(Var_sym,varDef(varNum));
                if pow>0
                    for powCount=1:pow
                        gravCons{dir} = gravCons{dir} +...
                            integral(@(x)(scalF(x)*varAlloc{varNum}(x)),...
                            0, L, 'ArrayValued', true);
                    end
                end
            end
        end

        %zero order terms.........discrete elements
        for numTrms=1:length(trms{dir}{1})
            %ident scaling fcators and contributing terms...
            [scalFac_sym, Var_sym] = coeffs(trms{dir}{1}(numTrms), varDef);
            scalFcn = matlabFunction(scalFac_sym, 'vars', par);
            for varNum=1:length(varDef)
                pow = polynomialDegree(Var_sym,varDef(varNum));
                if pow>0
                    for powCount=1:pow
                        for elem_j = 1:length(elem)
                            xp = elem(elem_j).xp;
                            scalF = @(x)(elem(elem_j).m*scalFcn(b(x)*(a(x)-elem(elem_j).e)));
                            gravCons{dir} = gravCons{dir} + scalF(xp)*varAlloc{varNum}(xp);
                        end
                    end
                end
            end
        end
        writeFunction(gravCons{dir}, 0, [fPath,'\+odr_',num2str(0),'\'], [fcnName,'_dir_',num2str(dir)], parIdx)

        % HOT terms............................................................
        for odr=1:3

            if odr<=min([3,length(trms{dir})-1])
                for numTrms = 1:length(trms{dir}{odr+1})

                    %ident scaling fcators and contributing terms...
                    [scalFac_sym, Var_sym] = coeffs(trms{dir}{odr+1}(numTrms), varDef);
                    scalFcn = matlabFunction(scalFac_sym, 'vars', par);
                    scalF = @(x)(m(x)*scalFcn(b(x)*(a(x)-e(x))));

                    %ident the comtributing shape functionals....
                    allocPos = 1;
                    for varNum=1:length(varDef)
                        pow = polynomialDegree(Var_sym,varDef(varNum));
                        if pow>0
                            for powCount=1:pow
                                var{allocPos} = varAlloc{varNum};
                                %record the dimensionality of the assigned vector...
                                if varNum==1, dim{allocPos} = 2; else, dim{allocPos} = 1; end
                                allocPos = allocPos+1;
                            end
                        end
                    end

                    argin.odr = odr; argin.scalF = scalF;
                    spaces = zeros(1,odr+1); % indicator to record free vectors argument spaces..
                    for j=1:length(var)
                        avail = find(spaces~=1); %identify free spaces..
                        if dim{j}==1
                            argin.vctr{avail(1)} = @(x,J)(var{j}(x));
                            spaces(avail(1))=1; % lock assigned space
                        end
                        if dim{j}==2
                            argin.vctr{avail(1)} = @(x,J)(ones(N_tot,1));
                            argin.vctr{avail(2)} = @(x,J)(dimDelim(var{j}(x),J));
                            spaces(1,avail(1:2)) = 1; % lock assigned space
                        end
                    end
                    gravMat{dir}{odr} = gravMat{dir}{odr} +...
                        matCalc2(argin,N_tot,L); %call integrator...
                    clear argin; clear var; clear dim;
                end


                %%discrete elemets...
                for numTrms = 1:length(trms{dir}{odr+1})

                    %ident scaling fcators and contributing terms...
                    [scalFac_sym, Var_sym] = coeffs(trms{dir}{odr+1}(numTrms), varDef);
                    scalFcn = matlabFunction(scalFac_sym, 'vars', par);

                    %ident the comtributing shape functionals....
                    allocPos = 1;
                    for varNum=1:length(varDef)
                        pow = polynomialDegree(Var_sym,varDef(varNum));
                        if pow>0
                            for powCount=1:pow
                                var{allocPos} = varAlloc{varNum};
                                %record the dimensionality of the assigned vector...
                                if varNum==1, dim{allocPos} = 2; else, dim{allocPos} = 1; end
                                allocPos = allocPos+1;
                            end
                        end
                    end

                    for elem_j = 1:length(elem)
                        xp = elem(elem_j).xp;
                        scalF = @(x)(elem(elem_j).m*scalFcn(b(x)*(a(x)-elem(elem_j).e)));
                        argin.odr = odr; argin.scalF = scalF;
                        spaces = zeros(1,odr+1); % indicator to record free vectors argument spaces..
                        for j=1:length(var)
                            avail = find(spaces~=1); %identify free spaces..
                            if dim{j}==1
                                argin.vctr{avail(1)} = @(x,J)(var{j}(x));
                                spaces(avail(1))=1; % lock assigned space
                            end
                            if dim{j}==2
                                argin.vctr{avail(1)} = @(x,J)(ones(N_tot,1));
                                argin.vctr{avail(2)} = @(x,J)(dimDelim(var{j}(x),J));
                                spaces(1,avail(1:2)) = 1; % lock assigned space
                            end
                        end
                        gravMat{dir}{odr} = gravMat{dir}{odr} +...
                            matCalc_desc(argin,N_tot,xp); %call integrator...
                    end
                    clear argin; clear var; clear dim;
                end
            end
            writeFunction(gravMat{dir}{odr}, odr, [fPath,'\+odr_',num2str(odr),'\'], [fcnName,'_dir_',num2str(dir)], parIdx)
        end

        inpt{dir}{itm}.fcnName = [fcnName,'_dir_',num2str(dir)];

    end

end

writeOperator('vector', fPath, 'combGravMats_1', inpt{1})
writeOperator('vector', fPath, 'combGravMats_2', inpt{2})

%%
function outVec = evalMatr(cons,matr, q)
    sz = size(matr{1}); N = sz(1);
    for J=1:N
        for I=1:N
            subElem{1}(J,I) = matr{2}(1,:,J,I)*q +...
                q'*matr{3}(:,:,J,I)*q;
        end
    end
    outVec = cons + (subElem{1} + matr{1})*q;
end
%% function to output selected columns from pre-integrated matrices...
function out = dimDelim(mat, dim)
    if length(dim)==1
        out = mat(1:N_tot,dim);
    end
    if length(dim)==3 %TO DO: still bRoKeN
        out = mat(1:N_tot,dim);
    end
end

function matr = matCalc_desc(argin,N,x)
odr = argin.odr;
scalF = argin.scalF;
vctr = argin.vctr;

if odr==1
    matr = scalF(x)*(vctr{1}(x,1).*vctr{2}(x,1)');
end

if odr==2
    for I=1:N
        for J=1:N
            Ipick = [zeros(1,I-1),1,zeros(1,N-I)];
            Jpick = [zeros(1,J-1),1,zeros(1,N-J)];
            pickFcn = scalF(x)*(Jpick*vctr{1}(x,1)).*(Ipick*vctr{2}(x,1))';
            matr(1,:,J,I) = pickFcn*(vctr{3}(x,J)');
        end
    end
end

if odr==3
    for I=1:N
        for J=1:N
            Ipick = [zeros(1,I-1),1,zeros(1,N-I)];
            Jpick = [zeros(1,J-1),1,zeros(1,N-J)];
            pickFcn = scalF(x)*(Jpick*vctr{1}(x,1)).*(Ipick*vctr{2}(x,1))';
            matr(:,:,J,I) = pickFcn*(vctr{3}(x,J).*vctr{4}(x,I)');
        end
    end
end

end


end