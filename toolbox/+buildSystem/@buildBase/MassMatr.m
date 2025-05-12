function obj = MassMatr(obj)

%clc;
fPath = ['+project\+inertia'];

%open saved basis info..

N_tot = obj.basis.Ntot;
%[ut, vt, wt, tt, w1t, v1t, w1, v1, t]
varAlloc = {@(x)(-project.basis.fw(x)-project.basis.fv(x)), @(x)(project.basis.V(x)), @(x)(project.basis.W(x)),@(x)(project.basis.T(x)),...
    @(x)(project.basis.W1(x)), @(x)(project.basis.V1(x)),...
    @(x)(project.basis.W1(x)), @(x)(project.basis.V1(x)), @(x)(project.basis.T(x))};

%extract saved symbolic expressions....
exprs = open('AlgExprs\kinetExpr.mat'); trms = exprs.T_keep;
varDef = exprs.varDef; parDef = exprs.parDef; 

%extract attachment positions from EoMs...
L = obj.geom.L;
b = obj.geom.b;
a = obj.geom.a;

for itm=1:length(obj.inertia)


    %item specific properties...
    m = obj.inertia(itm).m;
    mxx = obj.inertia(itm).mxx;
    myy = obj.inertia(itm).myy;
    mzz = obj.inertia(itm).mzz;
    e = obj.inertia(itm).e;

    elem = obj.inertia(itm).elem;

    %item object name for stiffness function...
    fcnName = obj.inertia(itm).name;
    if isempty(obj.inertia(itm).fctrId)
        parIdx = [];
    else
        for nameid=1:length(obj.par)
            if strcmp(obj.inertia(itm).fctrId, obj.par{nameid})
                parIdx = nameid;
                break;
            end
        end
    end

    disp(['Running mass model, item: ', fcnName, '\n'])

    %%

    %assign null matrices...
    M{1} = zeros(N_tot,N_tot);
    M{2} = zeros(1,N_tot,N_tot,N_tot);
    M{3} = zeros(N_tot,N_tot,N_tot,N_tot);

    G{1} = zeros(N_tot,N_tot);
    G{2} = zeros(1,N_tot,N_tot,N_tot);
    G{3} = zeros(N_tot,N_tot,N_tot,N_tot);

    for odr=1:3

        %distributed bits....
        for numTrms = 1:length(trms{odr})

            %ident scaling fcators and contributing terms...
            [scalFac_sym, Var_sym] = coeffs(trms{odr}(numTrms), varDef);
            scalFcn = matlabFunction(scalFac_sym, 'vars', parDef);
            scalF = @(x)(scalFcn(m(x),mxx(x),mzz(x),myy(x),b(x),a(x),e(x)));

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
                    spaces(avail(1))=1;
                end
                if dim{j}==2
                    if spaces(1)==0
                        argin.vctr{avail(1)} = @(x,J)(ones(N_tot,1));
                        argin.vctr{avail(3)} = @(x,J)(dimDelim(var{j}(x),J));
                        spaces(1,[avail(1),avail(3)]) = 1; %lock assigmend space
                    else
                        argin.vctr{avail(1)} = @(x,J)(ones(N_tot,1));
                        argin.vctr{avail(2)} = @(x,J)(dimDelim(var{j}(x),J));
                        spaces(1,avail(1:2)) = 1; % lock assigned space
                    end
                end
            end
            Matr0 = matCalc2(argin,N_tot,L); %call integrator...
            clear argin; clear var; clear dim;

            % permuting to obtain alternative arrangemnets
            if floor(odr/2)~=odr/2
                dimOdr = [1:odr-1,odr+1,odr];
            else
                dimOdr = [1:odr,odr+2,odr+1];
            end
            M{odr} = M{odr}  + Matr0 + permute(Matr0, dimOdr);

            %Coupling matrix......
            if odr>1
                %d/dt(d[]/dqti) terms
                for vrt=odr:odr+1
                    hrz = odr-1+find([odr,odr+1]~=vrt);
                    for elm=1:odr-1
                        rm_idx = find([1:odr-1]~=elm);
                        if floor(odr/2)~=odr/2
                            dimOdr = [rm_idx,elm,vrt,hrz];
                        else
                            dimOdr = [1,rm_idx+1,elm+1,vrt+1,hrz+1];
                        end
                        G{odr} = G{odr} + permute(Matr0, dimOdr);
                    end
                end

                %-d[]/dq terms
                for vrt=1:odr-1
                    rm_idx = find([1:odr-1]~=vrt);
                    if floor(odr/2)~=odr/2
                        dimOdr = [rm_idx,odr,vrt,odr+1];
                    else
                        dimOdr = [1,rm_idx+1,odr+1,vrt+1,odr+2];
                    end
                    G{odr} = G{odr} - permute(Matr0, dimOdr);
                end
            end
        end


        %% descrete bits...
        for numTrms = 1:length(trms{odr})

            %ident scaling fcators and contributing terms...
            [scalFac_sym, Var_sym] = coeffs(trms{odr}(numTrms), varDef);
            scalFcn = matlabFunction(scalFac_sym, 'vars', parDef);

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
                scalF = @(x)(scalFcn(elem(elem_j).m, elem(elem_j).mxx,...
                    elem(elem_j).mzz, elem(elem_j).myy, b(x), a(x),...
                    elem(elem_j).e));
                argin.odr = odr; argin.scalF = scalF;
                spaces = zeros(1,odr+1); % indicator to record free vectors argument spaces..
                for j=1:length(var)
                    avail = find(spaces~=1); %identify free spaces..
                    if dim{j}==1
                        argin.vctr{avail(1)} = @(x,J)(var{j}(x));
                        spaces(avail(1))=1;
                    end
                    if dim{j}==2
                        if spaces(1)==0
                            argin.vctr{avail(1)} = @(x,J)(ones(N_tot,1));
                            argin.vctr{avail(3)} = @(x,J)(dimDelim(var{j}(x),J));
                            spaces(1,[avail(1),avail(3)]) = 1; %lock assigmend space
                        else
                            argin.vctr{avail(1)} = @(x,J)(ones(N_tot,1));
                            argin.vctr{avail(2)} = @(x,J)(dimDelim(var{j}(x),J));
                            spaces(1,avail(1:2)) = 1; % lock assigned space
                        end
                    end
                end
                if elem_j==1
                    Matr0 = matCalc_desc(argin,N_tot,xp); %call integrator...
                else
                    Matr0 = Matr0 + matCalc_desc(argin,N_tot,xp);
                end
            end
            clear argin; clear var; clear dim;

            % permuting to obtain alternative arrangemnets
            if floor(odr/2)~=odr/2
                dimOdr = [1:odr-1,odr+1,odr];
            else
                dimOdr = [1:odr,odr+2,odr+1];
            end
            M{odr} = M{odr}  + Matr0 + permute(Matr0, dimOdr);

            %Coupling matrix......
            if odr>1
                %d/dt(d[]/dqti) terms
                for vrt=odr:odr+1
                    hrz = odr-1+find([odr,odr+1]~=vrt);
                    for elm=1:odr-1
                        rm_idx = find([1:odr-1]~=elm);
                        if floor(odr/2)~=odr/2
                            dimOdr = [rm_idx,elm,vrt,hrz];
                        else
                            dimOdr = [1,rm_idx+1,elm+1,vrt+1,hrz+1];
                        end
                        G{odr} = G{odr} + permute(Matr0, dimOdr);
                    end
                end

                %-d[]/dq terms
                for vrt=1:odr-1
                    rm_idx = find([1:odr-1]~=vrt);
                    if floor(odr/2)~=odr/2
                        dimOdr = [rm_idx,odr,vrt,odr+1];
                    else
                        dimOdr = [1,rm_idx+1,odr+1,vrt+1,odr+2];
                    end
                    G{odr} = G{odr} - permute(Matr0, dimOdr);
                end
            end
        end
         disp([' .....', num2str(odr,1), ' order mass terms complete\n']);

        %%
        writeFunction(M{odr}, odr, [fPath,'\+odr_',num2str(odr),'\'], fcnName, parIdx);
        writeFunction(G{odr}, odr, [fPath,'\+odr_',num2str(odr),'\'], ['gyr_',fcnName], parIdx);
    end
    inpt{itm}.fcnName = fcnName;
    inpt_gyr{itm}.fcnName = ['gyr_',fcnName];
end

%% function to combine stiffness items...

writeOperator('matrix', fPath, 'combMassMatr', inpt)
writeOperator('matrix', fPath, 'combGyrMatr', inpt_gyr)


%% function to output selected columns from pre-integrated matrices...
function out = dimDelim(mat, dim)
    if length(dim)==1
        out = mat(1:N_tot,dim);
    end
    if length(dim)==3
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