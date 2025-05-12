function  obj = AerMatr(obj)

fPath = ['+project\+aero\+flow'];

%extract attachment positions from EoMs...
a = obj.geom.a; b = obj.geom.b;
WT = obj.geom.WTWidth;

%open saved basis info..
N_tot = obj.basis.Ntot;

Aer = obj.aer;
xGrid = Aer.xGrid;
NA = Aer.NA;
N_gam = length(xGrid);
xTrail = Aer.xTrail;
%X = model.Aer.X;
bx = Aer.bAer;
bx = diag(bx);
dx = Aer.dx;
dx = diag(dx);

% varAlloc = {@(x)(-basis.u{1}(x)-basis.u{2}(x)) , basis.v{1}, basis.w{1}, basis.t{1},...
%     basis.w{2}, basis.v{2},...
%     basis.w{2}, basis.v{2},  basis.t{1}};

varAlloc = {@(x)(-project.basis.fv(x)-project.basis.fw(x)) , @(x)(project.basis.V(x)), @(x)(project.basis.W(x)), @(x)(project.basis.T(x)),...
    @(x)(project.basis.W1(x)), @(x)(project.basis.V1(x)),...
    @(x)(project.basis.W1(x)), @(x)(project.basis.V1(x)), @(x)(project.basis.T(x))};


%extract saved symbolic expressions....
exprs = open('AlgExprs\AerDownWash.mat');
varDef = exprs.varDef; parDef = exprs.parDef;

%create a combined symbolic function for RHS..
Q = sym('q',[N_tot,1]); assume(Q, 'real');
Qt = sym('qt',[N_tot,1]); assume(Qt, 'real');

U = sym('U'); alp = sym('alp');
assume([U,alp], 'real');

disp(['generating aerodynamic flow terms...\n'])

%% airspeed projection.......

%create a combined symbolic function for RHS..
for dir=1:2
    trms = exprs.K_crd{dir};

    %assign null matrices...
    KK{1} = zeros(N_gam,N_tot);
    KK{2} = zeros(1,N_tot,N_gam,N_tot);
    KK{3} = zeros(N_tot,N_tot,N_gam,N_tot);

    for odr=1:3%min([3,length(trms)])
        if odr<=min([3,length(trms)])
            for numTrms = 1:length(trms{odr})

                %ident scaling fcators and contributing terms...
                [scalFac_sym, Var_sym] = coeffs(trms{odr}(numTrms), varDef);
                scalFcn = matlabFunction(scalFac_sym, 'vars', parDef);
                scalF = @(x)(scalFcn(b(x),a(x)));

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
                spaces = zeros(1,odr); % indicator to record free vectors argument spaces..
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
                Matr0 = matCalc_Aer(argin,N_tot,xGrid); %call integrator...
                clear argin; clear var; clear dim;
                KK{odr} = KK{odr}  + Matr0;

            end
        end

        writeFunction(KK{odr}, odr, [fPath,'\+odr_',num2str(odr),'\'], ['Kcrd_',num2str(dir)], [])
    end

    %     if dir==1
    %         Kcrdsym{dir} = (U*cos(alp)*(-1+evalMatr_K(KK,Q,Qt)));
    %     end
    %     if dir==2
    %         Kcrdsym{dir} = (U*sin(alp)*(evalMatr_K(KK,Q,Qt)));
    %     end

end

%Kcrd = matlabFunction(Kcrdsym{1} + Kcrdsym{2}, 'File','AeroEqns\Kcrd', 'vars', {Q,Qt,U,alp});

%% D_W0 terms...

trms = exprs.D_crd;

%assign null matrices...
DD{1} = zeros(N_gam,N_tot);
DD{2} = zeros(1,N_tot,N_gam,N_tot);
DD{3} = zeros(N_tot,N_tot,N_gam,N_tot);

for odr=1:3
    for numTrms = 1:length(trms{odr})

        %ident scaling fcators and contributing terms...
        [scalFac_sym, Var_sym] = coeffs(trms{odr}(numTrms), varDef);
        scalFcn = matlabFunction(scalFac_sym, 'vars', parDef);
        scalF = @(x)(scalFcn(b(x),a(x)));

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
        spaces = zeros(1,odr); % indicator to record free vectors argument spaces..
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
        Matr0 = matCalc_Aer(argin,N_tot,xGrid); %call integrator...
        clear argin; clear var; clear dim;
        DD{odr} = DD{odr}  + Matr0;

    end
    writeFunction(DD{odr}, odr, [fPath,'\+odr_',num2str(odr),'\'], ['Dcrd'], [])
end

% disp('  Allocating symbolics..')
%
% Dcrdsym = evalMatr_D(DD,Q,Qt);
% Dcrd = matlabFunction(Dcrdsym, 'File','AeroEqns\Dcrd','vars', {Q,Qt,U,alp});


%% D_W0 terms...

trms = exprs.D_W0_trm;

%assign null matrices...
DW0{1} = zeros(N_gam,N_tot);
DW0{2} = zeros(1,N_tot,N_gam,N_tot);
DW0{3} = zeros(N_tot,N_tot,N_gam,N_tot);

GW0t{1} = zeros(N_gam,N_tot);
GW0t{2} = zeros(1,N_tot,N_gam,N_tot);
GW0t{3} = zeros(N_tot,N_tot,N_gam,N_tot);

for odr=1:3
    for numTrms = 1:length(trms{odr})

        %ident scaling fcators and contributing terms...
        [scalFac_sym, Var_sym] = coeffs(trms{odr}(numTrms), varDef);
        scalFcn = matlabFunction(scalFac_sym, 'vars', parDef);
        scalF = @(x)(scalFcn(b(x),a(x)));

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
        spaces = zeros(1,odr); % indicator to record free vectors argument spaces..
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
        Matr0 = matCalc_Aer(argin,N_tot,xGrid); %call integrator...
        clear argin; clear var; clear dim;
        DW0{odr} = DW0{odr}  + Matr0;

        if odr>1
            %-d[]/dq terms
            for hrz=1:odr-1
                rm_idx = find([1:odr-1]~=hrz);
                if floor(odr/2)~=odr/2
                    dimOdr = [rm_idx,hrz,odr,odr+1];
                else
                    dimOdr = [1,rm_idx+1,hrz+1,odr+1,odr+2];
                end
                GW0t{odr} = GW0t{odr} + permute(Matr0, dimOdr);
            end
        end
    end

    writeFunction(DW0{odr}, odr, [fPath,'\+odr_',num2str(odr),'\'], ['DW0'], [])
    writeFunction(GW0t{odr}, odr, [fPath,'\+odr_',num2str(odr),'\'], ['GW0t'], [])
end


% DW0_sym = evalMatr_D(DW0,Q,Qt);
% GW0t_sym = evalMatr_G(GW0t,Q,Qt);
% MW0t_sym = evalMatr_M(DW0,Q,Qt);
% DW0 = matlabFunction(DW0_sym, 'File','AeroEqns\DW0','vars', {Q,Qt,U,alp});
% GW0t = matlabFunction(GW0t_sym, 'File','AeroEqns\GW0t','vars', {Q,Qt,U,alp});
% MW0t = matlabFunction(MW0t_sym, 'File','AeroEqns\MW0t','vars', {Q,Qt,U,alp});

%%
trms = exprs.D_W1_trm;

%assign null matrices...
DW1{1} = zeros(N_gam,N_tot);
DW1{2} = zeros(1,N_tot,N_gam,N_tot);
DW1{3} = zeros(N_tot,N_tot,N_gam,N_tot);

GW1t{1} = zeros(N_gam,N_tot);
GW1t{2} = zeros(1,N_tot,N_gam,N_tot);
GW1t{3} = zeros(N_tot,N_tot,N_gam,N_tot);

for odr=1:3%min([3,length(trms)])
    if odr<=min([3,length(trms)])
        for numTrms = 1:length(trms{odr})

            %ident scaling fcators and contributing terms...
            [scalFac_sym, Var_sym] = coeffs(trms{odr}(numTrms), varDef);
            scalFcn = matlabFunction(scalFac_sym, 'vars', parDef);
            scalF = @(x)(scalFcn(b(x),a(x)));

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
            spaces = zeros(1,odr); % indicator to record free vectors argument spaces..
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
            Matr0 = matCalc_Aer(argin,N_tot,xGrid); %call integrator...
            clear argin; clear var; clear dim;
            DW1{odr} = DW1{odr}  + Matr0;

            if odr>1
                %-d[]/dq terms
                for hrz=1:odr-1
                    rm_idx = find([1:odr-1]~=hrz);
                    if floor(odr/2)~=odr/2
                        dimOdr = [rm_idx,hrz,odr,odr+1];
                    else
                        dimOdr = [1,rm_idx+1,hrz+1,odr+1,odr+2];
                    end
                    GW1t{odr} = GW1t{odr} + permute(Matr0, dimOdr);
                end
            end
        end
    end

    writeFunction(DW1{odr}, odr, [fPath,'\+odr_',num2str(odr),'\'], ['DW1'], [])
    writeFunction(GW1t{odr}, odr, [fPath,'\+odr_',num2str(odr),'\'], ['GW1t'], [])
end


%create a combined symbolic function for RHS..
% DW1_sym = evalMatr_D(DW1,Q,Qt);
% GW1t_sym = evalMatr_G(GW1t,Q,Qt);
% MW1t_sym = evalMatr_M(DW1,Q,Qt);
% DW1 = matlabFunction(DW1_sym, 'File','AeroEqns\DW1','vars', {Q,Qt,U,alp});
% GW1t = matlabFunction(GW1t_sym, 'File','AeroEqns\GW1t','vars', {Q,Qt,U,alp});
% MW1t = matlabFunction(MW1t_sym, 'File','AeroEqns\MW1t','vars', {Q,Qt,U,alp});

%%

%create a combined symbolic function for RHS..
for dir=1:2
    trms = exprs.K_W0_trm{dir};

    %assign null matrices...
    KW0{1} = zeros(N_gam,N_tot);
    KW0{2} = zeros(1,N_tot,N_gam,N_tot);
    KW0{3} = zeros(N_tot,N_tot,N_gam,N_tot);

    %assign null matrices...
    DW0t{1} = zeros(N_gam,N_tot);
    DW0t{2} = zeros(1,N_tot,N_gam,N_tot);
    DW0t{3} = zeros(N_tot,N_tot,N_gam,N_tot);

    for odr=1:3%min([3,length(trms)])
        if odr<=min([3,length(trms)])
            for numTrms = 1:length(trms{odr})

                %ident scaling fcators and contributing terms...
                [scalFac_sym, Var_sym] = coeffs(trms{odr}(numTrms), varDef);
                scalFcn = matlabFunction(scalFac_sym, 'vars', parDef);
                scalF = @(x)(scalFcn(b(x),a(x)));

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
                spaces = zeros(1,odr); % indicator to record free vectors argument spaces..
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
                Matr0 = matCalc_Aer(argin,N_tot,xGrid); %call integrator...
                clear argin; clear var; clear dim;
                KW0{odr} = KW0{odr}  + Matr0;

                %-d[]/dq terms
                selSet = [1:odr-1,odr+1];
                for hrz=selSet
                    rm_idx = find(selSet~=hrz);
                    rm_idx = selSet(rm_idx);
                    if floor(odr/2)~=odr/2
                        dimOdr = [rm_idx,odr,hrz];
                    else
                        dimOdr = [1,rm_idx+1,odr+1,hrz+1];
                    end
                    DW0t{odr} = DW0t{odr} + permute(Matr0, dimOdr);
                end
            end
        end
        writeFunction(KW0{odr}, odr, [fPath,'\+odr_',num2str(odr),'\'], ['KW0_',num2str(dir)], []);
        writeFunction(DW0t{odr}, odr, [fPath,'\+odr_',num2str(odr),'\'], ['DW0t_',num2str(dir)], [])
    end

    %     if dir==1
    %         KW0_sym{dir} = (U*cos(alp)*evalMatr_K(KW0,Q,Qt));
    %         DW0t_sym{dir} = (U*cos(alp)*evalMatr_D(DW0t,Q,Qt));
    %     end
    %     if dir==2
    %         DW0t_sym{dir} = (U*sin(alp)*evalMatr_D(DW0t,Q,Qt));
    %         KW0_sym{dir} = (U*sin(alp)*(1+evalMatr_K(KW0,Q,Qt)));
    %     end

end

%
% KW0 = matlabFunction(KW0_sym{1} + KW0_sym{2}, 'File','AeroEqns\KW0', 'vars', {Q,Qt,U,alp});
% DW0t = matlabFunction(DW0t_sym{1} + DW0t_sym{2}, 'File','AeroEqns\DW0t', 'vars', {Q,Qt,U,alp});

%% function assign coordinates to function...

    function outMat = evalMatr_M(matr, q, qt)
        sz = size(matr{1}); N = sz(2);
        for J=1:length(xGrid)
            for I=1:N
                subElem{1}(J,I) = matr{2}(1,:,J,I)*q +...
                    q'*matr{3}(:,:,J,I)*q;
            end
        end
        outMat = (subElem{1} + matr{1});
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function outVec = evalMatr_D(matr, q, qt)
        sz = size(matr{1}); N = sz(2);
        for J=1:length(xGrid)
            for I=1:N
                subElem{1}(J,I) = matr{2}(1,:,J,I)*q +...
                    q'*matr{3}(:,:,J,I)*q;
            end
        end
        outVec = (subElem{1} + matr{1})*qt;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function outVec = evalMatr_K(matr, q, qt)
        sz = size(matr{1}); N = sz(2);
        for J=1:length(xGrid)
            for I=1:N
                subElem{1}(J,I) = matr{2}(1,:,J,I)*q +...
                    q'*matr{3}(:,:,J,I)*q;
            end
        end
        outVec = (subElem{1} + matr{1})*q;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function outVec = evalMatr_G(matr, q, qt)
        sz = size(matr{1}); N = sz(2);
        for J=1:length(xGrid)
            for I=1:N
                subElem{1}(J,I) = matr{2}(1,:,J,I)*qt +...
                    q'*matr{3}(:,:,J,I)*qt;
            end
        end
        outVec = (subElem{1} + matr{1})*qt;
    end

%% function to output selected columns from pre-integrated matrices...
    function out = dimDelim(mat, dim)
        if length(dim)==1
            out = mat(1:N_tot,dim);
        end
        if length(dim)==3
            out = mat(1:N_tot,dim);
        end
    end
end