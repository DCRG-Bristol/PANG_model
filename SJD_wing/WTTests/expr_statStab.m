function [true_ang, statResp, Uf, beta_yf, beta_xf] = expr_statStab(ang)

Data = open('WTTests\WT_testData\SJD_wing_rawData.mat');% experimental data must be added to this folder
Data=Data.Data;

modalData = open('WTTests\WT_testData\SJD_wing_proc_modalData.mat');
modalData = modalData.modalData;

%sensors...
sen_idx = [1,2,3,5]; %sensors to use
sen_lab = {'M_y, [Nm]', 'S_z/g, [kg]', 'M_x, [Nm]','Acc, [g]'};

%set-up a matrix to perform caliberations (same order as sen_idx)
calib_matr = eye(4);
calib_matr(4,4) = 1000/98.7; % convert acc. from V to g's

EI_ref = 2.268;
GJ_ref = 3.3843;

%% processing settings - instability detection..

instabTestSen = 3; %sensor to test for stability (amplitude threshold for LCO)
instabAmpl = 0.07; %min ampl to identify LCO


frqTestSen = 4; %sensor to estimate frq content
nondimFac = 0.5*1.225*(0.65*0.15)*0.15;

%% find closest avaiable angle..

%find closest matching data from experimental data set...
[~,p1_idx] = min(abs(Data.P1.val - ang));
true_ang = Data.P1.val(p1_idx);

%%

v_expr = Data.P2{p1_idx}.val; %experimental velo. at this angle..
dBlock = Data.P2{p1_idx}.Data;
instabTest=false;

hopfSpeed = [0,0];
hopfMy = [0,0];
hopfMx = [0,0];

[~,velOdr] = sort(v_expr);
for v_idx=1:length(v_expr)
    p2_idx=velOdr(v_idx); %just to retain generality if a different ordering is applied
    REC_data = dBlock{p2_idx}.REC{1};
    Y = (calib_matr*REC_data.Y(:,sen_idx)')';
    T = REC_data.T;

    [amplid_frq,peak_frq,Mean,poinc,prd,stab_res,UCF] = procREC(Y,T);

    for snsr=1:length(poinc)
        expr_poinc{snsr}(:,v_idx)=poinc{snsr};
    end

    expr_ampl(v_idx) = amplid_frq;
    expr_peak_frq(v_idx) = peak_frq;
    expr_equib(v_idx,:) = Mean;
    expr_LCO_prd(v_idx) = prd;
    expr_uns(v_idx,:)=UCF;

    if instabTest
        if stab_res
        else
            hopfSpeed(1,2) = v_expr(p2_idx);
            hopfMy(1,2) = expr_equib(v_idx,1)-expr_equib(1,1);
            hopfMx(1,2) = expr_equib(v_idx,3)-expr_equib(1,3);
        end
    else
        if stab_res
            hopfSpeed(1,1) = v_expr(velOdr(v_idx-1));
            hopfMy(1,1) =...
                expr_equib(velOdr(v_idx-1),1)-expr_equib(1,1);
            hopfMx(1,1) = expr_equib(v_idx-1,3)-expr_equib(1,3);
        end
    end
    instabTest=stab_res;
end

hopfSpeed(hopfSpeed==0)=NaN;
hopfMy(hopfMy==0)=NaN;
hopfMx(hopfMx==0)=NaN;

statResp.U = v_expr(velOdr);
statResp.beta_x = (expr_equib(velOdr,3)-expr_equib(1,3))/(GJ_ref);
statResp.beta_y = (expr_equib(velOdr,1)-expr_equib(1,1))/(EI_ref);

for mode=1:4
    pls = modalData.poles{mode}{p1_idx}(1:length(v_expr)); pls=pls(:)';
    statResp.frqs(mode,:) = abs(pls);
    statResp.damp(mode,:) = -real(pls)./abs(pls);
end

Uf = hopfSpeed;
beta_yf = hopfMy/(EI_ref);
beta_xf = hopfMx/(GJ_ref);



%% Identify Poincare intersections

function [amplid_frq,peak_frq,Mean,Poinc,prd,tt,ucS] = procREC(Y,T)

%amplid: max FFT amplitude
%tt: stable equiibia? (true/false)
%Mean: avegerated equilibrium
%Poinc: [Min, Max] for bifurc. info

tt = false;

for sens=1:length(Y(1,:))

    Mean(1,sens) = mean(Y(:,sens)); %mean value

    [YY,TT] = moveAVG(Y(:,sens),T);

    %estimate first derivate
    v = YY(2:end) - YY(1:end-1);

    %find 1st derivative =0
    test = v(1:end-1).*v(2:end);
    idx = find(test<0);
    Prds = TT(idx(2:2:end));
    Prds = Prds(2:end)-Prds(1:end-1);

    t_int=T(idx);

    Yp = YY(idx); %sets of stationaly points
    [idx,Poinc{sens}] = kmeans(Yp(:),2); %clustering
    Poinc{sens}=sort(Poinc{sens});

    x0=idx(1); samp=1; r1=1;
    for ii=2:length(idx)
        if idx(ii)~=x0
            tf(samp)=mean(t_int(r1:ii-1));
            x0=idx(ii);
            r1=ii;
            samp=samp+1;
        end
    end
    dt=tf(2:end)-tf(1:end-1);
    dt=dt(1:2*round(length(dt)/2)-2);
    prd=mean(dt(1:2:end)+dt(2:2:end));

    if sens==instabTestSen %LCO ?
        tt =...
            (max(Poinc{sens})-min(Poinc{sens}))>instabAmpl;
    end

    if sens==frqTestSen
        frqRang = [5, 18];
        %[P1,amp,frq] = getFFT(T,Y(:,sens));

        fs = 1./(T(2)-T(1));
        wind=round(length(T)/5);
        overlap=round(0.8*wind);

        f=linspace(1,50,300);
        [amp,frq, uc] = pwelch(Y(:,sens),wind,overlap,f,fs,'ConfidenceLevel', 0.95);

        [~,f1_idx] = min(abs(frq-frqRang(1)));
        [~,f2_idx] = min(abs(frq-frqRang(2)));


        [amplid_frq,f_idx] = max(amp(f1_idx:f2_idx));
        peak_frq = frq(f1_idx+f_idx-1);
        ucS=uc(f_idx,:)';

        amplid_frq=std(Y(:,sens));
        %prd = 1./peak_frq;
    end
end

if tt %if LCO, replace equilibrium with NaNs
    Mean = NaN*Mean;
    amplid_frq = NaN;
    % peak_frq=NaN;
else
    prd = NaN;
end

end

%% Moving avregae function...

function [Yout,Tout] = moveAVG(Y,T)

%bend and twist...
tarFrq = 10;
tarPrd = (1./tarFrq)./25;
pcs = (T(end)-T(1))/tarPrd;
gap = round(length(T)./pcs);

for sen=1:length(Y(1,:))
    for j=1:round(length(T)-gap+1)
        idx = j:j+gap-1;
        Yout(j,sen) = mean(Y(idx,sen));
        Tout(j) = mean(T(idx));
    end
end
end


%% FFT tool

function [P1,amp,frq] = getFFT(t,q)
dt = t(2)-t(1);
Fs = 1/dt;                 %Sampling frequency
L_sig = (t(end)-t(1))/dt;  % Length of signal
t = (0:L_sig-1)*dt;        % Time vector

YY = fft(q);   %Fourier transform
P2 = (YY/L_sig);
P1 = P2(1:round(L_sig/2)+1);
P1(2:end-1) = 2*P1(2:end-1);

frq = Fs*(0:(L_sig/2))/L_sig; amp = abs(P1);

end



end