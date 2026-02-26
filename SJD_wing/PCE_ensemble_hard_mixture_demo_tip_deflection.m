function ypreds = PCE_ensemble_hard_mixture_demo_tip_deflection(X)

% Load the classification model
load tip_deflection_first_hopf_aoa_1_2_classification_and_local_models net_tip_deflection myPCE_KS_tip_deflection

[idl,clas] = net_tip_deflection.predict(X); % Predict labels at X

ypreds = zeros(size(X,1),1); % Pre-allocate predictions

for ii = 1:2
    [II] = find(idl == ii);
    ypreds(II,1) = uq_evalModel(myPCE_KS_tip_deflection{ii},X(II,:));
end


