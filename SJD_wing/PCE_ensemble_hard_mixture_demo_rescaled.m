function ypreds = PCE_ensemble_hard_mixture_demo_rescaled(X, mean_Y, std_Y)

% Load the classification model
load first_hopf_aoa_1_2_classification_and_local_models net myPCE_KS

[idl,clas] = net.predict(X); % Predict labels at X

ypreds = zeros(size(X,1),1); % Pre-allocate predictions

for ii = 1:2
    [II] = find(idl == ii);
    ypreds(II,1) = uq_evalModel(myPCE_KS{ii},X(II,:))*std_Y+mean_Y;
end




