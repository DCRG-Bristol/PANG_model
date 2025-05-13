function [W1_M, W1_G] = W1_dt(q,qt,p0,p)

%such that W1_dt = W1_M*q_dt2 + W1_G

persistent d_1 d_2 d_3 g_1 g_2 g_3

if isempty(d_1)

    d_1 = project.aero.flow.odr_1.DW1(p);
    d_2 = project.aero.flow.odr_2.DW1(p);
    d_3 = project.aero.flow.odr_3.DW1(p);

    g_1 = project.aero.flow.odr_1.GW1t(p);
    g_2 = project.aero.flow.odr_2.GW1t(p);
    g_3 = project.aero.flow.odr_3.GW1t(p);

end

% W1_M = d_1;
% for j=1:length(W1_M(1,:))
%     for i=1:length(W1_M(:,1))
%         W1_M(i,j) = W1_M(i,j) + d_2(1,:,i,j)*q + q'*d_3(:,:,i,j)*q;
%     end
% end

W1_M = d_1 + squeeze(pagemtimes(d_2,q)) +...
    squeeze(pagemtimes(q',pagemtimes(d_3,q)));


%%
% W1_G = g_1;
% for j=1:length(W1_G(1,:))
%     for i=1:length(W1_G(:,1))
%         W1_G(i,j) = W1_G(i,j) + g_2(1,:,i,j)*qt + q'*g_3(:,:,i,j)*qt;
%     end
% end

W1_G = g_1 + squeeze(pagemtimes(g_2,qt)) +...
    squeeze(pagemtimes(q',pagemtimes(g_3,qt)));

W1_G = W1_G*qt;


end