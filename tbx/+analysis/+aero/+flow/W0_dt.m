function [W0_M, W0_G] = W0_dt(q,qt,p0,p)

%such that W0_dt = W0_M*q_dt2 + W0_G

persistent d_1 d_2 d_3 g_1 g_2 g_3 d1_1 d1_2 d1_3 d2_1 d2_2 d2_3

if isempty(d_1)

    d_1 = project.aero.flow.odr_1.DW0(p);
    d_2 = project.aero.flow.odr_2.DW0(p);
    d_3 = project.aero.flow.odr_3.DW0(p);

    g_1 = project.aero.flow.odr_1.GW0t(p);
    g_2 = project.aero.flow.odr_2.GW0t(p);
    g_3 = project.aero.flow.odr_3.GW0t(p);

    d1_1 = project.aero.flow.odr_1.DW0t_1(p);
    d1_2 = project.aero.flow.odr_2.DW0t_1(p);
    d1_3 = project.aero.flow.odr_3.DW0t_1(p);

    d2_1 = project.aero.flow.odr_1.DW0t_2(p);
    d2_2 = project.aero.flow.odr_2.DW0t_2(p);
    d2_3 = project.aero.flow.odr_3.DW0t_2(p);

end

% W0_M = d_1;
% for j=1:length(W0_M(1,:))
%     for i=1:length(W0_M(:,1))
%         W0_M(i,j) = W0_M(i,j) + d_2(1,:,i,j)*q + q'*d_3(:,:,i,j)*q;
%     end
% end

W0_M = d_1 + squeeze(pagemtimes(d_2,q)) +...
    squeeze(pagemtimes(q',pagemtimes(d_3,q)));


%%
% W0_G = g_1;
% for j=1:length(W0_G(1,:))
%     for i=1:length(W0_G(:,1))
%         W0_G(i,j) = W0_G(i,j) + g_2(1,:,i,j)*qt + q'*g_3(:,:,i,j)*qt;
%     end
% end

W0_G = g_1 + squeeze(pagemtimes(g_2,qt)) +...
    squeeze(pagemtimes(q',pagemtimes(g_3,qt)));

% W0_D1 = d1_1;
% for j=1:length(W0_D1(1,:))
%     for i=1:length(W0_D1(:,1))
%         W0_D1(i,j) = W0_D1(i,j) + d1_2(1,:,i,j)*q + q'*d1_3(:,:,i,j)*q;
%     end
% end

W0_D1 = d1_1 + squeeze(pagemtimes(d1_2,q)) +...
    squeeze(pagemtimes(q',pagemtimes(d1_3,q)));

% W0_D2 = d2_1;
% for j=1:length(W0_D2(1,:))
%     for i=1:length(W0_D2(:,1))
%         W0_D2(i,j) = W0_D2(i,j) + d2_2(1,:,i,j)*q + q'*d2_3(:,:,i,j)*q;
%     end
% end

W0_D2 = d2_1 + squeeze(pagemtimes(d2_2,q)) +...
    squeeze(pagemtimes(q',pagemtimes(d2_3,q)));


W0_G = W0_G*qt + p0(1)*cos(p0(3))*(W0_D1*qt) + p0(1)*sin(p0(3))*(W0_D2*qt);




end