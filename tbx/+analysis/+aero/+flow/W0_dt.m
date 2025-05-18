function [W0_M, W0_G] = W0_dt(q,qt,p0,p,beta)

%such that W0_dt = W0_M*q_dt2 + W0_G

persistent d_1 d_2 d_3 g_1 g_2 g_3 d1_1 d1_2 d1_3 d2_1 d2_2 d2_3 d3_1 d3_2 d3_3

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

    d3_1 = project.aero.flow.odr_1.DW0t_3(p);
    d3_2 = project.aero.flow.odr_2.DW0t_3(p);
    d3_3 = project.aero.flow.odr_3.DW0t_3(p);

end

W0_M = d_1 + squeeze(pagemtimes(d_2,q)) +...
    squeeze(pagemtimes(q',pagemtimes(d_3,q)));

W0_G = g_1 + squeeze(pagemtimes(g_2,qt)) +...
    squeeze(pagemtimes(q',pagemtimes(g_3,qt)));

W0_D1 = d1_1 + squeeze(pagemtimes(d1_2,q)) +...
    squeeze(pagemtimes(q',pagemtimes(d1_3,q)));

W0_D2 = d2_1 + squeeze(pagemtimes(d2_2,q)) +...
    squeeze(pagemtimes(q',pagemtimes(d2_3,q)));

W0_D3 = d3_1 + squeeze(pagemtimes(d3_2,q)) +...
    squeeze(pagemtimes(q',pagemtimes(d3_3,q)));


W0_G = W0_G*qt + cos(beta)*p0(1)*cos(p0(3))*(W0_D1*qt) + p0(1)*sin(p0(3))*(W0_D2*qt) +...
    sin(beta)*p0(1)*cos(p0(3))*(W0_D3*qt);




end