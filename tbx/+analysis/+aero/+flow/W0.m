function vel = W0(q,qt,p0,p,beta)

persistent k1_1 k1_2 k1_3 k2_1 k2_2 k2_3 k3_1 k3_2 k3_3 d_1 d_2 d_3

if isempty(d_1)

    k1_1 = project.aero.flow.odr_1.KW0_1(p);
    k1_2 = project.aero.flow.odr_2.KW0_1(p);
    k1_3 = project.aero.flow.odr_3.KW0_1(p);

    k2_1 = project.aero.flow.odr_1.KW0_2(p);
    k2_2 = project.aero.flow.odr_2.KW0_2(p);
    k2_3 = project.aero.flow.odr_3.KW0_2(p);

    k3_1 = project.aero.flow.odr_1.KW0_3(p);
    k3_2 = project.aero.flow.odr_2.KW0_3(p);
    k3_3 = project.aero.flow.odr_3.KW0_3(p);

    d_1 = project.aero.flow.odr_1.DW0(p);
    d_2 = project.aero.flow.odr_2.DW0(p);
    d_3 = project.aero.flow.odr_3.DW0(p);

end

kw0_1 = k1_1 + squeeze(pagemtimes(k1_2,q)) +...
    squeeze(pagemtimes(q',pagemtimes(k1_3,q)));

kw0_2 = k2_1 + squeeze(pagemtimes(k2_2,q)) +...
    squeeze(pagemtimes(q',pagemtimes(k2_3,q)));

kw0_3 = k3_1 + squeeze(pagemtimes(k3_2,q)) +...
    squeeze(pagemtimes(q',pagemtimes(k3_3,q)));

dw0 = d_1 + squeeze(pagemtimes(d_2,q)) +...
    squeeze(pagemtimes(q',pagemtimes(d_3,q)));

vel = cos(beta)*p0(1)*cos(p0(3))*(kw0_1*q) + p0(1)*sin(p0(3))*(1 + kw0_2*q) +...
    sin(beta)*p0(1)*cos(p0(3))*(kw0_3*q) + dw0*qt; 

end