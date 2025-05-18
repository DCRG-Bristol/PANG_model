function vel = Vy(q,qt,p0,p,beta)

persistent k1_1 k1_2 k1_3 k2_1 k2_2 k2_3 k3_1 k3_2 k3_3 d_1 d_2 d_3

if isempty(d_1)

    k1_1 = project.aero.flow.odr_1.Kcrd_1(p);
    k1_2 = project.aero.flow.odr_2.Kcrd_1(p);
    k1_3 = project.aero.flow.odr_3.Kcrd_1(p);

    k2_1 = project.aero.flow.odr_1.Kcrd_2(p);
    k2_2 = project.aero.flow.odr_2.Kcrd_2(p);
    k2_3 = project.aero.flow.odr_3.Kcrd_2(p);

    k3_1 = project.aero.flow.odr_1.Kcrd_3(p);
    k3_2 = project.aero.flow.odr_2.Kcrd_3(p);
    k3_3 = project.aero.flow.odr_3.Kcrd_3(p);

    d_1 = project.aero.flow.odr_1.Dcrd(p);
    d_2 = project.aero.flow.odr_2.Dcrd(p);
    d_3 = project.aero.flow.odr_3.Dcrd(p);
end

kcrd1 = k1_1 + squeeze(pagemtimes(k1_2,q)) +...
    squeeze(pagemtimes(q',pagemtimes(k1_3,q)));

kcrd2 = k2_1 + squeeze(pagemtimes(k2_2,q)) +...
    squeeze(pagemtimes(q',pagemtimes(k2_3,q)));

kcrd3 = k3_1 + squeeze(pagemtimes(k3_2,q)) +...
    squeeze(pagemtimes(q',pagemtimes(k3_3,q)));

dcrd = d_1 + squeeze(pagemtimes(d_2,q)) +...
    squeeze(pagemtimes(q',pagemtimes(d_3,q)));


vel = cos(beta)*p0(1)*cos(p0(3))*(-1 + kcrd1*q) + p0(1)*sin(p0(3))*(kcrd2*q) +...
    sin(beta)*p0(1)*cos(p0(3))*(kcrd3*q) + dcrd*qt; 

end