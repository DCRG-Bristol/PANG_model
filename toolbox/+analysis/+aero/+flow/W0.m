function vel = W0(q,qt,p0,p)

persistent k1_1 k1_2 k1_3 k2_1 k2_2 k2_3 d_1 d_2 d_3

if isempty(d_1)

    k1_1 = project.aero.flow.odr_1.KW0_1(p);
    k1_2 = project.aero.flow.odr_2.KW0_1(p);
    k1_3 = project.aero.flow.odr_3.KW0_1(p);

    k2_1 = project.aero.flow.odr_1.KW0_2(p);
    k2_2 = project.aero.flow.odr_2.KW0_2(p);
    k2_3 = project.aero.flow.odr_3.KW0_2(p);

    d_1 = project.aero.flow.odr_1.DW0(p);
    d_2 = project.aero.flow.odr_2.DW0(p);
    d_3 = project.aero.flow.odr_3.DW0(p);

end

kw0_1 = k1_1;
for j=1:length(kw0_1(1,:))
    for i=1:length(kw0_1(:,1))
        kw0_1(i,j) = kw0_1(i,j) + k1_2(1,:,i,j)*q + q'*k1_3(:,:,i,j)*q;
    end
end

kw0_2 = k2_1;
for j=1:length(kw0_2(1,:))
    for i=1:length(kw0_2(:,1))
        kw0_2(i,j) = kw0_2(i,j) + k2_2(1,:,i,j)*q + q'*k2_3(:,:,i,j)*q;
    end
end

dw0 = d_1;
for j=1:length(dw0(1,:))
    for i=1:length(dw0(:,1))
        dw0(i,j) = dw0(i,j) + d_2(1,:,i,j)*q + q'*d_3(:,:,i,j)*q;
    end
end

vel = p0(1)*cos(p0(3))*(kw0_1*q) + p0(1)*sin(p0(3))*(1 + kw0_2*q) + dw0*qt; 

end