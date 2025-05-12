function vel = Vy(q,qt,p0,p)

persistent k1_1 k1_2 k1_3 k2_1 k2_2 k2_3 d_1 d_2 d_3

if isempty(d_1)

    k1_1 = project.aero.flow.odr_1.Kcrd_1(p);
    k1_2 = project.aero.flow.odr_2.Kcrd_1(p);
    k1_3 = project.aero.flow.odr_3.Kcrd_1(p);

    k2_1 = project.aero.flow.odr_1.Kcrd_2(p);
    k2_2 = project.aero.flow.odr_2.Kcrd_2(p);
    k2_3 = project.aero.flow.odr_3.Kcrd_2(p);

    d_1 = project.aero.flow.odr_1.Dcrd(p);
    d_2 = project.aero.flow.odr_2.Dcrd(p);
    d_3 = project.aero.flow.odr_3.Dcrd(p);
end

kcrd1 = k1_1;
for j=1:length(kcrd1(1,:))
    for i=1:length(kcrd1(:,1))
        kcrd1(i,j) = kcrd1(i,j) + k1_2(1,:,i,j)*q + q'*k1_3(:,:,i,j)*q;
    end
end

kcrd2 = k2_1;
for j=1:length(kcrd2(1,:))
    for i=1:length(kcrd2(:,1))
        kcrd2(i,j) = kcrd2(i,j) + k2_2(1,:,i,j)*q + q'*k2_3(:,:,i,j)*q;
    end
end

dcrd = d_1;
for j=1:length(dcrd(1,:))
    for i=1:length(dcrd(:,1))
        dcrd(i,j) = dcrd(i,j) + d_2(1,:,i,j)*q + q'*d_3(:,:,i,j)*q;
    end
end

vel = p0(1)*cos(p0(3))*(-1 + kcrd1*q) + p0(1)*sin(p0(3))*(kcrd2*q) + dcrd*qt; 

end