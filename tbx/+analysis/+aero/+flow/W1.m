function vel = W1(q,qt,p0,p)

persistent d_1 d_2 d_3

if isempty(d_1)
    d_1 = project.aero.flow.odr_1.DW1(p);
    d_2 = project.aero.flow.odr_2.DW1(p);
    d_3 = project.aero.flow.odr_3.DW1(p);
end

% dw1 = d_1;
% for j=1:length(dw1(1,:))
%     for i=1:length(dw1(:,1))
%         dw1(i,j) = dw1(i,j) + d_2(1,:,i,j)*q + q'*d_3(:,:,i,j)*q;
%     end
% end

dw1 = d_1 + squeeze(pagemtimes(d_2,q)) +...
    squeeze(pagemtimes(q',pagemtimes(d_3,q)));

vel = dw1*qt;

end