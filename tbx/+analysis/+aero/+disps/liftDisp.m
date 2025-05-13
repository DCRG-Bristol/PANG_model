function Disp = liftDisp(q,qt,p0,p,vy,vz)

%type - lift/drag/moment
%vy - strip wise velocities in O(y) direction, as a diagonal matrix
%vz - strip wise velocity in O(z) direction, as a diagonal matrix

persistent k1_1 k1_2 k1_3 k2_1 k2_2 k2_3 
if isempty(k1_1)
    k1_1 = project.aero.disps.odr_1.Lif_proj1(p);
    k1_2 = project.aero.disps.odr_2.Lif_proj1(p);
    k1_3 = project.aero.disps.odr_3.Lif_proj1(p);

    k2_1 = project.aero.disps.odr_1.Lif_proj2(p);
    k2_2 = project.aero.disps.odr_2.Lif_proj2(p);
    k2_3 = project.aero.disps.odr_3.Lif_proj2(p);

end

dx = project.aero.geom.dx_fcn(p);

% Disp1 = k1_1';
% for j=1:length(Disp1(:,1))
%     for i=1:length(Disp1(:,1))
%         Disp1(i,j) = Disp1(i,j) + k1_2(1,:,j,i)*q + q'*k1_3(:,:,j,i)*q;
%     end
% end

Disp1 = k1_1 + squeeze(pagemtimes(k1_2,q)) +...
    squeeze(pagemtimes(q',pagemtimes(k1_3,q)));

% Disp2 = k2_1';
% for j=1:length(Disp2(:,1))
%     for i=1:length(Disp2(:,1))
%         Disp2(i,j) = Disp2(i,j) + k2_2(1,:,j,i)*q + q'*k2_3(:,:,j,i)*q;
%     end
% end

Disp2 = k2_1 + squeeze(pagemtimes(k2_2,q)) +...
    squeeze(pagemtimes(q',pagemtimes(k2_3,q)));

Disp = Disp1'*dx*(-vy) + Disp2'*dx*(vz);

end