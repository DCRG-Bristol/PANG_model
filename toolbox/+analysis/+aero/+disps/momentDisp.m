function Disp = momentDisp(q,qt,p0,p,vy,vz)

%type - lift/drag/moment
%vy - strip wise velocities in O(y) direction, as a diagonal matrix
%vz - strip wise velocity in O(z) direction, as a diagonal matrix

persistent k_1 k_2 k_3

if isempty(k_1)
    k_1 = project.aero.disps.odr_1.Mom_Proj(p);
    k_2 = project.aero.disps.odr_2.Mom_Proj(p);
    k_3 = project.aero.disps.odr_3.Mom_Proj(p);
end

dx = project.aero.geom.dx_fcn(p);

Disp = k_1';
for j=1:length(Disp(:,1))
    for i=1:length(Disp(:,1))
        Disp(i,j) = Disp(i,j) + k_2(1,:,j,i)*q + q'*k_3(:,:,j,i)*q;
    end
end
Disp = Disp*dx;

end