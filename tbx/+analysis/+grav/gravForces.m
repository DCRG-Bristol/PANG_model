function F = gravForces(q,qt,p0,p)

%p0 is {'U', 'alpha0', 'alpha', 'g'}; %hardcoded essential parameters...

[k1_0,k1_1, k1_2, k1_3] = project.grav.combGravMats_1(p);
[k2_0,k2_1, k2_2, k2_3] = project.grav.combGravMats_2(p);

k_0 = -cos(p0(2)).*k1_0 + -sin(p0(2)).*k2_0;
k_1 = -cos(p0(2)).*k1_1 + -sin(p0(2)).*k2_1;
k_2 = -cos(p0(2)).*k1_2 + -sin(p0(2)).*k2_2;
k_3 = -cos(p0(2)).*k1_3 + -sin(p0(2)).*k2_3;

K = k_1;

for j=1:length(K(:,1))
    for i=1:length(K(:,1))
        K(i,j) = K(i,j) + k_2(1,:,i,j)*q + q'*k_3(:,:,i,j)*q;
    end
end

F = p0(4)*(K*q+k_0);

end