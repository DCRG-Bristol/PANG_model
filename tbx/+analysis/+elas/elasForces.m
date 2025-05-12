function F = elasForces(q,qt,p0,p)


[k_1, k_2, k_3] = project.elas.combElasMatr(p);

K = k_1;

for j=1:length(K(:,1))
    for i=1:length(K(:,1))
        K(i,j) = K(i,j) + k_2(1,:,i,j)*q + q'*k_3(:,:,i,j)*q;
    end
end

F = -K*q;

end