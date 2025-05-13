function M = massMatr(q,qt,p0,p)

[k_1, k_2, k_3] = project.inertia.combMassMatr(p);

% M = k_1;
% 
% for j=1:length(M(:,1))
%     for i=1:length(M(:,1))
%         M(i,j) = M(i,j) + k_2(1,:,i,j)*q + q'*k_3(:,:,i,j)*q;
%     end
% end

M = k_1 + squeeze(pagemtimes(k_2,q)) +...
    squeeze(pagemtimes(q',pagemtimes(k_3,q)));


end