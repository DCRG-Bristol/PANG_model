function F = gyrForces(q,qt,p0,p)

%if nargin==4
    [k_1, k_2, k_3] = project.inertia.combGyrMatr(p);
% else
%     [k_1, k_2, k_3] = feval([obj.name,'.project.inertia.combGyrMatr([]);
% end

M = k_1;

for j=1:length(M(:,1))
    for i=1:length(M(:,1))
        M(i,j) = M(i,j) + k_2(1,:,i,j)*qt + q'*k_3(:,:,i,j)*qt;
    end
end

F = -M*qt;

end