function [k_1, k_2, k_3] = combMassMatr(p)
% Auto-generated constant matrix
k_1 = project.inertia.odr_1.baff_struct(p)+project.inertia.odr_1.fuel(p);
k_2 = project.inertia.odr_2.baff_struct(p)+project.inertia.odr_2.fuel(p);
k_3 = project.inertia.odr_3.baff_struct(p)+project.inertia.odr_3.fuel(p);
