function [k_0, k_1, k_2, k_3] = combGravMats_1(p)
% Auto-generated constant matrix
k_0 = project.grav.odr_0.baff_struct_dir_1(p)+project.grav.odr_0.fuel_dir_1(p);
k_1 = project.grav.odr_1.baff_struct_dir_1(p)+project.grav.odr_1.fuel_dir_1(p);
k_2 = project.grav.odr_2.baff_struct_dir_1(p)+project.grav.odr_2.fuel_dir_1(p);
k_3 = project.grav.odr_3.baff_struct_dir_1(p)+project.grav.odr_3.fuel_dir_1(p);
