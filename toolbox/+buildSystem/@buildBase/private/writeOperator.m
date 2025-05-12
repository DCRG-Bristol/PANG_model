function writeOperator(type, fPath, fcnName, inpt)

%type: 'vector' or 'matrix': if vector outouts in 'RHS form': e.g: 

file = fullfile(fPath, [fcnName, '.m']);
fid = fopen(file, 'w');

switch type
    case 'vector'

        fprintf(fid, 'function [k_0, k_1, k_2, k_3] = %s(p)\n', fcnName);
        fprintf(fid, '%% Auto-generated constant matrix\n');

        %% collect terms..

        for odr=0:3
            fprintf(fid, 'k_%i = ', odr);
            for itm=1:length(inpt)-1
                fprintf(fid, '%s.odr_%i.%s(p)+', replace(fPath,{'\','+'},{'.',''}), odr, inpt{itm}.fcnName);
            end
            fprintf(fid, '%s.odr_%i.%s(p);\n', replace(fPath,{'\','+'},{'.',''}), odr, inpt{end}.fcnName);
        end

        fclose(fid);


    case 'matrix'

    fprintf(fid, 'function [k_1, k_2, k_3] = %s(p)\n', fcnName);
    fprintf(fid, '%% Auto-generated constant matrix\n');

    %% collect terms..

    for odr=1:3
        fprintf(fid, 'k_%i = ', odr);
        for itm=1:length(inpt)-1
            fprintf(fid, '%s.odr_%i.%s(p)+', replace(fPath,{'\','+'},{'.',''}), odr, inpt{itm}.fcnName);
        end
        fprintf(fid, '%s.odr_%i.%s(p);\n', replace(fPath,{'\','+'},{'.',''}), odr, inpt{end}.fcnName);
    end

    fclose(fid);

end

end