function writeFunction(Matr, odr, fPath, fcnName, parIdx)

file = fullfile(fPath, [fcnName, '.m']);
fid = fopen(file, 'w');

fprintf(fid, 'function K = %s(p)\n', fcnName);
fprintf(fid, '%% Auto-generated constant matrix\n');

switch odr

    case 0
        fprintf(fid, 'K = zeros(%i,1);\n', length(Matr));
        fprintf(fid, 'K = [\n');
        K_numeric = Matr;  % or however you computed the matrix

        for i = 1:size(K_numeric, 1)
            fprintf(fid, '  ');
            fprintf(fid, '%e;\n', K_numeric(i, 1));
        end
        fprintf(fid, '];\n');

    case 1
        fprintf(fid, 'K = zeros(%i,%i);\n', size(Matr));
        fprintf(fid, 'K = [\n');
        K_numeric = Matr;  % or however you computed the matrix

        for i = 1:size(K_numeric, 1)
            fprintf(fid, '  ');
            fprintf(fid, '%e, ', K_numeric(i, 1:end-1));
            fprintf(fid, '%e;\n', K_numeric(i, end));
        end
        fprintf(fid, '];\n');

    case 2
        sz = size(Matr);
        fprintf(fid, 'K = zeros(%i,%i,%i,%i);\n', size(Matr));
        for I=1:sz(3)
            for J=1:sz(4)
                fprintf(fid, 'K(:,:,%i,%i) = [\n', I,J);
                K_numeric = Matr(:,:,I,J);  % or however you computed the matrix

                for i = 1:size(K_numeric, 1)
                    fprintf(fid, '  ');
                    fprintf(fid, '%e, ', K_numeric(i, 1:end-1));
                    fprintf(fid, '%e;\n', K_numeric(i, end));
                end
                fprintf(fid, '];\n');
            end
        end

    case 3
        fprintf(fid, 'K = zeros(%i,%i,%i,%i);\n', size(Matr));
        sz = size(Matr);
        for I=1:sz(3)
            for J=1:sz(4)
                fprintf(fid, 'K(:,:,%i,%i) = [\n', I,J);
                K_numeric = Matr(:,:,I,J);  % or however you computed the matrix

                for i = 1:size(K_numeric, 1)
                    fprintf(fid, '  ');
                    fprintf(fid, '%e, ', K_numeric(i, 1:end-1));
                    fprintf(fid, '%e;\n', K_numeric(i, end));
                end
                fprintf(fid, '];\n');
            end
        end

end

if ~isempty(parIdx)
    fprintf(fid, 'K=p(%i)*K;\n', parIdx);
end
fclose(fid);

end