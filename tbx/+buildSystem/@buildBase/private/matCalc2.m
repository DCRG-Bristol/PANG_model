function matr = matCalc2(argin,N,L)
odr = argin.odr;
scalF = argin.scalF;
vctr = argin.vctr;

if odr==1
    matr = integral(@(x)(scalF(x)*(vctr{1}(x,1).*vctr{2}(x,1)')),...
        0, L, 'ArrayValued', true);
end

if odr==2
    for I=1:N
        for J=1:N
            Ipick = [zeros(1,I-1),1,zeros(1,N-I)];
            Jpick = [zeros(1,J-1),1,zeros(1,N-J)];
            pickFcn = @(x)(scalF(x)*(Jpick*vctr{1}(x,1)).*...
                (Ipick*vctr{2}(x,1))');
            matr(1,:,J,I) = integral(@(x)(pickFcn(x)*(vctr{3}(x,J)')), 0, L, 'ArrayValued', true);
        end
    end
end

if odr==3
    for I=1:N
        for J=1:N
            Ipick = [zeros(1,I-1),1,zeros(1,N-I)];
            Jpick = [zeros(1,J-1),1,zeros(1,N-J)];
            pickFcn = @(x)(scalF(x)*(Jpick*vctr{1}(x,1)).*...
                (Ipick*vctr{2}(x,1))');
            matr(:,:,J,I) = integral(@(x)(pickFcn(x)*...
                (vctr{3}(x,J).*vctr{4}(x,I)')), 0, L, 'ArrayValued', true);
        end
    end
end

if odr==4
    for I=1:N
        for J=1:N
            Ipick = [zeros(1,I-1),1,zeros(1,N-I)];
            Jpick = [zeros(1,J-1),1,zeros(1,N-J)];
            pickFcn_1 = @(x)(scalF(x)*(Jpick*vctr{1}(x,1)).*...
                (Ipick*vctr{2}(x,1))');
            for i=1:N
                for j=1:N
                    Ipick2 = [zeros(1,i-1),1,zeros(1,N-i)];
                    Jpick2 = [zeros(1,j-1),1,zeros(1,N-j)];
                    pickFcn_2 = @(x)((Jpick2*vctr{3}(x,1)).*...
                        (Ipick2*vctr{4}(x,1))');
                    matr(1,:,j,i,J,I) = integral(@(x)(pickFcn_1(x)*pickFcn_2(x)*...
                        (vctr{5}(x,j)')), 0, L, 'ArrayValued', true);
                end
            end
        end
    end
end

if odr==5
    for I=1:N
        for J=1:N
            Ipick = [zeros(1,I-1),1,zeros(1,N-I)];
            Jpick = [zeros(1,J-1),1,zeros(1,N-J)];
            pickFcn_1 = @(x)(scalF(x)*(Jpick*vctr{1}(x,1)).*...
                (Ipick*vctr{2}(x,1))');
            for i=1:N
                for j=1:N
                    Ipick2 = [zeros(1,i-1),1,zeros(1,N-i)];
                    Jpick2 = [zeros(1,j-1),1,zeros(1,N-j)];
                    pickFcn_2 = @(x)((Jpick2*vctr{3}(x,1)).*...
                        (Ipick2*vctr{4}(x,1))');
                    matr(:,:,j,i,J,I) = integral(@(x)(pickFcn_1(x)*pickFcn_2(x)*...
                        (vctr{5}(x,j).*vctr{6}(x,i)')), 0, L, 'ArrayValued', true);
                end
            end
        end
    end
end

end