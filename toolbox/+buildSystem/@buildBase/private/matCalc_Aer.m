function matr = matCalc_Aer(argin,N,xGrid)
odr = argin.odr;
scalF = argin.scalF;
vctr = argin.vctr;

if odr==1
    for J=1:length(xGrid)
        x = xGrid(J);
        matr(J,:) = scalF(x)*(vctr{1}(x,1)');
    end
end

if odr==2
    for I=1:N
        for J=1:length(xGrid)
            x = xGrid(J);
            Ipick = [zeros(1,I-1),1,zeros(1,N-I)];
            pickFcn = scalF(x)*(Ipick*vctr{1}(x,1));
            matr(1,:,J,I) = pickFcn*(vctr{2}(x,I)');
        end
    end
end

if odr==3
    for I=1:N
        for J=1:length(xGrid)
            x = xGrid(J);
            Ipick = [zeros(1,I-1),1,zeros(1,N-I)];
            pickFcn = scalF(x)*(Ipick*vctr{1}(x,1));
            matr(:,:,J,I) = pickFcn*(vctr{2}(x,I).*vctr{3}(x,I)');
        end
    end
end

end