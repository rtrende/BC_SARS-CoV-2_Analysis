function [cscsim] = CSCsim(data)
    data = data ./ sum(data);
    cscsim = zeros(size(data,2));

    for i = 1:size(data,2)
        for j = i:size(data,2)
            cscsim(i,j) =  1 - real(2/pi * sqrt(2*(1-sum(sqrt(data(:,i).*data(:,j))))));
            cscsim(j,i) = cscsim(i,j);
        end
    end
end