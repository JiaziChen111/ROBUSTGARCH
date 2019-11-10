function [MedMad] = medianB(r)
    k = 30;
    n = length(r);
    MedMad = zeros(n,2);
    for i = 1:n
        if i < k/2
            MedMad(i,1) = median(r(1:(k+1)));
            MedMad(i,2) = median(abs(r(1:(k+1)) - MedMad(i,1)));
        elseif i >n-k/2
            MedMad(i,1) = median(r((n-k):n));
            MedMad(i,2) = median(abs(r((n-k):n) - MedMad(i,1)));
        else
            MedMad(i,1) = median(r((i-k/2+1):(i+k/2)));
            MedMad(i,2) = median(abs(r((i-k/2+1):(i+k/2)) - MedMad(i,1)));
        end
    end
end