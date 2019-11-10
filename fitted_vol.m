function [vol] = fitted_vol(theta,r)
n = length(r)+1;
k = 3;
h(1) = theta(1)/(1-theta(2)-theta(3));
for t=2:n
    if abs(r(t-1)/sqrt(h(t-1))) < k
        h(t)= theta(1)+ theta(2)*r(t-1)^2+ theta(3)*h(t-1);
    else
        h(t)= theta(1)+ theta(2)*1.005018*h(t-1)+ theta(3)*h(t-1);
    end
end
vol = sqrt(h);
end