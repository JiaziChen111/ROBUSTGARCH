function [fun] = ROBUSTGARCHloss(theta,r,sigma2R)
n = length(r);
k = 3;
h(1) = sigma2R;
J(1) = r(1)/sqrt(h(1));
for i=2:n
    if abs(J(i-1))<k
        h(i) = sigma2R*(1-theta(1)-theta(2))+ theta(1)*r(i-1)^2+ theta(2)*h(i-1);
    else
       h(i)= sigma2R*(1-theta(1)-theta(2))+ theta(1)*1.005018*h(i-1)+ theta(2)*h(i-1); 
    end
    J(i) = r(i)/sqrt(h(i));
end

for j=1:n
    if r(j) == 0
        aux(j) = r(j) + 0.00001;
    else
        aux(j) = r(j);
    end
end
y = log((aux.^2)./h);
fun = mean(-y + 0.8260*5*log(1+exp(y)./2));
end
