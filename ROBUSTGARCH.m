function [theta] = ROBUSTGARCH(r)
AUX = medianB(r);
Med = AUX(:,1);
Mad = AUX(:,2);
I = ((r-Med).^2)./((1.486*Mad).^2) < 3.841459;
mu_R = sum(r.*I)/sum(I);
J = ((r-mu_R).^2)./((1.486*Mad).^2) < 3.841459;
sigma2R = 1.318 * sum(((r-mu_R).^2).*J)/sum(J);
ini = gridGARCH(r-mu_R, sigma2R);
estim = @(x) ROBUSTGARCHloss(x,r-mu_R,sigma2R);
param = fmincon(estim,ini,[1 1],[0.9999],[],[],[0.0001 0.0001],[1 1]);
theta = [sigma2R*(1-param(1)-param(2)) param];
end

