function [C,CC,B,u]=VARcov(x,covmatrix,T,k,qq,mm,H,nlagsimp)
A = zeros(qq*k,qq*k); B = zeros(qq*k,qq); xx = zeros(T-k,qq*k);
for j = 1:k
    for i = 1:k;
        A((j-1)*qq+1:qq*j,(i-1)*qq+1:qq*i) = covmatrix((mm-1)*qq+1:mm*qq, (mm-1)*qq+1:mm*qq,H + i - j);
    end
    B((j-1)*qq+1:qq*j,:) = covmatrix((mm-1)*qq+1:mm*qq, (mm-1)*qq+1:mm*qq,H - j);
    xx(:,(j-1)*qq+1:qq*j) =  x(k+1-j:T-j ,(mm-1)*qq+1:mm*qq);
end
C = inv(A)*B;                                                               % Aj(L)
u = x(k+1:T , (mm-1)*qq+1:mm*qq) - xx*C;                                    % Residual wj(t)
CC(:,:,1) = eye(qq); CC(:,:,2:k + 1) = - reshape(C',qq,qq,k);               % Reshaping Aj(L) s.t. it is invertible
B = InvPolMatrix(CC,nlagsimp);                                              % Building Gj(L)= inv(Aj(L))


function[med,mad] = medianB(r)
    k = 30
    n = length(r)
    for i = 1:n
        med(i) = ifelse(i<k/2,median(r(1:(k+1))),ifelse(i>n-k/2,median(r((n-k):n)),median(r((i-k/2):(i+k/2)))))
        mad(i) = ifelse(i<k/2,median(abs(r(1:(k+1)) - med(i))),ifelse(i>n-k/2,median(abs(r((n-k):n) - med(i))),median(abs(r((i-k/2):(i+k/2)) - med(i)))))
    end
end
        
        

function [theta] = RobGARCH(theta,sigma2,ret)
h[1] = sigma2



#' @noRd
#' @importFrom stats median
medianB = function(r){
  k = 30
  n = length(r)
  med = c()
  mad = c()
  for(i in 1:n){
    med[i] = ifelse(i<k/2,median(r[1:(k+1)]),ifelse(i>n-k/2,median(r[(n-k):n]),median(r[(i-k/2):(i+k/2)])))
    mad[i] = ifelse(i<k/2,median(abs(r[1:(k+1)] - med[i])),ifelse(i>n-k/2,median(abs(r[(n-k):n] - med[i])),median(abs(r[(i-k/2):(i+k/2)] - med[i]))))
  }
  return(list(med,mad))
}

#' @export
fitted_Vol = function(theta,r){
  n= length(r)+1
  h= c()
  k = 3
  h[1]= theta[1]/(1-theta[2]-theta[3])
  for (t in 2:n){
    if(abs(r[t-1]/sqrt(h[t-1]))<k){
      h[t]= theta[1]+ theta[2]*r[t-1]^2+ theta[3]*h[t-1]
    } else{
      h[t]= theta[1]+ theta[2]*1.005018*h[t-1]+ theta[3]*h[t-1]
    }
  }
  return(sqrt(h))
}

#' @export
#' @import Rcpp
#' @importFrom stats constrOptim
ROBUSTGARCH = function(y){
  AUX= medianB(y)
  Med = AUX[[1]]
  MAD = AUX[[2]]
  I = (y-Med)^2/(1.486*MAD)^2<= 3.841459
  mu_R = sum(y*I)/sum(I)
  J = (y-mu_R)^2/(1.486*MAD)^2<=3.841459
  sigma2R = 1.318 * sum((y-mu_R)^2*J)/sum(J)
  ini = grid_RCPP(y-mu_R, sigma2R)
  ra <- matrix(c(1,0,0,1,-1,-1),ncol=2,byrow=TRUE)
  rb <- c(0.00001, 0.00001,-0.9999)
  param <- constrOptim(theta = ini, f = ROBUSTGARCHloss_RCPP, grad = NULL, ui = ra, ci = rb, sigma2 = sigma2R, r = y, outer.iterations = 400, outer.eps = 1e-07)$par
  param = c(sigma2R*(1-param[1]-param[2]),param[1],param[2])
}

