function [vi] = gridGARCH(r, sigma2R)
alfa1min = 0.005;
alfa1max = 0.2;
beta1min = 0.65;
beta1max = 0.98;
nalfa1 = 5;
nbeta1 = 5;
ml = 100000000;
lmalfa1 = (alfa1max-alfa1min)/nalfa1;
lmbeta1 = (beta1max-beta1min)/nbeta1;

for nj = 0:nalfa1
    for nk = 0:nbeta1
      alfa1 = alfa1min+nj*lmalfa1;
      beta1 = beta1min+nk*lmbeta1; 
      if (alfa1 + beta1 < 0.999)
          coeff(1) = alfa1;
          coeff(2) = beta1;
          nml = ROBUSTGARCHloss(coeff, r, sigma2R);
          if nml < ml
              vi(1) = coeff(1);
              vi(2) = coeff(2);
          end
      end
    end
end
end
          
      
