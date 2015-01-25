function xi = rat_poles(H, K)
% RAT_POLES computes the poles of the pencil (H, K).
  
  m  = size(H, 2);
  xi = zeros(1, m);
  
  j  = 1;
  while j <= m
    if j+2 <= m+1
      if K(j+2, j) == 0
        xi(j) = H(j+1, j)/K(j+1, j);
      else
        cxi = K(j+1:j+2, j:j+1)\H(j+1:j+2, j:j+1);
        xi(j)   = cxi(1, 1) + 1i*cxi(1, 2);
        xi(j+1) = conj(xi(j));
        j = j+1;
      end
    else
      xi(j) = H(j+1, j)/K(j+1, j);
    end
    j = j+1;
  end
  
end