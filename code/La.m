function Lambda_a = La(N2)
    N = sqrt(N2);
    Lambda_a = 3*(1-(tan(N/4)/(N/4))+(tan(N/4)^2/3))/(16*N^4);
end