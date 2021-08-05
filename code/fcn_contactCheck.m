function contact = fcn_contactCheck(dx_c,dy_c)  
    %this function checks if the leg makes contact
    p = get_params();
    params = p.params;
    lengths = p.lengths;
    E = params(1);
    I = params(2);
    l = lengths(1);
    BC = lengths(2); 
    Bb = lengths(3); 



    syms alp P1 P2 F1 F2 C

    %find deltas
    dx_b = dx_c-BC*sin(alp);
    dy_b = dy_c+BC*(1-cos(alp));
    dx_b1=dx_b-Bb*(cos(pi/4-alp)-cos(pi/4));
    dy_b1=dy_b-Bb*(sin(pi/4-alp)-sin(pi/4));
    dx_b2=dx_b+Bb*(sin(pi/4-alp)-sin(pi/4));
    dy_b2=dy_b-Bb*(cos(pi/4-alp)-cos(pi/4));
    gamma1 = (-dx_b1+dy_b1)*sqrt(2)/2;
    gamma2 = (dx_b2+dy_b2)*sqrt(2)/2;
    delta1 = (dx_b1+dy_b1)*sqrt(2)/2;
    delta2 = (-dx_b2+dy_b2)*sqrt(2)/2;
    
    
    
    %find coefficients
    %beam A1B1
    mu1 = sqrt(P1/(E*I)); %(1/mm)
    beta1 = -alp;
    c11 = (csc(l*mu1/2))^2*((-l*beta1+gamma1)*mu1-(gamma1*mu1*cos(l*mu1))+beta1*sin(mu1*l))/(2*mu1*(-2+l*mu1*cot(l*mu1/2)));
    c12 = (beta1-beta1*cos(l*mu1)-gamma1*mu1*sin(l*mu1))/(mu1*(-2+2*cos(l*mu1)+l*mu1*sin(l*mu1)));
    
    %beam A2B2
    mu2 = sqrt(P2/(E*I)); %(1/mm)
    beta2 = alp;
    c21 = (csc(l*mu2/2))^2*((-l*beta2+gamma2)*mu2-(gamma2*mu2*cos(l*mu2))+beta2*sin(mu2*l))/(2*mu2*(-2+l*mu2*cot(l*mu2/2)));
    c22 = (beta2-beta2*cos(l*mu2)-gamma2*mu2*sin(l*mu2))/(mu2*(-2+2*cos(l*mu2)+l*mu2*sin(l*mu2)));    
    
    %find moments
    M_A1 = -P1*c11;
    M_A2 = -P2*c21;
    M_B1 = -P1*(c11*cos(mu1*l)+c12*sin(mu1*l));
    M_B2 = -P2*(c21*cos(mu2*l)+c22*sin(mu2*l));
    
    %solve for P1 and P2
    eqn1 = M_B1 - M_A1 - F1*l + P1*gamma1 == 0;
    eqn2 = M_B2 - M_A2 - F2*l + P2*gamma2 == 0;
    eqn3 = M_A2 - M_A1 - F1*(l+Bb-delta1) + F2*(l+Bb-delta2) == 0;
    eqn4 = P1+F2 == 0;
    eqn5 = P2+F1 == 0;


    S = vpasolve([eqn1,eqn2,eqn3,eqn4,eqn5],[P1,P2,alp,F1,F2], [-0.1 0.1;-0.1 0.1;-0.1 0.1;-0.1 0.1;-0.1 0.1]);
    alph = S.alp;
        
    %Check Contact
    contacts = fcn_contacts(dx_c,dy_c);
    alpha_min = contacts(2,1);
    alpha_max = contacts(2,2);
    
    if alph<=alpha_min
        alph = alpha_min;
        contact = contacts(1,1);
    elseif alph>=alpha_max
        contact = contacts(1,2);
        alph = alpha_max;
    else
        contact = 0;
    end
end