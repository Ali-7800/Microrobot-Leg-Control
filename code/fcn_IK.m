function IK = fcn_IK(dx_c,dy_c)
    %define parameters
    p = get_params();
    params = p.params;
    coords = p.coords;
    lengths = p.lengths;
    E = params(1);
    I = params(2);
    Q = params(3);
    a = params(4);
    Kp = params(5);
    rho = params(6);
    dv_max = params(7);
    N_beams = params(8);
    l = lengths(1);
    BC = lengths(2); 
    Bb = lengths(3); 
    h = lengths(6);
    t = lengths(7);
    L = lengths(8);
    b = lengths(9);
    K = 1e5;
    epsilon = 1e-3;
    V_max = dv_max*L*sqrt(a/(rho*Kp))/t;



    syms alp P1 P2 F1 F2 C

    %calculate deltas
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
    
    
    
    %calculate coefficients for beams
    %beam A1B1
    mu1 = sqrt(P1/(E*I)); 
    beta1 = -alp;
    c11 = (csc(l*mu1/2))^2*((-l*beta1+gamma1)*mu1-(gamma1*mu1*cos(l*mu1))+beta1*sin(mu1*l))/(2*mu1*(-2+l*mu1*cot(l*mu1/2)));
    c12 = (beta1-beta1*cos(l*mu1)-gamma1*mu1*sin(l*mu1))/(mu1*(-2+2*cos(l*mu1)+l*mu1*sin(l*mu1)));
    
    %beam A2B2
    mu2 = sqrt(P2/(E*I)); 
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
        
    contact = 0;
    %%%%%%check if there is contact (comment this section for without contact model)
    contacts = fcn_contacts(dx_c,dy_c);
    alpha_min = contacts(2,1);
    alpha_max = contacts(2,2);
    
    if alph<=alpha_min
        alph = alpha_min;
        contact = contacts(1,1);
    elseif alph>=alpha_max
        alph = alpha_max;
        contact = contacts(1,2);
    else
        contact = 0;
    end
    %%%%%%
    
    
    delta1 = subs(delta1,alph);
    delta2 = subs(delta2,alph);
    gamma1 = subs(gamma1,alph);
    gamma2 = subs(gamma2,alph);
    
    %recalculate coefficients for beams
    %beam A1B1
    mu1 = sqrt(P1/(E*I)); 
    beta1 = -alph;
    c11 = (csc(l*mu1/2))^2*((-l*beta1+gamma1)*mu1-(gamma1*mu1*cos(l*mu1))+beta1*sin(mu1*l))/(2*mu1*(-2+l*mu1*cot(l*mu1/2)));
    c12 = (beta1-beta1*cos(l*mu1)-gamma1*mu1*sin(l*mu1))/(mu1*(-2+2*cos(l*mu1)+l*mu1*sin(l*mu1)));
    
    %beam A2B2
    mu2 = sqrt(P2/(E*I)); 
    beta2 = alph;
    c21 = (csc(l*mu2/2))^2*((-l*beta2+gamma2)*mu2-(gamma2*mu2*cos(l*mu2))+beta2*sin(mu2*l))/(2*mu2*(-2+l*mu2*cot(l*mu2/2)));
    c22 = (beta2-beta2*cos(l*mu2)-gamma2*mu2*sin(l*mu2))/(mu2*(-2+2*cos(l*mu2)+l*mu2*sin(l*mu2)));    
    
    %find moments
    M_A1 = -P1*c11;
    M_A2 = -P2*c21;
    M_B1 = -P1*(c11*cos(mu1*l)+c12*sin(mu1*l));
    M_B2 = -P2*(c21*cos(mu2*l)+c22*sin(mu2*l));
    
    
    %apply the contact model if leg is making contact
    if  contact == 0
        P1 = S.P1;
        P2 = S.P2;
        %"No Contact"
    else
        %"contact"
        theta = abs(alph);
        if contact == 1
            Cx = C*cos(theta);
            Cy = C*sin(theta);
        elseif contact == 2
            Cx = -C*cos(theta);
            Cy = C*sin(theta);
        elseif contact == 3
            Cx = -C*cos(theta);
            Cy = -C*sin(theta);
        elseif contact == 4
            Cx = C*cos(theta);
            Cy = -C*sin(theta);
        end
        
        %define contact moment and force
        Vc=[sqrt(2)*(Cx+Cy)/2;sqrt(2)*(-Cx+Cy)/2;0];
        Vxy = [Cx;Cy;0];
        Rc = [coords(contact,1);coords(contact,2);0];
        Mc = cross(Rc,Vxy);

        
        %solve for P1 and P2
        eqn1 = M_B1 - M_A1 - F1*l + P1*gamma1 == 0;
        eqn2 = M_B2 - M_A2 - F2*l + P2*gamma2 == 0;
        eqn3 = M_A2 - M_A1 - F1*(l+Bb-delta1) + F2*(l+Bb-delta2) + Mc(3) == 0;
        eqn4 = P1+F2+Vc(1) == 0;
        eqn5 = P2+F1+Vc(2) == 0;
        
        S = vpasolve([eqn1,eqn2,eqn3,eqn4,eqn5],[P1,P2,F1,F2,C], [-0.1 0.1;-0.1 0.1;-0.1,0.1;-0.1 0.1;-0.1 0.1]);
        P1 = S.P1;
        P2 = S.P2;
    end
        
    %parameters
    syms N2 V
    delta1 = double(delta1);
    delta2 = double(delta2);
    P1 = double(P1);
    P2 = double(P2);
    

    %solve for N1 and V1
    eqn = P1/N_beams == E*b*t^3*(4*h*N2+(delta1/Le(N2)))/(12*L^3);
    N = vpasolve(eqn,N2, [-16*pi^2 16*pi^2]);
    if size(N,1)==0
        V1 = NaN;
    else
        eqn1 = delta1 == h*((2*(Le(N))^2/La(N))*(sqrt(1-(La(N)/(4*(Le(N))^2))*((N-V^2)/(12*Q^2)))-1));
        V1 = vpasolve(eqn1,V,[0 V_max]);
    end
    if size(V1,1)==0
        V = linspace(0,V_max,K);
        error = abs(delta1-h*((2*(Le(N))^2/La(N)).*(sqrt(1-((La(N)/(4*(Le(N))^2)).*((N-V.^2)./(12*Q^2))))-1)));
        [M,I] = min(error);
        if M<epsilon
            V1 = V_max*I/K;
        else
            V1 = NaN;
        end
     end
    
    

    %solve for N2 and V2
    syms N2 V
    eqn = P2/N_beams == E*b*t^3*(4*h*N2+(delta2/Le(N2)))/(12*L^3);
    N = vpasolve(eqn,N2, [-16*pi^2 16*pi^2]);
    if size(N,1)==0
        V2 = NaN;
    else
        eqn1 = delta2 == h*((2*(Le(N))^2/La(N))*(sqrt(1-(La(N)/(4*(Le(N))^2))*((N-V^2)/(12*Q^2)))-1));
        V2 = vpasolve(eqn1,V,[0 V_max]);
    end
    if size(V2,1)==0
        V = linspace(0,V_max,K);
        error = abs(delta2-h*((2*(Le(N))^2/La(N)).*(sqrt(1-((La(N)/(4*(Le(N))^2)).*((N-V.^2)./(12*Q^2))))-1)));
        [M,I] = min(error);
        if M<epsilon
            V2 = V_max*I/K;
        else
            V2 = NaN;
        end
    end
    
        
    %calculate dv1 and dv2 (delta V1 and delta V2)
    dv1 = (V1*t/L)*sqrt(rho*Kp/a);
    dv2 = (V2*t/L)*sqrt(rho*Kp/a);
    if size(dv1,1)==0
       dv1 = NaN;
    end
    if size(dv2,1)==0
       dv2 = NaN;
    end

    IK = double([dv1,dv2]);
end