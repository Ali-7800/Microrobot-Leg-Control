function tip_displacement = fcn_FK(dv1,dv2)
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
    N_beams = params(8);
    l = lengths(1);
    BC = lengths(2); 
    Bb = lengths(3); 
    h = lengths(6);
    t = lengths(7);
    L = lengths(8);
    b = lengths(9);
    r = lengths(5)/2;

    syms N1 N2 P1 P2 F1 F2 alph C delta1 delta2 gamma1 gamma2
    
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

    %calculate V1 and V2 from input variables
    V1 = (dv1*L/t)*sqrt(a/(rho*Kp));
    V2 = (dv2*L/t)*sqrt(a/(rho*Kp));
    
    
    %forward kinematic system of equations (No Contact)
    eqn1 = M_B1 - M_A1 - F1*l + P1*gamma1 == 0;
    eqn2 = M_B2 - M_A2 - F2*l + P2*gamma2 == 0;
    eqn3 = M_A2 - M_A1 - F1*(l+Bb-delta1) + F2*(l+Bb-delta2) == 0;
    eqn4 = P1+F2 == 0;
    eqn5 = P2+F1 == 0;
    eqn6 = delta1 == h*((2*(Le(N1))^2/La(N1))*(sqrt(1-(La(N1)/(4*(Le(N1))^2))*((N1-V1^2)/(12*Q^2)))-1));
    eqn7 = P1/N_beams == E*b*t^3*(4*h*N1+(delta1/Le(N1)))/(12*L^3);
    eqn8 = delta2 == h*((2*(Le(N2))^2/La(N2))*(sqrt(1-(La(N2)/(4*(Le(N2))^2))*((N2-V2^2)/(12*Q^2)))-1));
    eqn9 = P2/N_beams == E*b*t^3*(4*h*N2+(delta2/Le(N2)))/(12*L^3);
    eqn10 = 2*Bb*(1-cos(alph)) == delta1+delta2-gamma1-gamma2;
    eqn11 = 2*Bb*sin(alph) == delta1-delta2+gamma1-gamma2;

    
    S = vpasolve([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10,eqn11],[N1 N2 P1 P2 F1 F2 alph delta1 delta2 gamma1 gamma2],[-16*pi^2,16*pi^2;-16*pi^2,16*pi^2;-0.05,0.05;-0.05,0.05;-0.05,0.05;-0.05,0.05;-0.5,0.5;-200,200;-200,200;-200,200;-200,200]);
    
    %solve for dx_c and dy_c
    alph1 = S.alph;
    delta1 = S.delta1;
    delta2 = S.delta2;
    dx_c = (sqrt(2)*(delta1 - delta2 + sqrt(2)*BC*sin(alph1)))/2;
    dy_c = (sqrt(2)*(delta1 - 2*Bb + delta2 - sqrt(2)*BC + sqrt(2)*Bb*cos(pi/4-alph1) + sqrt(2)*BC*cos(alph1) + sqrt(2)*Bb*sin(pi/4-alph1)))/2;
    
    contact = 0;
    %%%%%%check if there is contact (comment this section for without contact model)
    contacts = fcn_contacts(dx_c,dy_c);
    alpha_min = contacts(2,1);
    alpha_max = contacts(2,2);
    
    if alph1<=alpha_min
        alph1 = alpha_min;
        contact = contacts(1,1);
    elseif alph1>=alpha_max
        alph1 = alpha_max;
        contact = contacts(1,2);
    else
        contact = 0;
    end
    %%%%%%
    
    
    if contact == 0
        tip_displacement = double([dx_c,dy_c]);
    else
        syms x y N1 N2 P1 P2 F1 F2 alph C delta1 delta2 gamma1 gamma2
        theta = abs(alph);
        if contact == 1
            Cx = C*cos(theta);
            Cy = C*sin(theta);
            mt = -tan(pi/2+alph);
        elseif contact == 2
            Cx = -C*cos(theta);
            Cy = C*sin(theta);
            mt = tan(pi/2-alph);
        elseif contact == 3
            Cx = -C*cos(theta);
            Cy = -C*sin(theta);
            mt = -tan(pi/2+alph);
        elseif contact == 4
            Cx = C*cos(theta);
            Cy = -C*sin(theta);
            mt = tan(pi/2-alph);
        end

        %define contact point and moment
        Vc=[sqrt(2)*(Cx+Cy)/2;sqrt(2)*(-Cx+Cy)/2;0];
        Rc = [coords(contact,1);coords(contact,2);0];
        Vxy = [Cx;Cy;0];
        Mc = cross(Rc,Vxy);
        
        
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
        
    
        x1 = coords(1,1);
        x2 = coords(2,1);
        y1 = coords(1,2);
        y4 = coords(4,2);
        
        %solve for P1 and P2
        dx_c = (sqrt(2)*(delta1 - delta2 + sqrt(2)*BC*sin(alph)))/2;
        dy_c = (sqrt(2)*(delta1 - 2*Bb + delta2 - sqrt(2)*BC + sqrt(2)*Bb*cos(pi/4-alph) + sqrt(2)*BC*cos(alph) + sqrt(2)*Bb*sin(pi/4-alph)))/2;

        %forward kinematic system of equations (Contact)
        eqn1 = M_B1 - M_A1 - F1*l + P1*gamma1 == 0;
        eqn2 = M_B2 - M_A2 - F2*l + P2*gamma2 == 0;
        eqn3 = M_A2 - M_A1 - F1*(l+Bb-delta1) + F2*(l+Bb-delta2) + Mc(3) == 0;
        eqn4 = P1+F2+Vc(1) == 0;
        eqn5 = P2+F1+Vc(2) == 0;
        eqn6 = delta1 == h*((2*(Le(N1))^2/La(N1))*(sqrt(1-(La(N1)/(4*(Le(N1))^2))*((N1-V1^2)/(12*Q^2)))-1));
        eqn7 = P1/N_beams == E*b*t^3*(4*h*N1+(delta1/Le(N1)))/(12*L^3);
        eqn8 = delta2 == h*((2*(Le(N2))^2/La(N2))*(sqrt(1-(La(N2)/(4*(Le(N2))^2))*((N2-V2^2)/(12*Q^2)))-1));
        eqn9 = P2/N_beams == E*b*t^3*(4*h*N2+(delta2/Le(N2)))/(12*L^3);
        eqn10 = y == mt*(x-dx_c)+dy_c+BC;
        eqn11 = (x-coords(contact,1))^2+(y-coords(contact,2))^2 == r^2;
        eqn12 = y == -((x-coords(contact,1))/mt)+coords(contact,2);
        eqn13 = 2*Bb*(1-cos(alph)) == delta1+delta2-gamma1-gamma2;
        eqn14 = 2*Bb*sin(alph) == delta1-delta2+gamma1-gamma2;

        S = vpasolve([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10,eqn11,eqn12,eqn13,eqn14],[N1 N2 P1 P2 F1 F2 C delta1 delta2 gamma1 gamma2 alph x y],[-16*pi^2,16*pi^2;-16*pi^2,16*pi^2;-0.05,0.05;-0.05,0.05;-0.05,0.05;-0.05,0.05;-0.05,0.05;-200,200;-200,200;-300,300;-300,300;-0.2,0.2;x1,x2;y4-r,y1+r]);

        %solve for dx_c and dy_c
        delta1 = S.delta1;
        delta2 = S.delta2;
        alph = S.alph;
        dx_c = (sqrt(2)*(delta1 - delta2 + sqrt(2)*BC*sin(alph)))/2;
        dy_c = (sqrt(2)*(delta1 - 2*Bb + delta2 - sqrt(2)*BC + sqrt(2)*Bb*cos(pi/4-alph) + sqrt(2)*BC*cos(alph) + sqrt(2)*Bb*sin(pi/4-alph)))/2;
        tip_displacement = double([dx_c,dy_c]);
    end
    

end