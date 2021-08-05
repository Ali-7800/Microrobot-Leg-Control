function contacts = fcn_contacts(dx_c,dy_c)
    %determines which contact point the leg is touching and the range of alphas
    p = get_params();
    coords = p.coords;
    lengths = p.lengths;
    BC = lengths(2);
    r = lengths(5)/2;

    %contact point coordinates
    c1 = coords(1,:);
    c2 = coords(2,:);
    c3 = coords(3,:);
    c4 = coords(4,:);
    
    x1 = c1(1);
    x2 = c2(1);
    x3 = c3(1);
    x4 = c4(1);
    
    y1 = c1(2);
    y2 = c2(2);
    y3 = c3(2);
    y4 = c4(2);

    xc = dx_c;
    yc = dy_c+BC;
    
    syms x y mt
    eqn1 = y == mt*(x-xc)+yc;
    eqn2 = x<x2;
    eqn3 = x>x1;
    
    if xc>x2-r
        %calculate alpha min
        eqn4 = (x-x2)^2+(y-y2)^2 == r^2;
        eqn5 = y == -((x-x2)/mt)+y2;
        S = solve([eqn1,eqn2,eqn3,eqn4,eqn5],[x,y,mt],'Real',true);
        M_min = S.mt;
        alpha_min = pi/2-atan(M_min);
        
        %calculate alpha max
        eqn4 = (x-x4)^2+(y-y4)^2 == r^2;
        eqn5 = y == -((x-x4)/mt)+y4;
        S = solve([eqn1,eqn2,eqn3,eqn4,eqn5],[x,y,mt],'Real',true);
        M_max = S.mt;
        alpha_max = pi/2-atan(M_max);
        
        contacts = [2,4;alpha_min,alpha_max];

    elseif xc<x1+r
        %calculate alpha min
        eqn4 = (x-x3)^2+(y-y3)^2 == r^2;
        eqn5 = y == -((x-x3)/mt)+y3;
        S = solve([eqn1,eqn2,eqn3,eqn4,eqn5],[x,y,mt],'Real',true);
        M_min = S.mt;
        alpha_min = -pi/2-atan(M_min);
        
        %calculate alpha max
        eqn4 = (x-x1)^2+(y-y1)^2 == r^2;
        eqn5 = y == -((x-x1)/mt)+y1;
        S = solve([eqn1,eqn2,eqn3,eqn4,eqn5],[x,y,mt],'Real',true);
        M_max = S.mt;
        alpha_max = -pi/2-atan(M_max);
        
        contacts = [3,1;alpha_min,alpha_max];

    elseif xc==x1+r
        %calculate alpha min
        eqn4 = (x-x3)^2+(y-y3)^2 == r^2;
        eqn5 = y == -((x-x3)/mt)+y3;
        S = solve([eqn1,eqn2,eqn3,eqn4,eqn5],[x,y,mt],'Real',true);
        M_min = S.mt;
        alpha_min = -pi/2-atan(M_min);
        
        %calculate alpha max
        alpha_max = 0;
        contacts = [3,1;alpha_min,alpha_max];
        
    elseif xc==x2-r
        %calculate alpha min
        alpha_min = 0;
        
        %calculate alpha max
        eqn4 = (x-x4)^2+(y-y4)^2 == r^2;
        eqn5 = y == -((x-x4)/mt)+y4;
        S = solve([eqn1,eqn2,eqn3,eqn4,eqn5],[x,y,mt],'Real',true);
        M_max = S.mt;
        alpha_max = pi/2-atan(M_max);
        
        contacts = [2,4;alpha_min,alpha_max];
        
    else
        %calculate alpha min
        eqn4 = (x-x4)^2+(y-y4)^2 == r^2;
        eqn5 = y == -((x-x4)/mt)+y4;
        S = solve([eqn1,eqn2,eqn3,eqn4,eqn5],[x,y,mt],'Real',true);
        M_max = S.mt;
        alpha_max = pi/2-atan(M_max);
        
        %calculate alpha max
        eqn4 = (x-x3)^2+(y-y3)^2 == r^2;
        eqn5 = y == -((x-x3)/mt)+y3;
        S = solve([eqn1,eqn2,eqn3,eqn4,eqn5],[x,y,mt],'Real',true);
        M_min = S.mt;
        alpha_min = -pi/2-atan(M_min);
        
        contacts = [3,4;alpha_min,alpha_max];

    end
    

end