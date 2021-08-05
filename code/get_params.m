function p = get_params()     
     
     unit = 1e6; %1 = m, 1e3 = mm, 1e6 = micron

     %lengths
     l = 1.555e-3*unit; %beam length (length between A1 and B1 or A2 and B2) (m)
     BC = 1.78887e-3*unit; %distance between B and C (m)
     Bb = 0.2e-3*unit; %distance between B and B1/B2 (m)
     BD = 0.401130e-3*unit; %distance between B and D (m)
     Dd = 0.1e-3*unit; %leg width (m)
     h = 0.017e-3*unit; %h for V-shape actuator (m)
     t = 0.017e-3*unit; %in plane thickness for V-shape actuator (m)
     b = 0.025e-3*unit; %depth for for V-shape actuator (m)
     L = 3e-3*unit; %V-shape actuator span (m)
     w = 0.015e-3*unit; %beam width (m)
     
     %coordinates
     d1_x = -0.055e-3*unit; %Contact Point 1 (x)
     d2_x = 0.055e-3*unit; %Contact Point 2 (x)
     d_y1 = -0.296130e-3*unit; %Contact Point Upper (y)
     d_y2 = -0.301130e-3*unit; %Contact Point Lower (y)
     
     %parameters for beam simulation
     E = 169e9/unit^2; %Young's Modulus (Pa = N/m^2)
     I = (b)*(w)^3/12; %moment of inertia (m^4)
     Q = h/t;
     a = 3.4e-6; %thermal expansion coefficient (1/C)
     Kp = 149; %Conduction Coefficient 
     rho = 2.258478955955251e-04; %0.0017; %Electrical Resistivity (Ohm*m)
     dv_max = 15; %maximum voltage
     N_beams = 13;
     
     p.coords = [d1_x,d_y1;d2_x,d_y1;d2_x,d_y2;d1_x,d_y2];
     p.params = [E,I,Q,a,Kp,rho,dv_max,N_beams];
     p.lengths = [l,BC,Bb,BD,Dd,h,t,L,b];


end