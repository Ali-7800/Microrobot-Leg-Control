Voltages = zeros([500 2]);
%% Format
% for i = 1:number_of_points
%     IK = fcn_IK(x(i),y(i));
%     Voltages(i,1) = double(IK(1,1)); %this is your V1 trajectory
%     Voltages(i,2) = double(IK(1,2)); %this is your V1 trajectory
% end
%% circle
for i = 1:500
    IK = fcn_IK(15*cos(2*pi*i/500-pi/2),15*sin(2*pi*i/500-pi/2)+20);
    Voltages(i,1) = double(IK(1,1));
    Voltages(i,2) = double(IK(1,2));
end
%% figure 8
for i = 1:500
    IK = fcn_IK(20*cos(2*pi*i/500),15*sin(4*pi*i/500)+20);
    Voltages(i,1) = double(IK(1,1));
    Voltages(i,2) = double(IK(1,2));
end
%% figure 8 V2
for i = 1:500
    IK = fcn_IK(16*cos(2*pi*i/500),15*sin(4*pi*i/500)+25);
    Voltages(i,1) = double(IK(1,1));
    Voltages(i,2) = double(IK(1,2));
end
%% straight line
for i = 1:250
    IK = fcn_IK2(i/12.5-10,10);
    Voltages(i,1) = double(IK(1,1));
    Voltages(i,2) = double(IK(1,2));
end
for i = 251:500
    IK = fcn_IK(-(i-250)/12.5+10,10);
    Voltages(i,1) = double(IK(1,1));
    Voltages(i,2) = double(IK(1,2));
end
%% ellipise
for i = 1:500
    IK = fcn_IK(16*cos(2*pi*i/500-pi/2),22*sin(2*pi*i/500-pi/2)+25);
    Voltages(i,1) = double(IK(1,1));
    Voltages(i,2) = double(IK(1,2));
end

%% figure 8 V3
for i = 1:500
    IK = fcn_IK(100*cos(2*pi*i/500),25*sin(4*pi*i/500)+46);
    Voltages(i,1) = double(IK(1,1));
    Voltages(i,2) = double(IK(1,2));
end