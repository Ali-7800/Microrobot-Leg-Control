%% Find workspace where beam makes no contact
space_no_contact = zeros(300);
middle = [];
lim = 0;
k = 1;
while(lim == 0)
    IK = fcn_IK(0,k);
    if isnan(IK(1,1))||isnan(IK(1,2))
        lim =1;
    else
        middle = [middle;1];
        k = k + 1
    end
end
%%
I = 1;
while(I<k)
    lim = 0;
    J = 1;
    while(J<30 && lim==0)
        a = 1;
        IK = fcn_IK(J*10,I);
        contact = fcn_contactCheck(J*10,I);
        if isnan(IK(1,1))||isnan(IK(1,2))||contact>0
           while(lim==0 && a<11)
                IK = fcn_IK((J-1)*10+a,I);
                contact = fcn_contactCheck((J-1)*10+a,I);
                if isnan(IK(1,1))||isnan(IK(1,2))||contact>0
                    lim = 1;
                else
                    space_no_contact(I,(J-1)*10+a)=1; 
                    a=a+1
                end
           end
        else
            space_no_contact(I,J*10)=1;
            J=J+1
        end
    end
    I = I + 1
end

%% find workspace (IK method slow)
%%steps to follow
%1-run this one time and save space as space_with_contacts
%2-disable contact model in IK and run again and save space as space_without_contacts
space = zeros(300);
I = 1;
while(I<k)
    lim = 0;
    J = 1;
    while(J<30 && lim==0)
        a = 1;
        IK = fcn_IK(J*10,I);
        if isnan(IK(1,1))||isnan(IK(1,2))
           while(lim==0 && a<11)
                IK = fcn_IK((J-1)*10+a,I);
                if isnan(IK(1,1))||isnan(IK(1,2))
                    lim = 1;
                else
                    space(I,(J-1)*10+a)=1; 
                    a=a+1
                end
           end
        else
            space(I,J*10)=1;
            J=J+1
        end
    end
    I = I + 1
end
%% find workspace (FK method fast and more accurate)
%%steps to follow
%1-run this one time and save deltas as deltas_with_contacts
%2-disable contact model in FK and run again and save deltas as deltas_without_contacts
p = get_params();
params = p.params;
dV_max = params(7);
dV1 = 0;
deltas = [[0,0,0,0]];
for dV2 = 1:dV_max
    FK = fcn_FK(dV1,dV2);
    deltas = [deltas;[dV1,dV2,FK(1),FK(2)]];
end
for dV1 = 1:dV_max
    FK = fcn_FK(dV1,dV2);
    deltas = [deltas;[dV1,dV2,FK(1),FK(2)]];
end




