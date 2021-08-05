%% plot workspace if workspace is calculated using IK method
no_contact = zeros([k 2]);
without_contact = zeros([k 2]);
with_contact = zeros([k 2]);

for i=1:k-1
    no_contact(i,2) = i;
    with_contact(i,2) = i;
    without_contact(i,2) = i;
    for j=1:300
        if space_no_contacts(i,j)==1
           no_contact(i,1) = j;
        end
        if space_with_contacts(i,j)==1
            with_contact(i,1) = j;
        end
        if space_without_contacts(i,j)==1
            without_contact(i,1) = j;
        end
    end
end
hold on

fill([-flip(with_contact(:,1));with_contact(:,1)],[flip(with_contact(:,2));with_contact(:,2)],[0 0.4470 0.7410])
fill([-flip(without_contact(:,1));without_contact(:,1)],[flip(without_contact(:,2));without_contact(:,2)],[0.8500 0.3250 0.0980])
fill([-flip(no_contact(:,1));no_contact(:,1)],[flip(no_contact(:,2));no_contact(:,2)],[0.9290 0.6940 0.1250])
xlim([-300 300])
% axis equal
ylim([0 100])
xlabel("\delta{x_c}(\mum)")
ylabel("\delta{y_c}(\mum)")
% legend("Contact point amplified workspace","Unamplified workspace","Workspace with no contact")
title("Workspace for \DeltaV_{max}=15V")


%% plot workspace if workspace is calculated using FK method
no_contact = zeros([k 2]);

for i=2:k-1
    no_contact(i,2) = i;
    for j=1:300
        if space_no_contact(i,j)==1
           no_contact(i,1) = j;
        end
    end
end
hold on
fill([flip(-deltas_with_contact(:,3));deltas_with_contact(:,3)],[flip(deltas_with_contact(:,4));deltas_with_contact(:,4)],[0 0.4470 0.7410])
fill([flip(-deltas_without_contact(:,3));deltas_without_contact(:,3)],[flip(deltas_without_contact(:,4));deltas_without_contact(:,4)],[0.8500 0.3250 0.0980])
fill([-flip(no_contact(:,1));no_contact(:,1)],[flip(no_contact(:,2));no_contact(:,2)],[0.9290 0.6940 0.1250])
xlabel("\delta{x_c}(\mum)")
ylabel("\delta{y_c}(\mum)")

%% draw trajectory
for i=1:500
    x = 100*cos(2*pi*i/500);
    y = 25*sin(4*pi*i/500)+46;
    plot(x,y,"*g")
    drawnow;
end