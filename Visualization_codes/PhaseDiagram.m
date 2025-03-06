clear;

PhaseDiag = struct('J0',{},'I0',{},'phase',{});


I0 = 1e-3*[(1:4), (-4:-1), 6, -6];

for it1 = (1:numel(I0))
    Temp = struct('J0',{},'I0',{},'phase',{});
    J0 = (-1: 2/10: 1);
    for it2 = (1:numel(J0))
        Temp(it2).I0 = I0(it1);
        Temp(it2).J0 = J0(it2);
        Temp(it2).phase = 'FL';
    end
    PhaseDiag = cat(2,PhaseDiag,Temp);
end
%--------------------------------------

I0 = 1e-3*[(1:3),-(1:3)];

for it1 = (1:numel(I0))
    Temp = struct('J0',{},'I0',{},'phase',{});
    J0 = (-0.2: 0.4/10: 0.2);
    for it2 = (1:numel(J0))
        Temp(it2).I0 = I0(it1);
        Temp(it2).J0 = J0(it2);
        Temp(it2).phase = 'FL';
    end
    PhaseDiag = cat(2,PhaseDiag,Temp);
end
%--------------------------------------

I0 = 1e-4*[1,-1];

for it1 = (1:numel(I0))
    Temp = struct('J0',{},'I0',{},'phase',{});
    J0 = (-1e-3: 2e-3/10: 1e-3);
    for it2 = (1:numel(J0))
        Temp(it2).I0 = I0(it1);
        Temp(it2).J0 = J0(it2);
        Temp(it2).phase = 'FL';
    end
    PhaseDiag = cat(2,PhaseDiag,Temp);
end
%--------------------------------------

I0 = 1e-5*[1,-1];

for it1 = (1:numel(I0))
    Temp = struct('J0',{},'I0',{},'phase',{});
    J0 = (-1e-4: 2e-4/10: 1e-4);
    for it2 = (1:numel(J0))
        Temp(it2).I0 = I0(it1);
        Temp(it2).J0 = J0(it2);
        Temp(it2).phase = 'FL';
    end
    PhaseDiag = cat(2,PhaseDiag,Temp);
end


% I0 width check
J0 = [0.005,0.01,0.015,0.1,0.2,0.3,0.4,0.5,0.6];
J0 = [J0,-J0];
J0 = [J0,0];
J0 = [J0,0.02,0.03,0.04,0.045,0.05,0.06];

for it1 = (1:numel(J0))
    Temp = struct('J0',{},'I0',{},'phase',{});
    I0 = [2.^(-24:0), -2.^(-24:0)];
    for it2 = (1:numel(I0))
        Temp(it2).J0 = J0(it1);
        Temp(it2).I0 = I0(it2);
        Temp(it2).phase = 'FL';
    end
    PhaseDiag = cat(2,PhaseDiag,Temp);
end

% NFL phases
J0 = [0.005*(1:3),0.005*(-3:-1),-1,-0.03125,0,0.03125];

Temp = struct('J0',{},'I0',{},'phase',{});
for it = (1:numel(J0))
    Temp(it).I0 = 0;
    Temp(it).J0 = J0(it);
    Temp(it).phase = 'NFL';
end
PhaseDiag = cat(2,PhaseDiag,Temp);


J0 = [1];

Temp = struct('J0',{},'I0',{},'phase',{});
for it = (1:numel(J0))
    Temp(it).I0 = 0;
    Temp(it).J0 = J0(it);
    Temp(it).phase = 'NFL';
end
PhaseDiag = cat(2,PhaseDiag,Temp);


figure;
hold on;
xlim([-1.1,1.1]);
ylim([-1.1,1.1]);
%ylim([-0.06,0.06])
ax = gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
title('Zero Temperature Phase Diagram','FontSize',30);
set(gca,'XScale','linear','YScale','linear','fontsize',30);
%symlog(gca,'x',-9);
xlabel('$I_{0}$','Interpreter','latex','FontSize',30);
ylabel('$J_{0}$','Interpreter','latex','FontSize',30);

for it = (1:numel(PhaseDiag))

    switch PhaseDiag(it).phase
        case 'FL'
            p(1) = scatter(PhaseDiag(it).I0, PhaseDiag(it).J0, [], 'red', 'x', 'LineWidth', 0.1);
        case 'NFL'
            p(2) = scatter(PhaseDiag(it).I0, PhaseDiag(it).J0, [], 'blue', 'x', 'LineWidth', 5);
        otherwise
            fprinf(['WRN: unknown phase ''',PhaseDiag(it).phase],'''');
    end
end

legend(p, {'FL','NFL'},'Location','northeast','FontSize',25);
hold off;