%file: skarn_post_process.m

%author: Evan J. Ramos
%date:   15 Jan 2018

%Description: Example code for post processing data from skarn_main.m or
%             skarn_generic.m. NOTE: you must run scripts before the post
%             processing can take place


close all
clc

%% Tracer distribution values

if exist('Pluton','var')
    all_ages = Ages_lime + Ages_plut + Ages_surf + Ages_volc;
    AGE_L = Ages_lime./all_ages; AGE_L(isnan(AGE_L)) = .25;
    AGE_P = Ages_plut./all_ages; AGE_P(isnan(AGE_P)) = .25;
    AGE_S = Ages_surf./all_ages; AGE_S(isnan(AGE_S)) = .25;
    AGE_V = Ages_volc./all_ages; AGE_V(isnan(AGE_V)) = .25;
    AGE_L_ROCK = AGE_L(Pluton.dof_rock,:);
    AGE_P_ROCK = AGE_P(Pluton.dof_rock,:); 
    AGE_S_ROCK = AGE_S(Pluton.dof_rock,:); 
    AGE_V_ROCK = AGE_V(Pluton.dof_rock,:);
end

%% Thermal evolution

T_sol        = 710; %[C] solidus temp
T_gtop       = 575; %[C] max garnet stability temp
T_gbot       = 375; %[C] min garnet stability temp

temp_ts = t;
T_plots = Ts-273;

figure('Units','Normalized','Position',[0 .6 .2 .8])%,'Renderer','painters');
for i = 1:length(temp_ts)
    colormap gray(20)
    contourf(Xc,Yc,reshape(T_plots(:,i),Grid.Ny,Grid.Nx),20,...
             'LineStyle','None')
    axis equal
    hold on
    plot_pluton(Pluton.dof,[],Grid)
    contour(Xc,Yc,reshape(T_plots(:,i),Grid.Ny,Grid.Nx),[T_sol T_sol],'LineColor','r','LineWidth',2)
    contour(Xc,Yc,reshape(T_plots(:,i),Grid.Ny,Grid.Nx),[T_gtop T_gtop],'LineColor',[1 .5 0],'LineWidth',2)
    contour(Xc,Yc,reshape(T_plots(:,i),Grid.Ny,Grid.Nx),[T_gbot T_gbot],'LineColor',[1 1 0],'LineWidth',2)
    contour(Xf,Yf,reshape(Psi(:,i),Grid.Ny+1,Grid.Nx+1),10,'LineColor','c','LineWidth',1)
    set(gca,'XTick',[],'YTick',[],'FontSize',18)
    colorbar
    caxis([min(T_plots(:)) max(T_plots(:))])
    drawnow
    hold off
end

%% Points A-F total

figure('Units','Normalized','Position',[.2 .2 .8 .6])
NTick = 4;

for i = 1:length(Pluton.dof_rock)
    
    %d18O
    subplot(length(Pluton.dof_rock),2,i*2 - 1)
    [hAx,hLine1,hLine2] = plotyy(t,delta_ws(Pluton.dof_rock(i),:),...
    t,Ts(Pluton.dof_rock(i),:)-273);
    hold on
    plot(t,delta_gt(Pluton.dof_rock(i),:),'Color',[0 0.4470 0.7410],...
        'LineWidth',1)
    xlabel('Time [kyr]')
    ylabel(hAx(1),'\delta^{18}O [VSMOW]')
    ylabel(hAx(2),'T [^{\circ}C]')
    title(char(64 + i))
    hLine1.LineStyle = '--';
    hLine2.LineStyle = '-';
    hLine1.LineWidth = 1;
    hLine2.LineWidth = 1;
    legend('H_2O','Garnet','Location','Northeast')
    minmin = min([delta_gt(Pluton.dof_rock(i),:) delta_ws(Pluton.dof_rock(i),:)]);
    maxmax = max([delta_gt(Pluton.dof_rock(i),:) delta_ws(Pluton.dof_rock(i),:)]);
    hAx(1).XLim = [min(t) max(t)];
    hAx(1).YLim = [minmin maxmax];
    hAx(1).YTick = floor(linspace(minmin,maxmax,NTick));
    hAx(2).XLim = [min(t) max(t)];
    minT        = min(Ts(Pluton.dof_rock(i),:)) - 273;
    maxT        = max(Ts(Pluton.dof_rock(i),:)) - 273;
    hAx(2).YTick = floor(linspace(minT,maxT,NTick));
    hAx(2).YLim = [min(Ts(Pluton.dof_rock(i),:)) max(Ts(Pluton.dof_rock(i),:))]-273;
    set(gca,'XGrid','on')
    
    %tracer
    subplot(length(Pluton.dof_rock),2,i*2)
    plot(t,Ages_surf(Pluton.dof_rock(i),:),'r',...
         t,Ages_volc(Pluton.dof_rock(i),:),'b',...
         t,Ages_plut(Pluton.dof_rock(i),:),'g',...
         t,Ages_lime(Pluton.dof_rock(i),:),'m',...
         'LineWidth',1)
    grid on
    xlabel('Time [kyr]')
    ylabel('C')
    title(char(64 + i))
    legend('Surface','Metavolcanic','Magmatic','Limestone','Location','Northeast')
    axis([0 t(end) 0 1])
    
end

%% Solely when garnet is stable

%Isolating conditions
T_top = 575 + 273;
T_bot = 375 + 273;
range = Ts > T_top | Ts < T_bot;
delta_ws(range) = NaN;
delta_gt(range) = NaN;

figure('Units','Normalized','Position',[0 0 1 1])
%delta_gt = delta_gt + .5;
for i = 1:length(Pluton.dof_rock)
    
    %d18O
    subplot(length(Pluton.dof_rock),2,i*2 - 1)
    h = plot(t,delta_ws(Pluton.dof_rock(i),:),'k-',...
             t,delta_gt(Pluton.dof_rock(i),:),'k--');
    h(1).LineWidth = 1;
    h(2).LineWidth = 1;
    title(char(64 + i))
    minx = min([delta_ws(Pluton.dof_rock(i),:) delta_gt(Pluton.dof_rock(i),:)]);
    maxx = max([delta_ws(Pluton.dof_rock(i),:) delta_gt(Pluton.dof_rock(i),:)]);
    xlim([0 10])
    set(gca,'XTick',[])
    if i == 6
        ylim([-20 0])
    elseif i == 4
        ylim([-15 5])
    elseif i == 5
        ylim([-15 0])
    end
    set(gca,'Xcolor',[0 0 0],'YColor',[0 0 0],'FontSize',15)
    if i == 1
        legend('H_2O','Garnet','Location','Northeast')
    elseif i == length(Pluton.dof_rock)
        xlabel('Time [kyr]')
    end
    
    %tracer
    subplot(length(Pluton.dof_rock),2,i*2)
    area(t,[AGE_L_ROCK(i,:)' AGE_P_ROCK(i,:)' AGE_S_ROCK(i,:)' AGE_V_ROCK(i,:)'],...
         'LineStyle','None')
    xlim([0 10])
    ylim([0 1])
    set(gca,'XTick',[],'YTick',[],'Xcolor',[0 0 0],'YColor',[0 0 0],'FontSize',15)
    if i == 1
        legend('Limestone','Magmatic','Surface','Metavolcanic')
    end
  
end


%% Generate video
%Stable isotope
figure('Units','Normalized','Position',[.2 .2 .4 .6])
for i = 1:length(t)
    colormap parula(20)
    contourf(Xc,Yc,reshape(delta_ws(:,i),Grid.Ny,Grid.Nx),...
             20,'Linestyle','None')
    hold on
    set(gca,'XTick',[],'YTick',[])
    plot_pluton(Pluton.dof,[],Grid)
    c = colorbar;
    caxis([min(delta_ws(:)) max(delta_ws(:))])
    axis equal
    title(sprintf('t = %.2f kyr',t(i)),'FontSize',20)
    if i > 1
        contour(Xf,Yf,reshape(Psi(:,i-1),Grid.Ny+1,Grid.Nx+1),20,'Color','k','LineWidth',1)
    end
    pause
    hold off

end