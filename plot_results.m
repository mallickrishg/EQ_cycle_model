% Script to plot results from running an earthquake cycle.
% Make 4 plots: 
% 1. Time/Time-step evolution along with time series of relevant VS/VW regions
% 2. Cumulative Slip coloured by EQ/post-seismic/inter-seismic
% 3. Frictional parameters - model setup
% 4. Surface displacement time-series 
% Rishav Mallick, EOS, 2018

% load data from results folder
clear
% load results/result_data_test.mat
load results/result_shallow_transients.mat
%% set parameters for plotting
% downsampling in time
down = 1;
plotmovie=0;
% threshold for earthquakes
Veq = 1e-3;
% make new colormap
cspec=[cmap('steelblue',100,10,40);...
%     cmap('skyblue',50,40,37);...
    (cmap('seagreen',50,60,10));...
%     cmap('lawngreen',20,40,49);...
%     cmap('greenyellow',50,63,25);...
    flipud(cmap('orange',100,57,10));...
    flipud(cmap('orangered',50,40,25))];

%% extract relevent timesteps
del_t = zeros(size(t));
del_t(2:end) = t(2:end) - t(1:end-1);

j = 1;
% initialize arrays for Yn = integrated array; 
% Vn = velocity extracted from array of derivatives
Yn = [];
Vn=[];
index =[];

for i = 1:length(t)-1
    if i == 1
        index(j) = i;
        j = j+1;
    elseif i > 1 && max(V(i,:)) >=Veq  
        if ((t(i) - t(index(j-1))) > 3)
            index(j) = i;
            j = j+1;
        end
    else
        if ((t(i) - t(index(j-1))) > 40/365*3.15e7)
            index(j) = i;
            j = j+1;
        end
    end
end
Yn = Y(index',:);
Vn = V(index',:);
%% calculate sliprate in VS regions treating it as a stress-free crack
%%% driven by back-slip and in the shadow of VW regions
sig0 = ss.K*ss.V_plate;
Ipin = false(ss.M,1);
Ipin(top+0:bottom-0) = true;
subK = ss.K(~Ipin,~Ipin);
subsig0 = sig0(~Ipin);

ss.sfc_v = zeros(ss.M,1);
dummysliprate = subK\subsig0;
ss.sfc_v(~Ipin) = dummysliprate;
vratio = ss.sfc_v./ss.V_plate;
%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                    F I G U R E S                     %
%                                                      %
%                   Function of Time                   %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
down = 10;
% plot Time and Time-step evolution
figure(1);clf;set(gcf,'name','Time Evolution')
subplot(3,2,[1 3])
% pcolor(t(1:down:end-1)./3.15e7,ss.y3f/1e3,log10(V(1:down:end,:)')), shading flat
% pcolor(t(1:down:end-1)./3.15e7,ss.y3f/1e3,(1/ss.V_plate(1).*V(1:down:end,:)')), shading flat
pcolor(t(1:down:end-1)./3.15e7,ss.y3f/1e3,log10(V(1:down:end,:)')), shading flat

set(gca,'YDir','reverse');h=colorbar('Location','NorthOutside');
colormap(cspec);
title(h,'log_{10} V'),ylabel('Depth (km)');
caxis([-12 -6])
% caxis([-10 -8])
set(gca,'XTickLabel','','FontSize',15)

subplot(3,2,5);cla;
plot(t(1:down:end-1)/3.15e7,(mean(V(1:down:end,topup:top),2)./ss.V_plate(1)),'b-','LineWidth',1),hold on
plot(t(1:down:end-1)/3.15e7,(mean(V(1:down:end,n:topup),2)./ss.V_plate(1)),'r-','LineWidth',2), hold on
plot(t(1:down:end-1)/3.15e7,(mean(V(1:down:end,top:bottom),2)./ss.V_plate(1)),'m-','LineWidth',1), hold on
% plot(t(1:down:end-1)/3.15e7,(min(V(1:down:end,:),[],2)./ss.V_plate(1)),'b-','LineWidth',1)

% plot(t(1:down:end-1)/3.15e7,repmat(mean(vratio(1:top)),1,length(t(1:down:end-1))),'b--','LineWidth',2)
% plot(t(1:down:end-1)/3.15e7,repmat(mean(vratio(bottom:end)),1,length(t(1:down:end-1))),'m--','LineWidth',2)

xlabel('Time (Yr)','Fontsize',15),ylabel('V/V_{pl}')
axis tight,grid on,ylim([1e-3 1e2])
%title('Mean Velocity for VS regions','Fontsize',15','Fontweight','normal')
set(gca,'FontSize',15,'Color','none','YScale','log')
% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                Function of Time Steps                %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
subplot(3,2,[2 4]);cla;
imagesc(1:down:length(t)-1,ss.y3f/1e3,log10(V(1:down:end,:)')), shading flat
% pcolor(index',ss.y3f/1e3,log10(Vn')), shading flat
% pcolor(3.75e5:length(t)-1,ss.y3f(1:end)/1e3,log10(V(3.75e5:end,1:end)')), shading interp
set(gca,'YDir','reverse');
h=colorbar('Location','NorthOutside');h.Label.String='log_{10}V (m/s)';
ylabel('Depth (km)');
caxis([-12 -1])
set(gca,'XTickLabel','','FontSize',15)

subplot(3,2,6);cla;
plot(1:down:length(t)-1,(max(V(1:down:end,1:top),[],2)),'b-','LineWidth',1),hold on
plot(1:down:length(t)-1,(max(V(1:down:end,top:end),[],2)),'m-','LineWidth',1)
% plot(1:down:length(t)-1,log10(max(V(1:down:end,top:bottom),[],2)),'r-','LineWidth',1)
xlabel('Time Steps','Fontsize',15),ylabel('V_{max}')
axis tight, grid on
% title('Maximum slip velocity for VS regions','Fontsize',15,'Fontweight','normal')
set(gca,'FontSize',15,'Color','none','YScale','log')

%% friction properties and model setup
figure(2);clf;
plot(ss.a-ss.b,ss.y3f/1000,'k-','LineWidth',3)
hold on
plot(ss.a,ss.y3f/1000,'b-','LineWidth',2)
plot(ss.b,ss.y3f/1000,'r-','LineWidth',2)
legend('a-b','a','b')
set(legend,'box','off')
box on
axis tight
grid on
set(gca,'YDir','reverse','FontSize',15,'Color','none')
ylabel('Depth (km)')
xlabel('a-b')
%% cumulative slip
figure(3),clf
skip = 10;
count = 1;
for i = 1:(length(Yn(:,1))/1)
% for i = round(length(Yn(:,1))/2):length(Yn(:,1))
    if max(Vn(i,:)) >= Veq 
        plot(ss.y3f/1e3,Yn(i,1:ss.dgfF:end),'-','LineWidth',0.01,'Color',rgb('orangered')), hold on
    elseif (max(Vn(i,:))<Veq && max(Vn(i,:))>1e1*ss.V_plate(1))
        plot(ss.y3f/1e3,Yn(i,1:ss.dgfF:end),'-','LineWidth',1,'Color',rgb('forestgreen')), hold on
    else
        if count>=skip
            plot(ss.y3f/1e3,Yn(i,1:ss.dgfF:end),'-','LineWidth',1,'Color',rgb('steelblue')), hold on
            count = 1;
        else
            count = count+1;
        end
    end    
end
xlabel('Distance along fault (km)')
ylabel('Cumulative Slip (m)')
axis tight, grid off
set(gca,'FontSize',15,'Color','none')

%% compute displacement field
ndisp = 50;
xplt = linspace(1e3,50e3,ndisp)';
ss.G = zeros(length(xplt),ss.M);
discol = jet(ndisp);

for i = 1:ss.M
    ss.G(:,i) = d1(xplt,ss.y3f(i)+Wf(i)/2,Wf(i));
end
% calculate displacement timeseries
disp = ss.G*Yn(:,1:ss.dgfF:end)';
for i =1:length(index)
    disp(:,i) = disp(:,i)+ss.V_plate(1)/pi*atan2(xplt,Transition).*t(index(i));
end

figure(4),clf
for i = 1:2:length(disp(:,1))
    plot(t(index)./3.15e7,disp(i,:),'-','LineWidth',1,'Color',discol(i,:)),hold on
end
axis tight, grid on
% xlim([110 193])
xlabel('Time (yr)'),ylabel('Displacement (m)')
set(gca,'FontSize',15,'Color','none')
% writetable(table([1:length(xplt)]',xplt),'results/station_locations.dat','Delimiter','\t','WriteVariableNames',0)
% writetable(table(t(index)./3.15e7,disp'),'results/displacement_timeseries.dat','Delimiter','\t','WriteVariableNames',0)

% %% experimental bit
% % plot stresses and show pertubation
% figure(6),clf
% 
% tpertplot = t>tperturb;
% tmod = t(tpertplot) - tperturb;
% tpertplot(end) = [];
% subplot(2,1,1)
% % plot(tmod(1:end-1)/3.15e7,(max(V(tpertplot,bottom:end-n),[],2)./ss.V_plate(1)),'b-','LineWidth',1),hold on
% plot(tmod(1:end-1),(max(V(tpertplot,1:bottom),[],2)./ss.V_plate(1)),'m-','LineWidth',1)
% % plot(tmod(1:end-1)/3.15e7,(mean(V(tpertplot,end-n:end),2)./ss.V_plate(1)),'r-','LineWidth',1)
% 
% xlabel('Time (s)','Fontsize',15),ylabel('V/V_{pl}')
% axis tight,grid on,
% % ylim([0 4])
% title('Velocity for VS regions','Fontsize',15','Fontweight','normal')
% set(gca,'FontSize',15,'Color','none','YScale','log','XScale','log')
% 
% subplot(2,1,2)
% dummy = (V(tpertplot,:)./ss.V_plate(1));
% dummy = dummy(1:35:end,:);
% colmap = jet(length(dummy(:,1)));
% for i = 1:length(dummy(:,1))
%     plot(ss.y3f./1e3,dummy(i,:),'LineWidth',.5,'Color',colmap(i,:)), hold on
%     axis tight, grid on
%     xlabel('X (on fault distance) (km)')
%     ylabel('V/V_{pl}')
%     set(gca,'FontSize',15,'Color','none','YScale','log')
% %     pause(0.02)
% end
% 
% axis tight, grid on