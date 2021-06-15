clear
addpath ~/Dropbox/scripts/unicycle/matlab/
import unicycle.*

%% 1. compute greens functions and stress kernels
% use unicycle to create the 'rcv' object and stress kernels
% setup fault solution type and G, nu
G = 30e3; nu = 0.25;
earthModel=greens.okada92(G,nu);

patchfname='flt_2D_testing.seg';
% rcv is the fault
rcv = geometry.receiver(patchfname,earthModel);
% self-stress kernels on the fault
% load stress kernels if they aren't already in the workspace
if ~exist('K','var')
    if exist('stresskernels/K.mat','file')
        disp('Loading K from available file')
        load stresskernels/K.mat
        if length(K)~=rcv.N
            disp('reloading appropriate kernels and fault file')                        
            [~,~,~,K,~,~]=earthModel.tractionKernels(rcv,rcv);
            save('stresskernels/K.mat','K')
        end
    else
        disp('No K in folder, will create traction kernels')
        [~,~,~,K,~,~]=earthModel.tractionKernels(rcv,rcv);
        save('stresskernels/K.mat','K')
    end
end
% define dip-distance in (m)
dipdist = -diag(rcv.xc*rcv.dv');


% minimum grid size
% fprintf(1,'grid size = %.2f m, minimum grid size = %.2f m\n',rcv.W(1),min([G*rcv.l(top)/rcv.b(top)/rcv.sigma(top) G*rcv.l(topvw-1)/rcv.b(topvw-1)/rcv.sigma(topvw-1)])) 
% fprintf(1,'a/b in VS regions = %.2f\n',rcv.a(1)/rcv.b(1))

% add time of perturbation
tperturb = 10*3.15e7;
% stress perturbation
deltau = [-4];
% total duration of simulation
tend = 1500*3.15e7;

%% expected slip rate

tau = K*rcv.Vpl;
subK = K(~rcv.pinnedPosition,~rcv.pinnedPosition);
subtau = tau(~rcv.pinnedPosition);
Vpred = zeros(rcv.N,1);
Vpred(~rcv.pinnedPosition) = (subK\subtau);
Vpred = Vpred./rcv.Vpl;

%% Vary parameter of interest and SOLVE (will be saved to directory)
thetalaw = 1; % 1 - ageing law, 2 - slip law

Vplfrac = 0.3;

% vary parameter of choice
nexp = 15;
lvec = logspace(-2,-1,nexp);
topvec = linspace(10e3,25e3,nexp);

% lvec = [0.015,0.02];
% topvec = [10e3,12e3];
[L,T] = meshgrid(lvec,topvec);

tduration = zeros(size(L));
tssestart = tduration; tsseend = tssestart;
Vpeak = tduration;

parfor i = 1:length(L(:))
    disp(['Iteration number ' num2str(i)])
    rcv = unicycle.geometry.receiver(patchfname,earthModel);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %                                                      %
    %         F R I C T I O N   P A R A M E T E R S        %
    %                                                      %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
    fault_width = max(dipdist);
    % create velocity weakening region
    VW = zeros(rcv.N,1);
    topvw = floor(50e3/(fault_width/rcv.N));
    botvw = ceil(60e3/(fault_width/rcv.N));
    VW(topvw:botvw) = 1;
    VW = logical(VW);
    % make VW region pinned
    rcv.pinnedPosition = VW;
    % effective confining pressure on fault (MPa) using depth dependent normal
    rcv.sigma = 100.*ones(rcv.N,1);
    
    
    % frictional parameters
    rcv.a = 1e-2.*ones(rcv.N,1);
    % if you want to set the entire thing as VS
    rcv.b = rcv.a - 9e-3;
    
    % static friction coefficient
    rcv.mu0 = 0.1*ones(rcv.N,1);
    % numer of state parameters
    rcv.dgf=4;
    % characteristic weakening distance (m)
    rcv.l = 0.05.*ones(rcv.N,1);
    % plate velocity (m/s)
    rcv.Vpl = 1e-9*ones(rcv.N,1);
    % reference slip rate (m/s)
    rcv.Vo = 1e-6*ones(rcv.N,1);
    % shear wave speed (m/s)
    rcv.Vs = 3e3.*ones(rcv.N,1);
    rcv.b = rcv.a - 9e-3;
    top    = floor(T(i)/(fault_width/rcv.N));
    bottom = ceil(30e3/(fault_width/rcv.N));
    
    rcv.b(top:bottom) = rcv.a(top:bottom) - 0.1e-3;
    rcv.l = L(i).*ones(rcv.N,1);
    
    [Y,V,t] = func_stressperturb_pinnedeqcycles(rcv,tperturb,deltau,top,bottom,tend,thetalaw);
    
    % find SSE
    Vmean_sse = mean(V(find(t==tperturb(1),1):end,top:bottom),2);
    Vpeak(i) = max(Vmean_sse);
    
    indexstart_sse = find(Vmean_sse >= Vplfrac*rcv.Vpl(1),1) + find(t==tperturb(1),1);
    if isempty(indexstart_sse)==1
        indexend_sse = NaN;
    else
        tssestart(i) = t(indexstart_sse)./3.15e7; 
        indexend_sse = find((Vmean_sse < Vplfrac*rcv.Vpl(1)) & (t(find(t==tperturb(1),1):end-1)>t(indexstart_sse)),1) + find(t==tperturb(1),1);
        
        if isempty(indexend_sse)==1
            tsseend(i) = NaN;
        else
            tsseend(i) = t(indexend_sse)./3.15e7;
        end
        
        % SSE duration
        tduration(i) = tsseend(i) - tssestart(i);
    end
    
    
end
%% plot results

figure(1),clf
% pcolor(lvec,topvec./1e3,tduration)
imagesc(log10(lvec),topvec./1e3,tduration)
xlabel('log_{10} D_c'),ylabel('top location (km)')
cb=colorbar;
cb.Label.String = 'SSE Duration (yrs)';
colormap(jet(20))
caxis([0 50])
set(gca,'YDir','normal','FontSize',15)
print('Figures/SSE_duration_2','-djpeg','-r300')

figure(2),clf
imagesc(log10(lvec),topvec./1e3,tssestart)
xlabel('log_{10} D_c'),ylabel('top location (km)')
cb=colorbar;
cb.Label.String = 'SSE Start Time (yrs)';
colormap jet
% caxis([0 2])
set(gca,'YDir','normal','FontSize',15)

figure(3),clf
imagesc(log10(lvec),topvec./1e3,tsseend)
xlabel('log_{10} D_c'),ylabel('top location (km)')
cb=colorbar;
cb.Label.String = 'SSE End Time (yrs)';
colormap jet
% caxis([0 2])
set(gca,'YDir','normal','FontSize',15)


figure(4),clf
imagesc(log10(lvec),topvec./1e3,Vpeak./1e-9)
xlabel('log_{10} D_c'),ylabel('top location (km)')
cb=colorbar;
cb.Label.String = 'V_{peak}/V_{pl}';
colormap(jet(20))
caxis([0 2])
set(gca,'YDir','normal','FontSize',15)
print('Figures/SSE_Vpeak_2','-djpeg','-r300')



