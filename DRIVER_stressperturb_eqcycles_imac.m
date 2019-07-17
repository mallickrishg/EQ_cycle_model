clear
close all
addpath ~/Dropbox/scripts/unicycle/matlab/
import unicycle.*

patchfname='flt_2D_testing.seg';
G = 30e3; nu = 0.25;
earthModel=greens.okada92(G,nu);
% rcv is the fault
rcv = geometry.receiver(patchfname,earthModel);
 
vwparams.a = 0.01;
vwparams.b = vwparams.a + 9e-3;
vwparams.bvs = vwparams.a - 5e-3;
vwparams.top = 60e3;
vwparams.bot = 80e3;

perturbparams.a = 0.01;
perturbparams.b = perturbparams.a - 0.2e-3;
perturbparams.top = 1e3;
perturbparams.bot = 30e3;
perturbparams.tperturb = 200*3.15e7;

nexp = 20;
down = 30;

deltauvec = linspace(-3,-2,nexp);

% save files
odir = '~/Documents/results_stressperturb_eqcycles';
flist = dir([odir '/*.mat']);
fnum = length(flist)+1;

for i = 1:nexp
    i
    perturbparams.deltau = deltauvec(i);
    
    tend = 700*3.15e7;% (yrs * yr2sec)
    
    [Y,V,t] = func_stressperturb_eqcycles(patchfname,vwparams,perturbparams,tend);
    % downsample spatially V and save it
    Vdown = V(:,1:down:end);
    % save files
    save([odir '/stressperturb_eqcycles_' num2str(fnum)],'rcv','vwparams','perturbparams','Vdown','t')
    fnum = fnum+1;
    
    clear Y V t
end

%% transfer all .mat to a new folder
flist = dir([odir '/run*']);
fnum = length(flist);
fol_out = [odir '/run' num2str(fnum+1)];
mkdir(fol_out)
movefile([odir '/*.mat'],fol_out);

%% print figures
odatadir = ['~/Documents/results_stressperturb_eqcycles/run' num2str(fnum+1)];
figodir = ['Figures/imac/perturbations/run' num2str(fnum+1)];
mkdir(figodir)
nexp = length(dir([fol_out '/*.mat']));

for i = 1:nexp
    datafnum = i;
    load([odatadir '/stressperturb_eqcycles_' num2str(datafnum) '.mat'])
    
    figure(1);clf
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

    plot_stressperturb_eqcycles(rcv,vwparams,perturbparams,[],Vdown,down,t)
    print([figodir '/eqcycle_deltau_' num2str(i)],'-djpeg','-r100')
    
    figure(2),clf
    VW = zeros(rcv.N,1);
    topvw = floor(vwparams.top/(rcv.W(1)));
    botvw = ceil(vwparams.bot/(rcv.W(1)));
    VW(topvw:botvw) = 1;
    VW = logical(VW);
    rcv.a = vwparams.a.*ones(rcv.N,1);
    rcv.b = vwparams.bvs.*ones(rcv.N,1); % for VS region
    rcv.b(VW) = rcv.a(VW).*0 + vwparams.b; % for VW region
    % location of perturbation
    top    = floor(perturbparams.top/(rcv.W(1)));
    bottom = ceil(perturbparams.bot/(rcv.W(1)));
    rcv.b(top:bottom) = perturbparams.b;
    %%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %          Friction properties and model setup       %  %
    plot(rcv.a-rcv.b,diag(rcv.xc*rcv.dv')./1e3,'k-','LineWidth',3)
    hold on
    plot(rcv.a,diag(rcv.xc*rcv.dv')/1000,'b-','LineWidth',2)
    plot(rcv.b,diag(rcv.xc*rcv.dv')/1000,'r-','LineWidth',2)
    legend('a-b','a','b')
    set(legend,'box','off')
    box on
    axis tight
    grid on
    set(gca,'FontSize',12,'Color','none')
    ylabel('Depth (km)')
    xlabel('a-b')
    print([figodir '/modelsetup_deltau_' num2str(i)],'-djpeg','-r100')
end