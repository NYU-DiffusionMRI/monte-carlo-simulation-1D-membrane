% In this example, we analyze Monte Carlo simualtion results in five 1d 
% micro-geometries consisting of permeable membranes (demo3_simulation.m)
%
% Author: Hong-Hsi Lee, December, 2019 (orcid.org/0000-0002-3663-6559)

% set up the root and target directories on your own computer
% root = '/directory/to/this/file';
root = '.';
root_input = fullfile(root,'hpc_code/input');
root_lib = fullfile(root,'hpc_code/lib');

%% Calculate diffusivity and kurtosis based on moments of diffusion displacement

files = dir(fullfile(root_input,'membrane_*'));
dx2 = 0;                            % Second order moment <x^2>, micron^2
dx4 = 0;                            % Fourth order moment <x^4>, micron^4
NPar = 0;                           % # random walker
sig = 0;                            % Signal
for i = 1:numel(files)
    rooti = fullfile(root_input,files(i).name);
    dx2i = load(fullfile(rooti,'dx2_diffusion.txt'));
    dx4i = load(fullfile(rooti,'dx4_diffusion.txt'));
    dx2 = dx2+dx2i;
    dx4 = dx4+dx4i;
    sigi = load(fullfile(rooti,'sig_diffusion.txt'));
    sig = sig+sigi;
    para = load(fullfile(rooti,'sim_para.txt'));
    NPar = NPar+para(3);
end
t = load(fullfile(rooti,'diff_time.txt'));  % Diffusion time, ms
bval = load(fullfile(rooti,'bval.txt'));    % b-value, ms/micron^2
dt = para(1);                       % Time-step, ms
D0 = para(4);                       % Intrinsic diffusivity, micron^2/ms
kappa0 = para(5);                   % Input value of permeability, micron/ms

dx2 = dx2/NPar;
dx4 = dx4/NPar;
D = dx2/2./t;                       % Diffusivity given by moments
K = dx4./dx2.^2-3;                  % Kurtosis given by moments
sig = sig/NPar;

% The genuine permeability in simulation is always larger than the input
% value. Here, we corrected the permeability based on Eq. (B3) in Appendix
% B of (Lee and Papaioannou, et al., NeuroImage, 2020).
dx = sqrt(2*D0*dt);
kappa = kappa0/(1-kappa0*dx/D0);

% Calculate the bead distance in numerical phantoms
at = [];                            % Membrane distance
for i = 1:numel(files)
    rooti = fullfile(root_input,files(i).name);
    xm = load(fullfile(rooti,'phantom_xMem.txt'));  % Normalized membrane position
    vox = load(fullfile(rooti,'phantom_res.txt'));  % Axonal length, micron
    xm = xm*vox;
    ai = [xm(1)*2; diff(xm)];       % Membrane distance of the i-th axon
    at = cat(1,at,ai);
end
abar = mean(at);                    % Mean membrane distance
astd = std(at);                     % Standard deviation of membrane distance

% Calculate theoretical values of time-dependence parameters in Eqs. (17-18)
% of (Lee and Papaioannou, et al., NeuroImage, 2020).
tr = abar/2/kappa;                  % Mean residence time within a typical membrane distance
zeta = D0/kappa/abar;
Dinf = D0/(1+zeta);                 % Bulk diffusivity in t -> infinity limit
A = Dinf*sqrt(tr/2/pi)*astd^2/abar^2*( zeta/(1+zeta) )^(3/2);
cD = 2*A;                           % Time-dependence amplitude of diffusivity
cK = 4*A/Dinf;                      % Time-dependence amplitude of kurtosis

% Plot figure
figure('unit','inch','position',[0 0 15 5]);
h1 = subplot(131); p1 = h1.Position;
h2 = subplot(132); p2 = h2.Position;
h3 = subplot(133); p3 = h3.Position;
close all

figure('unit','inch','position',[0 0 15 5]);
% D(t) versus t/tr
axes('position',p1);
hold on;
plot(t/tr,D,'b-','linewidth',1);
xlim([0 250]); ylim([0.9 2]);
set(gca,'fontsize',12,'xtick',0:50:1000,'ytick',[0.9 1:0.2:2]);
xlabel('$t/\tau_r$','interpreter','latex','fontsize',20);
ylabel('$D(t)$, $\mu$m$^2$/ms','interpreter','latex','fontsize',20);
pbaspect([1 1 1]); box on; grid on;

% inset: ( D(t)-D_infty )/D_infty versus t/tr, log scale
p1i = [0.18 0.39 0.16 0.43];
h1i = axes('Position',p1i);
hold on;
plot(t/tr,D/Dinf-1,'b-','linewidth',1);
plot(t/tr,cD/Dinf./sqrt(t),'b:','linewidth',1);
set(h1i,'xscale','log','yscale','log');
xlim([1e-1 250]); ylim([1e-3 1e0]);
set(gca,'fontsize',10,'xtick',10.^[-1:2]);
xlabel('$t/\tau_r$','interpreter','latex','fontsize',14);
ylabel('$(D(t)-D_\infty)/D_\infty$','interpreter','latex','fontsize',14);
box on; grid on;

% K(t) versus t/tr
axes('position',p2);
hold on;
plot(t/tr,K,'r-','linewidth',1);
xlim([0 250]); ylim([0 0.5]);
set(gca,'fontsize',12,'xtick',0:50:1000,'ytick',0:0.1:1);
xlabel('$t/\tau_r$','interpreter','latex','fontsize',20);
ylabel('$K(t)$','interpreter','latex','fontsize',20);
pbaspect([1 1 1]); box on; grid on;

% inset: K(t) versus t/tr, log scale
p2i = [0.46 0.39 0.16 0.43];
h2i = axes('Position',p2i);
hold on;
plot(t/tr,K,'r-','linewidth',1);
plot(t/tr,cK./sqrt(t),'r:','linewidth',1);
set(h2i,'xscale','log','yscale','log');
xlim([1e-1 250]); ylim([1e-3 1e0]);
set(gca,'fontsize',10,'xtick',10.^[-1:2]);
xlabel('$t/\tau_r$','interpreter','latex','fontsize',14);
ylabel('$K(t)$','interpreter','latex','fontsize',14);
box on; grid on;

% Diffusion metrics versus sqrt(tr/t)
axes('position',p3);
hold on;
tlist = 2:10000;
hd = plot(1./sqrt(t(tlist)/tr),D(tlist)/Dinf-1,'b-','linewidth',1);
hk = plot(1./sqrt(t(tlist)/tr),K(tlist),'r-','linewidth',1);
xlim([0 2.5]); ylim([0 0.2]);

hdr = refline(cD/Dinf/sqrt(tr),0);
set(hdr,'color','b','linestyle',':','linewidth',1);
hkr = refline(cK/sqrt(tr),0);
set(hkr,'color','r','linestyle',':','linewidth',1);

xlim([0 2.5]); ylim([0 0.5]);
set(gca,'xtick',0:0.5:10,'ytick',0:0.1:10,'fontsize',12);
xlabel('$\sqrt{\tau_r/t}$','interpreter','latex','fontsize',20);
ylabel('Diffusion Metrics','interpreter','latex','fontsize',20);
legend([hd hk],{'$(D(t)-D_\infty)/D_\infty$','$K(t)$'},'interpreter','latex','fontsize',16,'location','northwest');
pbaspect([1 1 1]); grid on; box on;

%% Compare diffusivity and kurtosis based on cumulant and signal

% Fit signals upto bval^2, bval^3, and bval^4 respectively, and estimate
% diffusivity, kurtosis, and time-dependence parameters
Di = zeros(3,1);                    % Bulk diffusivity in t -> infinity limit
cDi = zeros(3,1);                   % Time-dependence amplitude of diffusivity
cKi = zeros(3,1);                   % Time-dependence amplitude of kurtosis

A_bval = [bval bval.^2 bval.^3 bval.^4];
[~,It] = min(abs(t/tr-4));
flist = It:10000;
A_time = [ones(numel(flist),1) 1./sqrt(t(flist)/tr)];

figure('unit','inch','position',[0 0 15 10])
for i = 1:3
    X = A_bval(:,1:(i+1))\log(sig.');
    X = X.';
    Ds = -X(:,1);                   % Diffusivity given by signals
    Ws = X(:,2);
    Ks = 6*Ws./Ds.^2;               % Kurtosis given by signals
    
    Xd = A_time\Ds(flist);
    Xk = A_time(:,2)\Ks(flist);
    Di(i) = Xd(1);
    cDi(i) = Xd(2);
    cKi(i) = Xk;
    
    subplot(2,3,i);
    hold on;
    hc = plot(1./sqrt(t(tlist)/tr),D(tlist),'-','linewidth',1); 
    hs = plot(1./sqrt(t(tlist)/tr),Ds(tlist),'--','linewidth',1);
    set(gca,'fontsize',12,'xtick',0:0.5:10,'ytick',[0.8 1:0.2:2]);
    pbaspect([1 1 1]); box on; grid on;
    xlabel('$\sqrt{\tau_r/t}$','interpreter','latex','fontsize',20);
    ylabel('$D(t)$, $\mu$m$^2$/ms','interpreter','latex','fontsize',20);
    legend([hc hs],{'Cumulant',sprintf('Signal fit up to $b^%u$',i+1)},'interpreter','latex','fontsize',18);
    xlim([0 2.5]); ylim([0.8 1.8]);

    subplot(2,3,i+3);
    hold on;
    hc = plot(1./sqrt(t(tlist)/tr),K(tlist),'-','linewidth',1);
    hs = plot(1./sqrt(t(tlist)/tr),Ks(tlist),'--','linewidth',1);
    set(gca,'fontsize',12,'xtick',0:0.5:10,'ytick',0:0.1:1);
    xlabel('$\sqrt{\tau_r/t}$','interpreter','latex','fontsize',20);
    ylabel('$K(t)$','interpreter','latex','fontsize',20);
    pbaspect([1 1 1]); box on; grid on;
    legend([hc hs],{'Cumulant',sprintf('Signal fit up to $b^%u$',i+1)},'interpreter','latex','fontsize',18,'location','southeast');
    xlim([0 2.5]); ylim([0 0.5]);
end

% Calculate the error of time-dependence parameters for the signal fit upto
% bval^2, bval^3, and bval^4 respectively
Di_err = Di/Dinf-1;                 % Error in D_infty
cDi_err = cDi/(2*A/sqrt(tr))-1;     % Error in c_D
cKi_err = cKi/(4*A/Dinf/sqrt(tr))-1;% Error in c_K

for i = 1:3
    fprintf('The error of time-dependence parameters for the signal fit upto b^%u:\n',i+1);
    fprintf('D_infty: %.2f%%\n',Di_err(i)*100);
    fprintf('c_D: %.2f%%\n',cDi_err(i)*100);
    fprintf('c_K: %.2f%%\n\n',cKi_err(i)*100);
end



