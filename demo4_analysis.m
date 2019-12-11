clear
close all

%% analyze packing
rootpck = '/Users/hhl309/Documents/GitHub/monte-carlo-1D-membrane/hpc_code/input_cuda/membrane_v4';
files = dir(fullfile(rootpck,'membrane*'));
abar = zeros(numel(files),1);
avar = zeros(numel(files),1);
Lt = zeros(numel(files),1);
for i = 1:numel(files)
    filename = files(i).name;
    mx = load(fullfile(rootpck,filename,'phantom_xMem.txt'));
    L = load(fullfile(rootpck,filename,'phantom_res.txt'));
    mx = mx*L;
    ai = [mx(1)-mx(end)+L; diff(mx)];
    abar(i) = mean(ai);
    avar(i) = var(ai);
    Lt(i) = L;
end

%%
root = '/Volumes/labspace/Honghsi/projects/simulation1D_cuda/proj000004';

dx2 = 0; dx4 = 0; NPar = 0; sig = 0;
Di = []; Ki = [];
for i = 1:5
    rooti = fullfile(root,sprintf('setup%06u',i));
    dx2i = load(fullfile(rooti,'dx2_diffusion.txt'));
    dx4i = load(fullfile(rooti,'dx4_diffusion.txt'));
%     dx2 = dx2+fi(i)*dx2i;
%     dx4 = dx4+fi(i)*dx4i;
    sigi = load(fullfile(rooti,'sig_diffusion.txt'));
%     sig = sig+fi(i)*sigi;
    para = load(fullfile(rooti,'sim_para.txt'));
%     NPar = NPar+fi(i)*para(3);
    ti = load(fullfile(rooti,'diff_time.txt'));
    Di = cat(2,Di,dx2i/para(3)/2./ti);
    Ki = cat(2,Ki,dx4i/para(3)./(dx2i/para(3)).^2-3);
end
t = load(fullfile(rooti,'diff_time.txt'));
bval = load(fullfile(rooti,'bval.txt'));
dt = para(1);
D0 = para(4);
kappa = para(5);
% kappa = kappa/(1-kappa*sqrt(dt/2)*(1/sqrt(D0)+1/sqrt(D0)));

%%

tri = abar/2/kappa;
zeta = D0/kappa./abar;
Dinfi = D0./(1+zeta);
Ai = Dinfi.*sqrt(tri/2/pi).*avar./abar.^2.*( zeta./(1+zeta) ).^(3/2);

Xi = zeros(5,2);
for i = 1:5
    xx = t/tri(i);
    [~,I1] = min(abs(xx-10));
    flist = I1:10000;
    Xi(i,:) = [ones(numel(flist),1) 1./sqrt(xx(flist))]\Di(flist,i);
end

figure('unit','inch','position',[0 0 15 6]);
for i = 1:5
    xx = t/tri(i);
    subplot(2,5,i)
    hold on;
    plot(xx,Di(:,i)/Dinfi(i)-1,'b-','linewidth',1);
    plot(xx,2*Ai(i)/Dinfi(i)./sqrt(t),'b:','linewidth',1);
    set(gca,'xscale','log','yscale','log');
    xlim([1e-1 250]); ylim([1e-3 1e0]);
    set(gca,'fontsize',10,'xtick',10.^[-1:2]);
    xlabel('$t/\tau_r$','interpreter','latex','fontsize',14);
    ylabel('$(D(t)-D_\infty)/D_\infty$','interpreter','latex','fontsize',14);
    box on; grid on;
    title(sprintf('$L=%u \\mu$m',round(Lt(i))),'interpreter','latex','fontsize',14);
    
    subplot(2,5,i+5)
    hold on;
    plot(xx,Ki(:,i),'r-','linewidth',1);
    plot(xx,4*Ai(i)/Dinfi(i)./sqrt(t),'r:','linewidth',1);
    set(gca,'xscale','log','yscale','log');
    xlim([1e-1 250]); ylim([1e-3 1e0]);
    set(gca,'fontsize',10,'xtick',10.^[-1:2]);
    xlabel('$t/\tau_r$','interpreter','latex','fontsize',14);
    ylabel('$K(t)$','interpreter','latex','fontsize',14);
    box on; grid on;
end

%%

Dinf = sum(fi.*Dinfi);
A = sum(Ai.*fi);
Kinf = 3*(sum(fi.*Dinfi.^2)-Dinf^2)/Dinf^2;
X(1) = Dinf;
X(2) = 2*A/sqrt(tr);
Xk(1) = Kinf;
Xk(2) = 4*A/Dinf/sqrt(tr);

figure('unit','inch','position',[0 0 15 5]);
h1 = subplot(131); p1 = h1.Position;
h2 = subplot(132); p2 = h2.Position;
h3 = subplot(133); p3 = h3.Position;
close all

figure('unit','inch','position',[0 0 15 5]);
axes('position',p1);
hold on;
plot(xx,D,'b-','linewidth',1);
xlim([0 250]); ylim([0.9 2]);
set(gca,'fontsize',12,'xtick',0:50:1000,'ytick',[0.9 1:0.2:2]);
xlabel('$t/\tau_r$','interpreter','latex','fontsize',20);
ylabel('$D(t)$, $\mu$m$^2$/ms','interpreter','latex','fontsize',20);
pbaspect([1 1 1]); box on; grid on;

p1i = [0.18 0.39 0.16 0.43];
h1i = axes('Position',p1i);
hold on;
plot(xx,D/X(1)-1,'b-','linewidth',1);
plot(xx,X(2)/X(1)*1./sqrt(xx),'b:','linewidth',1);
set(h1i,'xscale','log','yscale','log');
xlim([1e-1 250]); ylim([1e-3 1e0]);
set(gca,'fontsize',10,'xtick',10.^[-1:2]);
xlabel('$t/\tau_r$','interpreter','latex','fontsize',14);
ylabel('$(D(t)-D_\infty)/D_\infty$','interpreter','latex','fontsize',14);
box on; grid on;

axes('position',p2);
hold on;
plot(xx,K,'r-','linewidth',1);
xlim([0 250]); ylim([0 0.5]);
set(gca,'fontsize',12,'xtick',0:50:1000,'ytick',0:0.1:1);
xlabel('$t/\tau_r$','interpreter','latex','fontsize',20);
ylabel('$K(t)$','interpreter','latex','fontsize',20);
pbaspect([1 1 1]); box on; grid on;

p2i = [0.46 0.39 0.16 0.43];
h2i = axes('Position',p2i);
hold on;
plot(xx,K-Xk(1),'r-','linewidth',1);
% plot(xx,2*X(2)/X(1)*1./sqrt(xx),'r:','linewidth',1);
plot(xx,Xk(2)./sqrt(xx),'r:','linewidth',1);
set(h2i,'xscale','log','yscale','log');
xlim([1e-1 250]); ylim([1e-3 1e0]);
set(gca,'fontsize',10,'xtick',10.^[-1:2]);
xlabel('$t/\tau_r$','interpreter','latex','fontsize',14);
ylabel('$K(t)-K_\infty$','interpreter','latex','fontsize',14);
box on; grid on;

axes('position',p3);
hold on;
tlist = 2:10000;
hd = plot(1./sqrt(xx(tlist)),D(tlist)/X(1)-1,'b-','linewidth',1);
hk = plot(1./sqrt(xx(tlist)),K(tlist)-Xk(1),'r-','linewidth',1);
xlim([0 2.5]); ylim([0 0.2]);
hdr = refline(X(2)/X(1),0); set(hdr,'color','b','linestyle',':','linewidth',1);
% hkr = refline(2*X(2)/X(1),0); set(hkr,'color','r','linestyle',':','linewidth',1);
hkr2 = refline(Xk,0); set(hkr2,'color','r','linestyle',':','linewidth',1);
xlim([0 2.5]); ylim([0 0.5]);
set(gca,'xtick',0:0.5:10,'ytick',0:0.1:10,'fontsize',12);
xlabel('$\sqrt{\tau_r/t}$','interpreter','latex','fontsize',20);
ylabel('Diffusion Metrics','interpreter','latex','fontsize',20);
legend([hd hk],{'$(D(t)-D_\infty)/D_\infty$','$K(t)-K_\infty$'},'interpreter','latex','fontsize',16,'location','northwest');
pbaspect([1 1 1]); grid on; box on;

%%
SV = sum(fi.*2./abar);
Dm = D0*(1-SV*(4/3*sqrt(D0*t/pi)));
figure; plot(1./sqrt(t),D,'-b'); hold on;
% plot(t(1:10),Dm(1:10),'--r');
plot(1./sqrt(t),Dinf+2*A./sqrt(t),'--r');

%% Plot D and K from cumulant and signal

figure('unit','inch','position',[0 0 15 10])
Ax = [bval bval.^2 bval.^3 bval.^4];

Di = zeros(3,1); cDi = zeros(3,1); cKi = zeros(3,1);
[~,It] = min(abs(xx-4));
flist = It:10000;
At = [ones(numel(flist),1) 1./sqrt(xx(flist))];
for i = 1:3
    X = Ax(:,1:(i+1))\log(sig.');
    X = X.';
    Ds = -X(:,1);
    Ws = X(:,2);
    Ks = 6*Ws./Ds.^2;
    
    Xd = At\Ds(flist);
    Xk = At(:,2)\Ks(flist);
    Di(i) = Xd(1);
    cDi(i) = Xd(2);
    cKi(i) = Xk;
    
    subplot(2,3,i);
    hold on;
    hc = plot(1./sqrt(xx(tlist)),D(tlist),'-','linewidth',1); 
    hs = plot(1./sqrt(xx(tlist)),Ds(tlist),'--','linewidth',1);
    set(gca,'fontsize',12,'xtick',0:0.5:10,'ytick',[0.8 1:0.2:2]);
    pbaspect([1 1 1]); box on; grid on;
    xlabel('$\sqrt{\tau_r/t}$','interpreter','latex','fontsize',20);
    ylabel('$D(t)$, $\mu$m$^2$/ms','interpreter','latex','fontsize',20);
    legend([hc hs],{'Cumulant',sprintf('Signal fit up to $b^%u$',i+1)},'interpreter','latex','fontsize',18);
    xlim([0 2.5]); ylim([0.8 1.8]);

    subplot(2,3,i+3);
    hold on;
    hc = plot(1./sqrt(xx(tlist)),K(tlist),'-','linewidth',1);
    hs = plot(1./sqrt(xx(tlist)),Ks(tlist),'--','linewidth',1);
    set(gca,'fontsize',12,'xtick',0:0.5:10,'ytick',0:0.1:1);
    xlabel('$\sqrt{\tau_r/t}$','interpreter','latex','fontsize',20);
    ylabel('$K(t)$','interpreter','latex','fontsize',20);
    pbaspect([1 1 1]); box on; grid on;
    legend([hc hs],{'Cumulant',sprintf('Signal fit up to $b^%u$',i+1)},'interpreter','latex','fontsize',18,'location','southeast');
    xlim([0 2.5]); ylim([0 0.5]);
end

Di_err = Di/Dinf-1;
cDi_err = cDi/(2*A/sqrt(tr))-1
cKi_err = cKi/(4*A/Dinf/sqrt(tr))-1


%% analyze packing

rootpck = '/Users/hhl309/Documents/GitHub/monte-carlo-1D-membrane/packing';

for i = 1:5
    xm = load(fullfile(rootpck,sprintf('membrane_%u',i),'phantom_xMem.txt'));
    L = load(fullfile(rootpck,sprintf('membrane_%u',i),'phantom_res.txt'));
    ai = diff(xm);
    3*var(ai.^2)/mean(ai.^2)^2
end




