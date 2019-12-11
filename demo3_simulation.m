% This example demonstrates the short-time limit of diffusion in
% extra-cylindrical space of randomly packed impermeable cylinders in 2d,
% shown in Figure 4, point 2 in (Fieremans and Lee, NeuroImage 2018), with
% more details in supplementary information.
%
% Author: Hong-Hsi Lee, September, 2018 (orcid.org/0000-0002-3663-6559)

% ********** Setup the directory on your computer **********
root = '/Users/hhl309/Documents/GitHub/monte-carlo-1D-membrane/hpc_code/input_cuda';
rootpck = '/Users/hhl309/Documents/GitHub/monte-carlo-1D-membrane/packing/membrane_v4';
target = fullfile(root,'membrane_v4'); mkdir(target)

% fileID = fopen(fullfile(target,'root.txt'),'w');
% fprintf(fileID,'%s',target);
% fclose(fileID);

% Create files for simualtion

% Simulation parameters
dt = 2e-3;  % time of each step, ms
TN = 7e5;   % # steps
NPar = 1e7; % # random walkers
Din = 2;    % IAS diffusivity, µm^2/ms
kappa = 0.4154;% permeability, µm/ms
pinit = 1;  % initial position, 1=ICS, 2=center
threadpb = 512;  % Thread per block for cuda
fileID = fopen(fullfile(target,'simParamInput.txt'),'w');
fprintf(fileID,sprintf('%g\n%u\n%u\n%g\n%g\n%u\n%u\n',dt,TN,NPar,Din,kappa,pinit,threadpb));
fclose(fileID);

Nbval = 4;
bval = [0.1 0.4 1 1.5].';
fileID = fopen(fullfile(target,'Nbval.txt'),'w');
fprintf(fileID,sprintf('%u\n',Nbval));
fclose(fileID);

fileID = fopen(fullfile(target,'bval.txt'),'w');
fprintf(fileID,sprintf('%g\n',bval));
fclose(fileID);

files = dir(fullfile(rootpck,'membrane_*'));
for i = 1:5
    copyfile(fullfile(rootpck,files(i).name),fullfile(target,files(i).name));
    copyfile(fullfile(target,'simParamInput.txt'),fullfile(target,files(i).name,'simParamInput.txt'));
    copyfile(fullfile(target,'Nbval.txt'),fullfile(target,files(i).name,'Nbval.txt'));
    copyfile(fullfile(target,'bval.txt'),fullfile(target,files(i).name,'bval.txt'));
end

%% Have a look for the microstructure
% The lookup table A saves two axon labels in one integer. If the first
% and the second axon labels are ax1 and ax2, ax1 = mod(A,Nmax), and ax2 =
% floor(A/Nmax).
% Other parameters:
%   fov: field of view of the entire gemoetry and lookup table A in µm
%   Nax: # axons
%   rCir: axon's outer radius
%   gratio: g-ratio, the ratio of inner to outer radius
%   [xCir,yCir]: axon's center position
root_packing = '/Users/hhl309/Documents/GitHub/monte-carlo-3D-coaxial-cylinder/packing/sphere_1';

A = load(fullfile(root_packing,'phantom_APix.txt'));
Nmax = load(fullfile(root_packing,'phantom_Nmax.txt'));
NPix = load(fullfile(root_packing,'phantom_NPix.txt'));
fov = load(fullfile(root_packing,'phantom_res.txt'));

Nax = load(fullfile(root_packing,'phantom_NAx.txt'));
% gratio = load(fullfile(root_packing,'phantom_gratio.txt')); 
rCir = load(fullfile(root_packing,'phantom_rCir.txt'));
xCir = load(fullfile(root_packing,'phantom_xCir.txt'));
yCir = load(fullfile(root_packing,'phantom_yCir.txt'));
zCir = load(fullfile(root_packing,'phantom_zCir.txt'));

%% Plot the microstructure and the lookup table
% Plot the microstructure
[xs, ys, zs] = sphere;

% show initial configuration
figure; set(gcf,'unit','inch','position',[0 0 12 5])
subplot(121);
hold on;
for i = 1:numel(rs)
    surf(xs*rCir(i)+xCir(i), ys*rCir(i)+yCir(i), zs*rCir(i)+zCir(i));
end
xlim([0 1]); ylim([0 1]); zlim([0 1]);
pbaspect([1 1 1]); box on
title('Axon Packing','interpreter','latex','fontsize',20)
set(gca,'xtick',[],'ytick',[])

% Plot the lookup table
subplot(122);
A2 = zeros(NPix,NPix,NPix);
for i = 1:NPix
    for j = 1:NPix
        A2(i,j,:) = A(NPix*(i-1)+j,:);
    end
end
imagesc(rot90(A2(:,:,50)));
% cmap = colormap('parula');
% Ibg = A==0;                     % background region
% Iol = A>Nax;                    % two-axon region
% A2 = ceil(single(A)/Nax*64);    % rescale the colormap for # axons
% A2(Ibg) = 1;
% A2(Iol) = 1;
% [nx,ny,nz] = size(A2);
% imgc = cmap(uint16(A2(:)),:);
% imgc(Ibg,:) = 0;                % background region is black
% imgc(Iol,:) = 1;                % two-axon region is white
% imgc = reshape(imgc,[nx,ny,nz,3]);
% 
% image(rot90(squeeze(A(:,:,end)))); %caxis([0 Nax]);
% box on; axis off
% title('Lookup Table','interpreter','latex','fontsize',20)

%% Run the simulation in the extra-cylindrical space in 2d
% Parameters are defined in the 'simParamInput.txt', and the # particle is
% defined in main.cpp, line 50.
clear
clc
root_code = '/Users/hhl309/Documents/GitHub/monte-carlo-3D-sphere/hpc_code/lib';
cd(root_code)
system('g++ -std=c++11 RNG.cpp diffusion_lib.cpp main.cpp -o my_code')
system('./my_code')

%% Read simulation results

% load geometry
root_packing = '/Users/hhl309/Documents/GitHub/monte-carlo-3D-sphere/input/sphere_2';
rc = load(fullfile(root_packing,'phantom_rCir.txt'));
fin = 4/3*pi*sum(rc.^3);

% read simulation parameters
fileID = fopen(fullfile(root_code,'sim_para.txt'),'r');
n=1; C = {};
tline=fgetl(fileID);
while ischar(tline)
    C{n}=tline;
    tline=fgetl(fileID);
    n=n+1;
end
fclose(fileID);

dt = str2double(C{1});        % time step
TN = str2double(C{2});        % # step
Np = str2double(C{3});        % # particles
Din = str2double(C{4});       % IAS diffusivity
Dex = str2double(C{5});      % EAS diffusivity
kappa = str2double(C{6});     % permeability

% read simulation result
dx2 = load(fullfile(root_code,'dx2_diffusion.txt'))/Np;
dx4 = load(fullfile(root_code,'dx4_diffusion.txt'))/Np;

% calculate radial diffusivity (RD) and radial kurtosis (RK)
t = dt*(1:TN/1e3:TN).';
Dx = dx2(:,1)/2./t;
Kx = dx4(:,1)./dx2(:,1).^2 - 3;

dx2mean = mean(dx2(:,[1,4,6]),2);
MD = dx2mean/2./t;
tinst = t(3:end-2);
MDinst = (-dx2mean(5:end)+8*dx2mean(4:end-1)-8*dx2mean(2:end-3)+dx2mean(1:end-4))/12/mean(diff(t))/2;

figure;
subplot(121);
tlist = 501:1000;
temp = [ones(numel(tlist),1),1./t(tlist)]\MD(tlist);
Dinf = temp(1); c=temp(2);
plot(1./t,MD-Dinf,'.');
xlim([0 2]);
hr = refline(c,0); set(hr,'color','r');

subplot(122);
% temp = [ones(numel(tlist),1),t(tlist).^(-3/2)]\MDinst(tlist);
% Dinf = temp(1); c=temp(2);
Ipos = (MDinst>Dinf).*(tinst<0.1); Ipos = logical(Ipos(:));
temp = [ones(nnz(Ipos),1) -log(tinst(Ipos))]\log(MDinst(Ipos)-Dinf);

plot(tinst,MDinst-Dinf); hold on;
% plot(tinst,exp(temp(1))*tinst.^(-temp(2)),'-r');
set(gca,'yscale','log','xscale','log');

% temp(2)


% % read particle
% Nin = load(fullfile(root_code,'NParICS.txt'));
% Nex = load(fullfile(root_code,'NParECS.txt'));
% figure; hold on; plot(t,Nin/Np,'r.'); hin = refline(0,fin); set(hin,'color','r');
% plot(t,Nex/Np,'b.'); hex = refline(0,1-fin); set(hex,'color','b');
% ylim([0 1]);

%% Load diffusion trajectory
xt = load(fullfile(root_code,'x_diffusion.txt'));
yt = load(fullfile(root_code,'y_diffusion.txt'));
zt = load(fullfile(root_code,'z_diffusion.txt'));
xt = single(xt); yt = single(yt); zt = single(zt);
save(fullfile(root_code,'diffusion_trajectory.mat'),'xt','yt','zt');

%% Plot diffusion trajectory
figure; hold on;
% for i = 5%size(xt,2)
%     hi = plot3(xt(1,i),yt(1,i),zt(1,i),'ro');
%     set(hi,'markersize',10);
%     plot3(xt(:,i),yt(:,i),zt(:,i),'.-');
% end
plot3(xt(1000,:),yt(1000,:),zt(1000,:),'.');
xlim([0 1]); ylim([0 1]); zlim([0 1]);
pbaspect([1 1 1]);
view(3);

xc = load(fullfile('/Users/hhl309/Documents/GitHub/monte-carlo-3D-sphere/input/sphere_2','phantom_xCir.txt'));
yc = load(fullfile('/Users/hhl309/Documents/GitHub/monte-carlo-3D-sphere/input/sphere_2','phantom_yCir.txt'));
zc = load(fullfile('/Users/hhl309/Documents/GitHub/monte-carlo-3D-sphere/input/sphere_2','phantom_zCir.txt'));
rc = load(fullfile('/Users/hhl309/Documents/GitHub/monte-carlo-3D-sphere/input/sphere_2','phantom_rCir.txt'));
res = load(fullfile('/Users/hhl309/Documents/GitHub/monte-carlo-3D-sphere/input/sphere_2','phantom_res.txt'));
[xs,ys,zs] = sphere;
ii = [0 2*(xc<1/2)-1];
jj = [0 2*(yc<1/2)-1];
kk = [0 2*(zc<1/2)-1];
for i = 1:2
    for j = 1:2
        for k = 1:2
            h = surf(rc*xs+xc+ii(i), rc*ys+yc+jj(j), rc*zs+zc+kk(k));
            set(h,'facealpha',0.5);
        end
    end
end
xlabel('x');
ylabel('y');
zlabel('z');
%% Calculate permeability by using particle density
rt = sqrt((xt-0.5).^2+(yt-0.5).^2+(zt-0.5).^2);
nbin = 48*2;
edges = linspace(0,0.5,nbin+1);
centers = edges(1:end-1)/2 + edges(2:end)/2;
Nt = zeros(size(rt,1),nbin);
for i = 1:size(rt,1)
    rti = rt(i,:);
    Ni = histcounts(rti,edges);
    Nt(i,:) = Ni./(4*pi*centers.^2);
end
%% Plot permeability over time points
list = nbin/4+1:nbin/4*3;
ti = 200;
figure('unit','inch','position',[0 0 12 6]); subplot(121);
plot(centers(list)*res,Nt(ti,list),'.-');
hold on;
plot([rc*res,rc*res],[0 max(Nt(ti,list))]);
pbaspect([1 1 1]);
xlabel('$x$ ($\mu$m)','interpreter','latex','fontsize',20);
ylabel('$n(x)$ a.u.','interpreter','latex','fontsize',20);
title(sprintf('$t=$%.2f ms',t(ti)),'interpreter','latex','fontsize',20)

kappa_fit = zeros(size(rt,1),1);
npt_in = 5; npt_ex = round(npt_in/sqrt(Din/Dex));
for i = 1:size(rt,1)
    Nti = Nt(i,:); Nti = Nti(:);
    listin = nbin/2-npt_in:nbin/2-1; listin = listin(:);
    temp = [ones(numel(listin),1) listin*mean(diff(centers))*res]\Nti(listin);
%     nin = temp(1);
    Jin = Din*abs(temp(2));
    
    listex = nbin/2+2:nbin/2+npt_ex+1; listex = listex(:);
    temp = [ones(numel(listex),1) listex*mean(diff(centers))*res]\Nti(listex);
%     nex = temp(1);
    Jex = Dex*abs(temp(2));
    
    nin = mean(Nti(listin));
    nex = mean(Nti(listex));
    dn = nin-nex;
    
    Jmean = (Jin+Jex)/2;
    kappa_fit(i) = Jmean/dn;
end

subplot(122); plot(t,kappa_fit)
hr = refline(0,kappa); set(hr,'color','r');
ylim([0 0.6]);
pbaspect([1 1 1]);
xlabel('$t$ (ms)','interpreter','latex','fontsize',20);
ylabel('$\kappa$ ($\mu$m/ms)','interpreter','latex','fontsize',20);
title(['$\kappa=$' sprintf('%.3f',mean(kappa_fit(end/2:end))) ' $\mu$m/ms'],...
    'interpreter','latex','fontsize',20);

%% Plot short-time limit for randomly packed impermeable cylinders (2d)
%  Particles diffuse in extra-axonal space only.

figure; set(gcf,'unit','inch','position',[1 1 14 7])
subplot(121)
hold on
h = plot(sqrt(t),RD/Dout,'-'); axis square
set(h,'linewidth',3)

% calculate surface S, volume V, surfave-to-volume ratio SoV
S = fov*2*pi*sum(rCir);
V = fov^2*(1-pi*sum(rCir.^2));
SoV= S/V;
D_theory = Dout*(1-SoV/2*(4/3)*sqrt(Dout*t)/sqrt(pi));
K_theory = (3/8)*(8/5)*SoV*sqrt(Dout*t)/sqrt(pi);

h_theory = plot(sqrt(t),D_theory/Dout,'--r'); set(h_theory,'linewidth',2);
legend([h,h_theory],{'Simulation','Mitra Limit'},'interpreter','latex','fontsize',30,'location','southwest')

set(gca,'xtick',0:0.05:1,'ytick',0:0.5:1,'fontsize',20);
xlim([0 0.3]); 
ylim([0 1.25]); box on; grid on
xlabel('$\sqrt{t}$ (ms$^{1/2}$)','interpreter','latex','fontsize',30);
ylabel('$D(t)/D_0$','interpreter','latex','fontsize',30);

subplot(122);
hold on
h = plot(sqrt(t),RK,'-'); axis square
set(h,'linewidth',3)
h_theory = plot(sqrt(t),K_theory,'--r');
set(h_theory,'linewidth',2);
hl = refline(0,0); set(hl,'color','k')
legend([h,h_theory],{'Simulation','Mitra limit'},'interpreter','latex','fontsize',30,'location','southeast')

set(gca,'xtick',0:0.05:1,'ytick',-2:0.5:1,'fontsize',20);
xlim([0 0.3]); 
ylim([-2 1]); box on; grid on
xlabel('$\sqrt{t}$ (ms$^{1/2}$)','interpreter','latex','fontsize',30);
ylabel('$K(t)$','interpreter','latex','fontsize',30);
