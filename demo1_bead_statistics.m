% In this example, we calculate (1) the power spectrum of restrictions
% along axons based on the position of axonal beads in Fig. 6 in Hellwig et
% al., 1994 (https://doi.org/10.1007/BF00198906) and (2) the mean and 
% variance of the bead number within a sliding window.
%
% Author: Hong-Hsi Lee, December, 2019 (orcid.org/0000-0002-3663-6559)

% set up the root directory on your own computer
% root = '/directory/to/this/file';
root = '.';

%% Data captured from Fig. 6 in Hellwig, et al., 1994

% fiber.pos(:,1), The position of axonal beads on the 33 primay collaterals
file = load(fullfile(root,'hpc_code/input/Hellwig2014_bead.mat'));
fiber = file.fiber;

figure('unit','inch','position',[0 0 10 8]); hold on;
for i = 1:numel(fiber)
    plot(fiber(i).pos(:,1),fiber(i).pos(:,2),'.','markersize',8);
end
box on;
set(gca,'fontsize',12,'xtick',0:100:400,'ytick',1:33)
xlabel('length ($\mu$m)','interpreter','latex','fontsize',20);
ylabel('collateral no.','interpreter','latex','fontsize',20);
title('Fig. 6 in Hellwig et al., 1994','interpreter','latex','fontsize',20)
axis ij;
ylim([0 34])

%% Calculate the Power Spectrum of restrictions based on the bead position

% Step 1: Concatenated axon

% fiber.pos(:,1), The position of axonal beads on the 33 primay collaterals
file = load(fullfile(root,'hpc_code/input/Hellwig2014_bead.mat'));
fiber = file.fiber;

% Calculate the theoretical value of plateau of power spectrum
Lt = zeros(numel(fiber),1);         % Axonal length for each axon
plateau = zeros(numel(fiber),1);    % Plateau of the power spectrum for each axon
at = [];                            % Bead distance of all axons
for i = 1:numel(fiber)
    Lt(i) = range(fiber(i).pos(:,1));
    ai = diff(fiber(i).pos(:,1));
    plateau(i) = var(ai)/mean(ai)^3;
    at = cat(1,at,ai);
end
Lmin = min(Lt);                     % Axonal length of the shortest axon
abar = mean(at);                    % Mean bead distance
plateau = plateau * abar;           % Plateau of the power spectrum normalized with the mean bead distance

% Calculate the power spectrum of concatenated axon
rng(0);
nite = 200;                         % Randomly concatenate axons 200 times
Gammak = 0;                         % Power spectrum of restrictions
% Here, we choose the (shortest bead distance)/37 as the length unit to
% digitize the bead position of concatenated axon
vox = min(at)/37;
for i = 1:nite
    I = randperm(numel(fiber));     % Randomly shuffle the axon order
    ai = [];
    for j = I
        aj = diff(fiber(I(j)).pos(:,1));
        ai = cat(1,ai,aj);
    end
    pos = cumsum(ai);               % Bead position of concatenated axon
    pos = round(pos/vox);           % Digitize the bead position
    rho = false(round(max(pos)),1); % Concatenated axon
    rho(pos) = 1;
    rhok = fft(rho);                % FT of concatenated axon
    Gammak = Gammak + abs(rhok).^2/(length(rho)*vox);   % Power spectrum
end
Gammak = Gammak/nite;               % Power spectrum averaged over 200 times
L = size(Gammak,1)*vox;             % Axonal length of concatenated axon
kL = 0:(size(Gammak,1)-1);          % k*L/(2*pi)
ka = kL/L*abar;                     % k*abar/(2*pi)

% The smallest k*abar/(2*pi) is limited by the length of the shortest axon
ka_min = abar/Lmin;
Ik = find(ka>ka_min,1,'first');

% Plot figure
figure('unit','inch','position',[0 0 5 5]);
hold on;
ha = plot(ka(Ik:end),Gammak(Ik:end)*abar,'r-');
set(gca,'xscale','log','yscale','log');
hr = refline(0,mean(plateau));      % Theoretical value of plateau at low k
set(hr,'color','r','linestyle','--','linewidth',1);
set(gca,'fontsize',12);
xlabel('$k\bar{a}/2\pi$','interpreter','latex','fontsize',20);
ylabel('$\Gamma(k)\cdot\bar{a}$','interpreter','latex','fontsize',20);
xlim([1e-2 1e1]);
ylim([1e-1 1e1]);
box on; grid on;

% Step 2: Poisson distributed beads

% Calculate the power spectrum of Poisson distributed beads
rng(0);
nite = 200;                         % Randomly generate axons 200 times
N = 1e5;                            % The matrix size of axon
Gammak = zeros(N,1);                % Power spectrum
abar = 317;                         % Mean bead distance
for i = 1:nite
    pos = rand(N,1)<(1/abar);       % Bead position
    rho = zeros(N,1);               % Axon with Poisson distributed beads
    rho(pos) = 1;
    rhok = fft(rho);                % FT of the axon
    Gammak = Gammak + abs(rhok).^2/N;   % Power spectrum
end
Gammak = Gammak/nite;               % Power spectrum averaged over 200 times
L = N*1;                            % Axonal length
kL = 0:(N-1);                       % k*L/(2*pi)
ka = kL/L*abar;                     % k*abar/(2*pi)

% Plot figure
hp = plot(ka,Gammak*abar,'b-');
set(gca,'xscale','log','yscale','log');
hr = refline(0,1);
set(hr,'color','b','linestyle','--','linewidth',1);
box on; grid on;
legend([ha hp],{'Hellwig et al., 1994','Poisson'},'interpreter','latex','fontsize',20);
pbaspect([1 1 1]);

%% Calculate the mean and variance of the bead number within a sliding window

% Step 1: Concatenated axon

% fiber.pos(:,1), The position of axonal beads on the 33 primay collaterals
file = load(fullfile(root,'hpc_code/input/Hellwig2014_bead.mat'));
fiber = file.fiber;

rng(0);
% Here, we choose the (shortest bead distance)/37 as the length unit to
% digitize the bead position of concatenated axon
vox = min(at)/37;
sw = struct([]);                    % Sliding window result
for i = 1:numel(fiber)
    pos = fiber(i).pos(:,1);        % Bead position of the i-th axon
    pos = round(pos/vox);           % Digitize the bead position
    rho = false(round(max(pos)),1); % The i-th axon
    rho(pos) = 1;
    nmeani = [];                    % Mean of bead number in the window
    nvari = [];                     % Variance of bead number in the window
    % Here, we impose the sliding window in a width between 3*vox to L/10.
    for j = 3:40:ceil(numel(rho)/10)
        nj = conv(rho,ones(j,1),'valid');
        % We choose only 5% sliding windows to calculate mean and variance
        % of bead number to avoid correlations between chosen windows
        nj = nj(randperm(numel(nj),ceil(numel(rho)*0.05)));
        nmeani = cat(1,nmeani,mean(nj));
        nvari = cat(1,nvari,var(nj));
    end
    [nmeani,Ii] = sort(nmeani);
    sw(i).nmean = nmeani;
    sw(i).nvar = nvari(Ii);
end

% Step 2: Poisson distributed beads

rng(0);
nite = 25;                          % Randomly generate axons 25 times
N = 1e5;                            % The matrix size of axon
abar = 317;                         % Mean bead distance
swp = struct([]);                   % Sliding window result
swp(1).nmean = 0; swp(1).nvar = 0;  % Initialization
for i = 1:nite
    pos = rand(N,1)<(1/abar);       % Bead position
    rho = zeros(N,1);               % Axon with Poisson distributed beads
    rho(pos) = 1;
    nmeani = [];                    % Mean of bead number in the window
    nvari = [];                     % Variance of bead number in the window
    % Here, we impose the sliding window in a width between 3*vox to L/25.
    for j = 3:100:ceil(numel(rho)/25)
        nj = conv(rho,ones(j,1),'valid');
        % We choose only 2% sliding windows to calculate mean and variance
        % of bead number to avoid correlations between chosen windows
        nj = nj(randperm(numel(nj),ceil(numel(rho)*0.02)));
        nmeani = cat(1,nmeani,mean(nj));
        nvari = cat(1,nvari,var(nj));
    end
    swp.nmean = swp.nmean + nmeani;
    swp.nvar = swp.nvar + nvari;
end
% Average over 25 times
swp.nmean = swp.nmean/nite;
swp.nvar = swp.nvar/nite;

% Plot figure
figure('unit','inch','position',[0 0 5 5]);
hold on;
slp = zeros(numel(fiber),1);
list = 1:3;
% Calculate the slope of <N> vs <N^2> - <N>^2 for each axon
for i = 1:numel(fiber)
    ha = plot(sw(i).nmean,sw(i).nvar,'r-');
    slp(i) = sw(i).nmean(list)\sw(i).nvar(list);
end
hp = plot(swp.nmean,swp.nvar,'b.-','linewidth',1,'markersize',10);
hr = refline(1,0); set(hr,'color','k','linewidth',0.5);
hr2 = refline(mean(slp),0); set(hr2,'color','r','linewidth',1,'linestyle','--');
xlim([0 12]); ylim([0 12])

box on; grid on;
set(gca,'xtick',0:2:20,'ytick',0:2:20,'fontsize',12);
xlabel('$\langle N\rangle$','interpreter','latex','fontsize',20);
ylabel('$\langle N^2\rangle - \langle N\rangle^2$','interpreter','latex','fontsize',20);
legend([ha hp],{'Hellwig et al., 1994','Poisson'},'interpreter','latex','fontsize',20);
pbaspect([1 1 1]);

