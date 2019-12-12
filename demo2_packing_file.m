% In this example, we will create five 1d micro-geometries consisting of
% permeable membranes, whose position is designed based on the axonal bead
% distance in Fig. 6 in Hellwig et al., 1994 
% (https://doi.org/10.1007/BF00198906).
%
% Author: Hong-Hsi Lee, December, 2019 (orcid.org/0000-0002-3663-6559)

% set up the root and target directories on your own computer
% root = '/directory/to/this/file';
root = '.';
target = fullfile(root,'hpc_code/input');

% fiber.pos(:,1), The position of axonal beads on the 33 primay collaterals
file = load(fullfile(root,'hpc_code/input/Hellwig2014_bead.mat'));
fiber = file.fiber;

%% Create 1d microstructure consisting of permeable membranes

rng(0);
ncat = 5;                           % Five micro-geometries
for i = 1:ncat
    rooti = fullfile(target,sprintf('membrane_%u',i));
    mkdir(rooti);
    at = [];
    for j = 1:numel(fiber)
        ai = diff(fiber(j).pos(:,1));
        at = cat(1,at,ai);
    end
    at = at(randperm(numel(at)));
    xm = cumsum(at);
    L = xm(end);
    xm = xm-xm(1)/2;
    xm = xm/L;
    fid = fopen(fullfile(rooti,'phantom_xMem.txt'),'w');
    fprintf(fid,sprintf('%.8f\n',xm));
    fclose(fid);
    
    fid = fopen(fullfile(rooti,'phantom_NMem.txt'),'w');
    fprintf(fid,sprintf('%u\n',numel(xm)));
    fclose(fid);
    
    fid = fopen(fullfile(rooti,'phantom_res.txt'),'w');
    fprintf(fid,sprintf('%.8f\n',L));
    fclose(fid);
end

%% create lookup table
ncat = 5;
NPix = 1e4;
for i = 1:ncat
    rooti = fullfile(target,sprintf('membrane_%u',i));
    xm = load(fullfile(rooti,'phantom_xMem.txt'));
    
    A = zeros(NPix,1);
    am = ceil(xm*NPix);
    if numel(unique(am)) == numel(am)
        A(am) = 1:numel(xm);
    else
        fprintf('Use larger NPix for fiber %u\n',i)
    end
    fid = fopen(fullfile(rooti,'phantom_APix.txt'),'w');
    fprintf(fid,sprintf('%u\n',A));
    fclose(fid);
    
    fid = fopen(fullfile(rooti,'phantom_NPix.txt'),'w');
    fprintf(fid,sprintf('%u\n',NPix));
    fclose(fid);
end

%% Plot power spectrum
ps = 0;
pl = [];
for i = 1:5
    rooti = fullfile(target,sprintf('membrane_%u',i));
    mi = load(fullfile(rooti,'phantom_APix.txt'));
    ai = diff(find(mi)); pl(i) = var(ai)/mean(ai)^3;
    km = fft(logical(mi));
    ps = ps+abs(km).^2;
end
ps = ps/5/numel(ps);
xx = 0:(numel(ps)-1);
figure; plot(xx,ps); set(gca,'xscale','log','yscale','log')
refline(0,mean(pl));








