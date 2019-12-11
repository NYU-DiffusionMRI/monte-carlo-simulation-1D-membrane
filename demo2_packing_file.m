clear
root = '/Users/hhl309/Documents/GitHub/monte-carlo-1D-membrane/packing';
target = fullfile(root,'membrane_v4'); mkdir(target)
load(fullfile(root,'Hellwig2014_bead.mat'));

%% create packing based on histology
ncat = 5;
Lt = [100 400 1000 2000 3000];
for i = 1:ncat
    rng(0);
    rooti = fullfile(target,sprintf('membrane_%u',i));
    mkdir(rooti);
    at = [];
    for j = 1:numel(fiber)
        ai = diff(fiber(j).pos(:,1));
        at = cat(1,at,ai);
    end
    at = at(randperm(numel(at)));
    xm = cumsum(at);
    [~,If] = min(abs(xm-Lt(i)));
    xm = xm(1:If);
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
for i = 1:ncat
    rooti = fullfile(target,sprintf('membrane_%u',i));
    xm = load(fullfile(rooti,'phantom_xMem.txt'));
    
    L = load(fullfile(rooti,'phantom_res.txt'));
    D0 = 2; dt = 2e-3;
    dx = sqrt(2*D0*2e-3);
    NPix_min = 1/min(xm);
    NPix_max = L/dx;
    NPix = round(NPix_min/2+NPix_max/2);
    fprintf('Set NPix between %d and %d.\n',floor(NPix_min),ceil(NPix_max));
    
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








