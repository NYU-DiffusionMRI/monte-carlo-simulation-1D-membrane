% In this example, we create five 1d micro-geometries consisting of
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

%% Create 1d micro-geometry consisting of permeable membranes

rng(0);
ncat = 5;                           % Five micro-geometries
for i = 1:ncat
    rooti = fullfile(target,sprintf('membrane_%u',i));
    mkdir(rooti);
    at = [];                        % Membrane distance of all axons
    for j = 1:numel(fiber)
        ai = diff(fiber(j).pos(:,1));
        at = cat(1,at,ai);
    end
    at = at(randperm(numel(at)));   % Randomly shffule the membrane distance
    xm = cumsum(at);                % Membrane position
    L = xm(end);                    % Axonal length
    xm = xm-xm(1)/2;
    xm = xm/L;                      % Membrane position normalized with the axonal length
    
    % Create files for simulations
    fid = fopen(fullfile(rooti,'phantom_xMem.txt'),'w');    % Membrane position
    fprintf(fid,sprintf('%.8f\n',xm));
    fclose(fid);
    
    fid = fopen(fullfile(rooti,'phantom_NMem.txt'),'w');    % Membrane number
    fprintf(fid,sprintf('%u\n',numel(xm)));
    fclose(fid);
    
    fid = fopen(fullfile(rooti,'phantom_res.txt'),'w');     % Axonal length
    fprintf(fid,sprintf('%.8f\n',L));
    fclose(fid);
end

%% Create lookup table for 1d micro-geometry

ncat = 5;                           % Five micro-geometries
NPix = 1e4;                         % Matrix size of the lookup table
for i = 1:ncat
    rooti = fullfile(target,sprintf('membrane_%u',i));
    xm = load(fullfile(rooti,'phantom_xMem.txt'));  % Load membrane position
    at = [xm(1)*2; diff(xm)];       % Membrane distance
    
    % Check the matrix size of the lookup table
    L = load(fullfile(rooti,'phantom_res.txt'));    % Load axonal length
    D0 = 2;                         % Intrinsic diffusivity, micron^2/ms
    dt = 2e-3;                      % Time-step, ms
    dx = sqrt(2*D0*2e-3);           % Step size, micron
    NPix_max = ceil(L/dx);          % Upper bound of the matrix size
    fprintf('Membrane_%u: NPix shoudl be less than %d.\n',i,NPix_max);
    
    A = zeros(NPix,1);              % Lookup table
    am = ceil(xm*NPix);             % Membrane position in lookup table
    if numel(unique(am)) == numel(am)
        A(am) = 1:numel(xm);
    else
        fprintf('Use larger NPix for Membrane_%u.\n',i)
    end
    
    % Create files for simulations
    fid = fopen(fullfile(rooti,'phantom_APix.txt'),'w');    % Lookup table
    fprintf(fid,sprintf('%u\n',A));
    fclose(fid);
    
    fid = fopen(fullfile(rooti,'phantom_NPix.txt'),'w');    % Matrix size
    fprintf(fid,sprintf('%u\n',NPix));
    fclose(fid);
end
