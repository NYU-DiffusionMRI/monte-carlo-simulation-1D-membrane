% In this example, we perform Monte Carlo simualtions in five 1d 
% micro-geometries consisting of permeable membranes (demo2_packing_file.m)
%
% Author: Hong-Hsi Lee, December, 2019 (orcid.org/0000-0002-3663-6559)

% set up the root and target directories on your own computer
% root = '/directory/to/this/file';
root = '.';
root_input = fullfile(root,'hpc_code/input');
root_lib = fullfile(root,'hpc_code/lib');

%% Create files for simualtion

% Simulation parameters
dt = 2e-3;                          % time-step, ms
TN = 7e5;                           % # steps
NPar = 1e8;                         % # random walkers
Din = 2;                            % Intrinsic diffusivity, micron^2/ms
kappa = 0.4154;                     % Membrane permeability, micron/ms
pinit = 1;                          % Initial position, 1=everywhere, 2=center
threadpb = 512;                     % Thread-per-block in CUDA
fileID = fopen(fullfile(root_input,'simParamInput.txt'),'w');
fprintf(fileID,sprintf('%g\n%u\n%u\n%g\n%g\n%u\n%u\n',dt,TN,NPar,Din,kappa,pinit,threadpb));
fclose(fileID);

% Create b-table
Nbval = 4;                          % # b-value
bval = [0.1 0.4 1 1.5].';           % b-values, ms/micron^2
fileID = fopen(fullfile(root_input,'Nbval.txt'),'w');
fprintf(fileID,sprintf('%u\n',Nbval));
fclose(fileID);

fileID = fopen(fullfile(root_input,'bval.txt'),'w');
fprintf(fileID,sprintf('%g\n',bval));
fclose(fileID);

% Copyfile to each micro-geometry
files = dir(fullfile(root_input,'membrane_*'));
for i = 1:numel(files)
    copyfile(fullfile(root_input,'simParamInput.txt'),fullfile(root_input,files(i).name,'simParamInput.txt'));
    copyfile(fullfile(root_input,'Nbval.txt'),fullfile(root_input,files(i).name,'Nbval.txt'));
    copyfile(fullfile(root_input,'bval.txt'),fullfile(root_input,files(i).name,'bval.txt'));
end

%% Run simulations in CUDA C++, you must have an Nvidia GPU to run the code

files = dir(fullfile(root_input,'membrane_*'));
for i = 1:numel(files)
    target = fullfile(root_input,files(i).name);
    
    % The directory to nvcc might be different on your own computer.
    % It might be easier to compile the code in the terminal.
    system(sprintf('/usr/local/cuda/bin/nvcc -arch=sm_70 -rdc=true %s/main.cu -o %s/my_code',root_lib,target));
    system(sprintf('%s/my_code',target));
end

%% Exercise: Re-write the CUDA C++ code to a C++ code.
