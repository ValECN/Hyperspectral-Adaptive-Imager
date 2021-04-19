function [FCS, I, panchro, DMD_conf, IC] = rebuild_aqucube_Val(path_directory)

addpath(path_directory);
addpath(genpath(path_directory));

original_files = dir([path_directory '/*.h5']); 
filename_ex = [path_directory '/' original_files(1).name];
S = length(original_files);

IC = permute(h5read(filename_ex,'/IC'),[3 2 1]);

% Pre-allocating space
FCS = zeros([size(h5read(filename_ex,'/FC')),S]);
FCS = permute(FCS,[3 2 1 4]);
FCS = permute(FCS,[2 1 3 4]);

I = zeros([size(h5read(filename_ex,'/I')),S]); 

DMD_conf = zeros([size(h5read(filename_ex,'/DMD')),S]); 
DMD_conf = permute(DMD_conf, [2 1 3]);

panchro = zeros([size(h5read(filename_ex,'/panchro')),S]);

for k = 1:S
    
    filename        = [path_directory '/' original_files(k).name]; 
    FCS(:,:,:,k)    = permute(h5read(filename,'/FC'),[2 3 1]);
    I(:,:,k)        = h5read(filename,'/I');
    DMD_conf(:,:,k) = permute(h5read(filename,'/DMD'),[2,1]);
    panchro(:,:,k)  = h5read(filename,'/panchro');
    
end
end