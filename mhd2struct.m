%% Change the following parameters for each step
loadfolder = 'step5';
savename = 'voxel_step5.mat';

%% Change the following to your directories on your computer
loaddir = ['/Users/nykan/Documents/McGill Grad/Matlab/ACDF/Data/SceneObjects/' loadfolder];
savedir = ['/Users/nykan/Documents/McGill Grad/Matlab/ACDF/Data/SceneObjects/' savename];

%% This converts the voxels into a structure and saves it as a .mat file
subdir = dir(loaddir);

for io=[3:2:length(subdir);1:9]
    i=io(1);
    o=io(2);
    filename_bin = [loaddir '/' subdir(i).name];
    filename_mhd = [loaddir '/' subdir(i+1).name];
    
    f = fopen(filename_bin);
    data = fread(f);
    
    mhd = textscan(fopen(filename_mhd),'%*s %*s %f %f %f', 2);
    ndim = [mhd{1}(2), mhd{2}(2), mhd{3}(2)];
    
    name = compose(subdir(i+1).name(16:end-4));
    name = strrep(name{1},'-','_');
    vol = reshape(data, ndim);
    voxel.(name) = vol;
end

save(savedir,'voxel');