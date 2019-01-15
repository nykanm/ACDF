function vol = voxelread()
    [fbin,pbin,~] = uigetfile('*.bin');
    [fmhd,pmhd,~] = uigetfile('*.mhd');
    
    filebin = [pbin fbin];
    filemhd = [pmhd fmhd];
    
    f = fopen(filebin);
    data = fread(f);
    
    mhd = textscan(fopen(filemhd), '%*s %*s %f %f %f', 2);
    ndim = [mhd{1}(2), mhd{2}(2), mhd{3}(2)];
    
    vol = reshape(data, ndim);
    volumeViewer(vol);
end