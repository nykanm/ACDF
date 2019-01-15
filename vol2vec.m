function [xf, yf, zf, c] = vol2vec(vol, res, pos)
% This function is used by voxel_figure.m to visualise voxel structures.
% res is the ElementSpacing (e.g. [0.02 0.02 0.02]) 
% pos is the Position (e.g. [3.42573 -1.56917 -8.35422])
%% DiscRemoval/CervicalC5-Disc_voxels

[xx, yy, zz] = ind2sub(size(vol), 1:numel(vol));
idx = find(vol(:));

xi = xx(idx)';
yi = yy(idx)';
zi = zz(idx)';
c = vol(idx);

xf = (xi-min(xi))*res(1) + pos(1);
yf = (yi-min(yi))*res(2) + pos(2);
zf = (zi-min(zi))*res(3) + pos(3);

