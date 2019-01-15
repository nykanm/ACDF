%% Visualise the Changes in the structure from beginning to end
% Load the procedure beginning file
V_start = voxelread;
% Load the procedure ending file
V_final = voxelread;

% Calculate voxels that were removed
V_diff = V_final - V_start;
% Visualise the voxels that were removed
%volumeViewer(V_diff);

%% Visualise the structures that were removed
figure;
res = [0.02 0.02 0.02];
pos = [3.42573 -1.56917 -8.35422];
[xf, yf, zf, c] = vol2vec(V_start,res,pos);
scatter3(-xf,-yf,zf,20,c,'s','filled');
colormap(jet);
axis equal;
hold on;
[xf, yf, zf, c] = vol2vec(V_final,res,pos);
scatter3(xf,yf,zf,20,'g','s','filled');
%colormap(bone);
legend('Removed','Not removed');

%% Track the movement of a tool relative to the structure
% Load csv file
% Only use virtual position data
hold on;
% Restrict to only positions when tip in contact with disk annulus
contact = find(x.ContactVoxelsC4C5DiscNucleus>0);
coorInContact = [x.VirtualToolTipPosition_X(contact) x.VirtualToolTipPosition_Y(contact) x.VirtualToolTipPosition_Z(contact)];
scatter3(coorInContact(:,1),coorInContact(:,2),coorInContact(:,3));
