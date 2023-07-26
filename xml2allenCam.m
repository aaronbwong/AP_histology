function [cam_vector, cam_target] = xml2allenCam(xml_fn)
% adapted from codes of chongtianyifa
% see: https://github.com/petersaj/AP_histology/issues/14
%converts template slice in the xml file output from deepslice to allen
%coordinate and then find the vector of camera view and camtarget points
% check https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0216796
xml2allen_M = [0,0,25,0;-25,0,0,0;0,-25,0,0;13175,7975,0,1]; % voxel atlas (RAS axis orientation) to allen coordinate (PIR axis orientation)
xml_data=readstruct(xml_fn);
nslice = length(xml_data.slice);
anchoring = nan([nslice,9]);
for ii = 1:nslice
    anchoring(ii,:) = sscanf(xml_data.slice(ii).anchoringAttribute,'ox=%f&oy=%f&oz=%f&ux=%f&uy=%f&uz=%f&vx=%f&vy=%f&vz=%f');
end
slice_edge=zeros(size(anchoring,1),12);
cam_vector=zeros(size(anchoring,1),3); % each row for a camvector
cam_target=zeros(size(anchoring,1),3); % each row for a cametarget
for id_slice = 1:size(anchoring,1)
    o = anchoring(id_slice,1:3);
    u = anchoring(id_slice,4:6);
    v = anchoring(id_slice,7:9);
    xml_M= [u;v;o]; % convert pixels in template slice to voxel atlas coordinate system
    topLeft = o; % the most top left of template slice
    topRight = [1,0,1]*xml_M;
    bottomLeft = [0,1,1]*xml_M;
    bottomRight = [1,1,1]*xml_M;
    
    topLeft_a = [topLeft,1]*xml2allen_M/10; % coordinate in allen space, 10 um/voxel
    topRight_a = [topRight,1]*xml2allen_M/10;
    bottomLeft_a = [bottomLeft,1]*xml2allen_M/10;
    bottomRight_a = [bottomRight,1]*xml2allen_M/10;
    
    slice_edge(id_slice,:)=[topLeft_a(1:3),topRight_a(1:3),bottomLeft_a(1:3),bottomRight_a(1:3)];
    
    cam_vector(id_slice,:)=cross(topRight_a(1:3)-topLeft_a(1:3),bottomLeft_a(1:3)-topLeft_a(1:3)); % view (norm) of the slice
    cam_target(id_slice,:)=(topLeft_a(1:3)+bottomRight_a(1:3))/2; % or (topRight_a(1:3)+bottomLeft_a(1:3))/2 to find the center of slice
end

