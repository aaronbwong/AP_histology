function AP_grab_histology_ccf(tv,av,st,slice_im_path,slice_plane_fn)
% Grab CCF slices corresponding to histology slices
% loading of slice_plane_fn from deepslice is not necessary
% Andy Peters (peters.andrew.j@gmail.com)

% Initialize guidata
gui_data = struct;
gui_data.tv = tv;
gui_data.av = av;
gui_data.st = st;


% Load in slice images
gui_data.slice_im_path = slice_im_path;
slice_im_dir = dir([slice_im_path filesep '*.tif']);
slice_im_dir = [slice_im_dir;dir([slice_im_path filesep '*.jpg'])];
slice_im_dir = [slice_im_dir;dir([slice_im_path filesep '*.png'])];
slice_im_fn = natsortfiles(cellfun(@(path,fn) [path filesep fn], ...
    {slice_im_dir.folder},{slice_im_dir.name},'uni',false));
gui_data.slice_im = cell(length(slice_im_fn),1);
for curr_slice = 1:length(slice_im_fn)
    gui_data.slice_im{curr_slice} = imread(slice_im_fn{curr_slice});
end

% Create figure, set button functions
gui_fig = figure( ...
    'WindowScrollWheelFcn',@scroll_atlas_slice, ...
    'KeyPressFcn',@keypress);

% Set up axis for histology image
gui_data.histology_ax = subplot(1,2,1,'YDir','reverse'); 
hold on; axis image off;
gui_data.histology_im_h = image(gui_data.slice_im{1},'Parent',gui_data.histology_ax);
gui_data.curr_histology_slice = 1;
title(gui_data.histology_ax,'No saved atlas position');

% Set up 3D atlas axis
gui_data.atlas_ax = subplot(1,2,2, ...
    'ZDir','reverse','color','k', ...
    'XTick',[1,size(av,1)],'XTickLabel',{'Front','Back'}, ...
    'YTick',[1,size(av,3)],'YTickLabel',{'Left','Right'}, ...
    'ZTick',[1,size(av,2)],'ZTickLabel',{'Top','Bottom'});
hold on
axis vis3d equal manual
view([90,0]);
[ap_max,dv_max,ml_max] = size(tv);
xlim([1,ap_max]);
ylim([1,ml_max]);
zlim([1,dv_max]);
colormap(gui_data.atlas_ax,'gray');
caxis([0,400]);
gui_data.atlas_title = title(sprintf('Slice position: %d',0));

% Create CCF colormap
% (copied from cortex-lab/allenCCF/setup_utils
ccf_color_hex = st.color_hex_triplet;
ccf_color_hex(cellfun(@numel,ccf_color_hex)==5) = {'019399'}; % special case where leading zero was evidently dropped
ccf_cmap_c1 = cellfun(@(x)hex2dec(x(1:2)), ccf_color_hex, 'uni', false);
ccf_cmap_c2 = cellfun(@(x)hex2dec(x(3:4)), ccf_color_hex, 'uni', false);
ccf_cmap_c3 = cellfun(@(x)hex2dec(x(5:6)), ccf_color_hex, 'uni', false);
gui_data.ccf_cmap = ...
    horzcat(vertcat(ccf_cmap_c1{:}),vertcat(ccf_cmap_c2{:}),vertcat(ccf_cmap_c3{:}))./255;

% Set mode for atlas view (can be either TV, AV, or TV-AV)
gui_data.atlas_mode = 'TV';

% Create slice object and first slice point
gui_data.atlas_slice_plot = surface(gui_data.atlas_ax,'EdgeColor','none'); % Slice on 3D atlas
gui_data.atlas_slice_point = camtarget;

% Set up atlas parameters to save for histology
gui_data.slice_vector = nan(1,3);
gui_data.slice_points = nan(length(gui_data.slice_im),3);

% load the slice plance
if nargin>4
%     [gui_data.Cam_Vector, gui_data.Cam_Target] = csv2allenCam(slice_plane_fn); % to initinize the slice plane found by deepslice
    [gui_data.Cam_Vector, gui_data.Cam_Target] = xml2allenCam(slice_plane_fn); % to initinize the slice plane found by deepslice
    disp('make sure the order of slice planes corresponds to the histology images!!')
    for curr_slice = 1:length(slice_im_fn)
        [~,f_n,ext]=fileparts(slice_im_fn{curr_slice});
        disp([f_n,ext])
    end
    view(gui_data.Cam_Vector(1,:))
    gui_data.atlas_slice_point = gui_data.Cam_Target(1,:);
end

% Upload gui data
guidata(gui_fig,gui_data);

% Draw the first slice
update_atlas_slice(gui_fig);

% Print controls
CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'non-modal';
msgbox( ...
    {'\fontsize{12}' ...
    '\bf Controls: \rm' ...
    '1,2 : move histology slice' ...
    'm : change atlas display mode (TV/AV/TV-AV overlay)' ...
    'Arrow keys: rotate CCF atlas', ...
    'Scroll wheel: move CCF slice in/out of plane', ...
    'Enter: set current histology and CCF slice pair', ...
    'Escape: save and close'}, ...
    'Controls',CreateStruct);

end 

function keypress(gui_fig,eventdata)

% Get guidata
gui_data = guidata(gui_fig);

switch eventdata.Key
    
    % Arrow keys: rotate atlas slice
    case 'leftarrow'
        set(gui_data.atlas_ax,'View',get(gui_data.atlas_ax,'View') + [1,0]);
        update_atlas_slice(gui_fig)
    case 'rightarrow'
        set(gui_data.atlas_ax,'View',get(gui_data.atlas_ax,'View') + [-1,0]);
        update_atlas_slice(gui_fig)
    case 'uparrow'
        set(gui_data.atlas_ax,'View',get(gui_data.atlas_ax,'View') + [0,-1]);
        update_atlas_slice(gui_fig)
    case 'downarrow'
        set(gui_data.atlas_ax,'View',get(gui_data.atlas_ax,'View') + [0,1]);
        update_atlas_slice(gui_fig)
    
    % 1/2 keys: cycle through histology slices
    % (if there's a saved plane point, move atlas to that position)
    case '1'
        gui_data.curr_histology_slice = max(gui_data.curr_histology_slice - 1,1);            
        guidata(gui_fig,gui_data);
        update_histology_slice(gui_fig);
        
    case '2'
        gui_data.curr_histology_slice = ...
            min(gui_data.curr_histology_slice + 1,length(gui_data.slice_im));
        guidata(gui_fig,gui_data);
        update_histology_slice(gui_fig);

    % M key: switch atlas display mode
    case 'm'
        atlas_slice_modes = {'TV','AV','TV-AV'};
        curr_atlas_mode_idx = strcmp(gui_data.atlas_mode,atlas_slice_modes);
        gui_data.atlas_mode = atlas_slice_modes{circshift(curr_atlas_mode_idx,[0,1])};
        guidata(gui_fig,gui_data);
        update_atlas_slice(gui_fig);

    % Enter: save slice coordinates
    case 'return'        
        % Store camera vector and point
        % (Note: only one camera vector used for all slices, overwrites)
        gui_data.slice_vector = get_camera_vector(gui_data);
        gui_data.slice_points(gui_data.curr_histology_slice,:) = ...
            gui_data.atlas_slice_point;
        guidata(gui_fig,gui_data);
                
        update_histology_slice(gui_fig);
        title(gui_data.histology_ax,'New saved atlas position');
        
    % Escape: save and exit
    case 'escape'
        opts.Default = 'Yes';
        opts.Interpreter = 'tex';
        user_confirm = questdlg('\fontsize{15} Save and quit?','Confirm exit',opts);
        if strcmp(user_confirm,'Yes')
            
            % Check that a CCF slice point exists for each histology slice
            if any(isnan(gui_data.slice_points(:)))
                createmode = struct;
                createmode.Interpreter = 'tex';
                createmode.WindowStyle = 'modal';
                msgbox('\fontsize{12} Some histology slice(s) not assigned CCF slice', ...
                    'Not saving','error',createmode);
                return
            end
            
            % Go through each slice, pull full-resolution atlas slice and
            % corrsponding coordinates       
            histology_ccf_init = cell(length(gui_data.slice_im),1);
            histology_ccf = struct( ...
                'tv_slices',histology_ccf_init, ...
                'av_slices',histology_ccf_init, ...
                'plane_ap',histology_ccf_init, ...
                'plane_ml',histology_ccf_init, ...
                'plane_dv',histology_ccf_init);
            
            h = waitbar(0,'Saving atlas slices...');
            for curr_slice = 1:length(gui_data.slice_im)
                gui_data.atlas_slice_point = gui_data.slice_points(curr_slice,:);
                [histology_ccf(curr_slice).tv_slices, ...
                    histology_ccf(curr_slice).av_slices, ...
                    histology_ccf(curr_slice).plane_ap, ...
                    histology_ccf(curr_slice).plane_ml, ...
                    histology_ccf(curr_slice).plane_dv] = ...
                    grab_atlas_slice(gui_data,1);
                waitbar(curr_slice/length(gui_data.slice_im),h, ...
                    ['Saving atlas slices (' num2str(curr_slice) '/' num2str(length(gui_data.slice_im)) ')...']);
            end                     
            close(h);
            slice_points = gui_data.slice_points;
            slice_vector = gui_data.slice_vector;
            save_fn = [gui_data.slice_im_path filesep 'histology_ccf.mat'];
            save(save_fn,'histology_ccf','slice_points','slice_vector','-v7.3');
            close(gui_fig);            
        end
end

end

function update_histology_slice(gui_fig)
% Draw histology slice (and move atlas if saved position)

% Get guidata
gui_data = guidata(gui_fig);

% Set next histology slice
set(gui_data.histology_im_h,'CData',gui_data.slice_im{gui_data.curr_histology_slice})

% If there's a saved atlas position, move atlas to there
if all(~isnan(gui_data.slice_points(gui_data.curr_histology_slice,:)))
    gui_data.atlas_slice_point = ...
        gui_data.slice_points(gui_data.curr_histology_slice,:);
    title(gui_data.histology_ax,'Saved atlas position')
else
    title(gui_data.histology_ax,'No saved atlas position')
    if isfield(gui_data,'Cam_Vector')
        view(gui_data.Cam_Vector(gui_data.curr_histology_slice,:)) % If there's no saved atlas position, loading atlas to plane found by deepslice
        gui_data.atlas_slice_point = gui_data.Cam_Target(gui_data.curr_histology_slice,:);
    end
end

% Upload gui data
guidata(gui_fig,gui_data);
update_atlas_slice(gui_fig);

end

function cam_vector = get_camera_vector(gui_data)
% Get the camera viewing vector to define atlas slice plane

% Grab current camera angle

% (Old way: more confusing, easily messed up by axes directions)
% [cam_az,cam_el] = view(gui_data.atlas_ax);
% 
% % Camera azimuth is 90 degrees offset from spherical standard (?!)
% cam_az_sphere = cam_az - 90;
% % Camera elevation is reversed (because of CCF orientation)
% cam_el_sphere = -cam_el;
% 
% [cam_vector_x,cam_vector_y,cam_vector_z] = ...
%     sph2cart(deg2rad(cam_az_sphere),deg2rad(cam_el_sphere),1);
% cam_vector = [cam_vector_x,cam_vector_y,cam_vector_z];

% (New way: just a normalized line from the camera to the center)
curr_campos = campos(gui_data.atlas_ax);
curr_camtarget = camtarget(gui_data.atlas_ax);
cam_vector = (curr_camtarget - curr_campos)./norm(curr_camtarget - curr_campos);

end

function scroll_atlas_slice(gui_fig,eventdata)
% Move point to draw atlas slice perpendicular to the camera

% Get guidata
gui_data = guidata(gui_fig);

% Move slice point along camera -> center axis
cam_vector = get_camera_vector(gui_data);

% Move slice point
gui_data.atlas_slice_point = gui_data.atlas_slice_point + ...
    eventdata.VerticalScrollCount*cam_vector;  % A positive or negative number that indicates the direction and number of scroll wheel clicks. The vertical scroll count is the sum of all scroll wheel clicks that occurred since the last time the callback executed.

% Upload gui data
guidata(gui_fig, gui_data);

% Update slice
update_atlas_slice(gui_fig)

end

function update_atlas_slice(gui_fig)
% Draw atlas slice through plane perpendicular to camera through set point

% Get guidata
gui_data = guidata(gui_fig);

% Get slice (larger spacing for faster pulling)
[tv_slice,av_slice,plane_ap,plane_ml,plane_dv] = grab_atlas_slice(gui_data,3);

% Update the slice display (depending on display mode)
switch gui_data.atlas_mode
    case 'TV'
        atlas_slice = tv_slice;
        colormap(gray);caxis([0,516]);
    case 'AV'
        av_boundaries = round(conv2(av_slice,ones(2)./4,'same')) ~= av_slice;
        atlas_slice = imoverlay(mat2gray(tv_slice,[0,516]),av_boundaries,'r');
        caxis([0,1]);
    case 'TV-AV'
        atlas_slice = av_slice;
        colormap(gui_data.ccf_cmap)
        caxis([1,size(gui_data.ccf_cmap,1)])
end
set(gui_data.atlas_slice_plot,'XData',plane_ap,'YData',plane_ml,'ZData',plane_dv,'CData',atlas_slice);

% Upload gui_data
guidata(gui_fig, gui_data);

end

function [tv_slice,av_slice,plane_ap,plane_ml,plane_dv] = grab_atlas_slice(gui_data,slice_px_space)
% Grab anatomical and labelled atlas within slice

% Get plane normal to the camera -> center axis, grab voxels on plane
cam_vector = get_camera_vector(gui_data);
plane_offset = -(cam_vector*gui_data.atlas_slice_point');

% Define a plane of points to index
% (the plane grid is defined based on the which cardinal plan is most
% orthogonal to the plotted plane. this is janky but it works)

[~,cam_plane] = max(abs(cam_vector./norm(cam_vector)));

switch cam_plane
    
    % Note: ML and DV directions are flipped to match 2D histology and 3D
    % atlas axes, so make ML and DV coordinates go backwards for true CCF
    % coordinates
    
    case 1
        [plane_ml,plane_dv] = ...
            meshgrid(1:slice_px_space:size(gui_data.tv,3), ...
            1:slice_px_space:size(gui_data.tv,2));
        plane_ap = ...
            (cam_vector(2)*plane_ml+cam_vector(3)*plane_dv + plane_offset)/ ...
            -cam_vector(1);
        
    case 2
        [plane_ap,plane_dv] = ...
            meshgrid(1:slice_px_space:size(gui_data.tv,1), ...
            1:slice_px_space:size(gui_data.tv,2));
        plane_ml = ...
            (cam_vector(1)*plane_ap+cam_vector(3)*plane_dv + plane_offset)/ ...
            -cam_vector(2);
        
    case 3
        [plane_ap,plane_ml] = ...
            meshgrid(size(gui_data.tv,1):-slice_px_space:1, ...
            1:slice_px_space:size(gui_data.tv,3));
        plane_dv = ...
            (cam_vector(1)*plane_ap+cam_vector(2)*plane_ml + plane_offset)/ ...
            -cam_vector(3);
        
end

% Get the coordiates on the plane
ap_idx = round(plane_ap);
ml_idx = round(plane_ml);
dv_idx = round(plane_dv);

% Find plane coordinates in bounds with the volume
% (CCF coordinates: [AP,DV,ML])
use_ap = ap_idx > 0 & ap_idx < size(gui_data.tv,1);
use_dv = dv_idx > 0 & dv_idx < size(gui_data.tv,2);
use_ml = ml_idx > 0 & ml_idx < size(gui_data.tv,3);
use_idx = use_ap & use_ml & use_dv;

curr_slice_idx = sub2ind(size(gui_data.tv),ap_idx(use_idx),dv_idx(use_idx),ml_idx(use_idx));

% Find plane coordinates that contain brain
curr_slice_isbrain = false(size(use_idx));
curr_slice_isbrain(use_idx) = gui_data.av(curr_slice_idx) > 0;

% Index coordinates in bounds + with brain
grab_pix_idx = sub2ind(size(gui_data.tv),ap_idx(curr_slice_isbrain),dv_idx(curr_slice_isbrain),ml_idx(curr_slice_isbrain));

% Grab pixels from (selected) volume
tv_slice = nan(size(use_idx));
tv_slice(curr_slice_isbrain) = gui_data.tv(grab_pix_idx);

av_slice = nan(size(use_idx));
av_slice(curr_slice_isbrain) = gui_data.av(grab_pix_idx);

% Update slice position title
plane_offset_mm = plane_offset/100; % CCF = 10um voxels
set(gui_data.atlas_title,'string', ...
    sprintf('Slice position: %.2f mm',plane_offset_mm));

end













