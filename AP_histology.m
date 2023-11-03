function AP_histology
% Toolbar GUI for running histology pipeline

% Set up the gui
screen_size_px = get(0,'screensize');
gui_aspect_ratio = 2; % width/length
gui_width_fraction = 0.3; % fraction of screen width to occupy
gui_border = 50; % border from gui to screen edge
gui_width_px = screen_size_px(3).*gui_width_fraction;
gui_height_px = gui_width_px/gui_aspect_ratio;
gui_position = [...
    gui_border, ... % left x
    screen_size_px(4)-(gui_height_px+gui_border+50), ... % bottom y
    gui_width_px,gui_height_px]; % width, height

histology_toolbar_gui = figure('Toolbar','none','Menubar','none','color','w', ...
    'Name','AP Histology', ...
    'Units','pixels','Position',gui_position);

% Set up the text to display coordinates
gui_data.gui_text = annotation('textbox','String','','interpreter','tex', ...
    'Units','normalized','Position',[0,0,1,1],'VerticalAlignment','top', ...
    'FontSize',12,'FontName','Consolas','PickableParts','none');

% File menu
gui_data.menu.file = uimenu(histology_toolbar_gui,'Text','File selection');
uimenu(gui_data.menu.file,'Text','Set raw image path','MenuSelectedFcn',{@set_image_path,histology_toolbar_gui});
uimenu(gui_data.menu.file,'Text','Set processing save path','MenuSelectedFcn',{@set_save_path,histology_toolbar_gui});

% Preprocessing menu
gui_data.menu.preprocess = uimenu(histology_toolbar_gui,'Text','Image preprocessing');
uimenu(gui_data.menu.preprocess,'Text','Create slice images','MenuSelectedFcn', ...
    {@ap_histology.create_slice_images,histology_toolbar_gui});
uimenu(gui_data.menu.preprocess,'Text','Rotate & center slices','MenuSelectedFcn', ...
    {@ap_histology.rotate_center_slices,histology_toolbar_gui});
uimenu(gui_data.menu.preprocess,'Text','Flip & re-order slices','MenuSelectedFcn', ...
    {@ap_histology.flip_reorder_slices,histology_toolbar_gui});

% Atlas menu
gui_data.menu.atlas = uimenu(histology_toolbar_gui,'Text','Atlas alignment');
uimenu(gui_data.menu.atlas,'Text','Choose histology atlas slices','MenuSelectedFcn', ...
    {@ap_histology.match_histology_atlas,histology_toolbar_gui});
uimenu(gui_data.menu.atlas,'Text','Auto-align histology/atlas slices','MenuSelectedFcn', ...
    {@ap_histology.align_auto_histology_atlas,histology_toolbar_gui});
uimenu(gui_data.menu.atlas,'Text','Manual align histology/atlas slices','MenuSelectedFcn', ...
    {@ap_histology.align_manual_histology_atlas,histology_toolbar_gui});

% Annotation menu
gui_data.menu.annotation = uimenu(histology_toolbar_gui,'Text','Annotation');
uimenu(gui_data.menu.annotation,'Text','Neuropixels probes','MenuSelectedFcn', ...
    {@ap_histology.annotate_neuropixels,histology_toolbar_gui});

% View menu
gui_data.menu.view = uimenu(histology_toolbar_gui,'Text','View');
uimenu(gui_data.menu.view,'Text','View aligned histology','MenuSelectedFcn', ...
    {@ap_histology.view_aligned_histology,histology_toolbar_gui});

% Create GUI variables
gui_data.image_path = char;
gui_data.save_path = char;

% Store guidata
guidata(histology_toolbar_gui,gui_data);

% Update GUI text
update_gui(histology_toolbar_gui);

end

function set_image_path(h,eventdata,histology_toolbar_gui)

% Get guidata
gui_data = guidata(histology_toolbar_gui);

% Pick image path
gui_data.image_path = uigetdir([],'Select path with raw images');

% Store guidata
guidata(histology_toolbar_gui,gui_data);

% Update GUI text
update_gui(histology_toolbar_gui);

end

function set_save_path(h,eventdata,histology_toolbar_gui)

% Get guidata
gui_data = guidata(histology_toolbar_gui);

% Pick image path
gui_data.save_path = uigetdir([],'Select path to save processing');

% Store guidata
guidata(histology_toolbar_gui,gui_data);

% Update GUI text
update_gui(histology_toolbar_gui);

end


function update_gui(histology_toolbar_gui)

% Get guidata
gui_data = guidata(histology_toolbar_gui);

% Check for files
raw_images_present = ~isempty(gui_data.image_path) && ...
    ~isempty(dir(fullfile(gui_data.image_path,'*.tif')));
processed_images_present = ~isempty(gui_data.save_path) && ...
    ~isempty(dir(fullfile(gui_data.save_path,'*.tif')));
atlas_slices_present = ~isempty(gui_data.save_path) && ...
    exist(fullfile(gui_data.save_path,'histology_ccf.mat'),'file');
atlas_alignment_present = ~isempty(gui_data.save_path) && ...
    exist(fullfile(gui_data.save_path,'atlas2histology_tform.mat'),'file');
neuropixels_annotations_present = ~isempty(gui_data.save_path) && ...
    exist(fullfile(gui_data.save_path,'probe_ccf.mat'),'file');

% Enable/disable appropriate menu options
gui_data.menu.preprocess.Enable = raw_images_present;
[gui_data.menu.preprocess.Children(2:end).Enable] = deal(processed_images_present);

gui_data.menu.atlas.Enable = processed_images_present;
[gui_data.menu.atlas.Children(2:end).Enable] = deal(processed_images_present & atlas_slices_present);

gui_data.menu.annotation.Enable = atlas_alignment_present;
gui_data.menu.view.Enable = atlas_alignment_present;


% Set text
image_path_text = sprintf('Raw image path:       %s',strrep(gui_data.image_path,filesep,repmat(filesep,1,2)));
save_path_text = sprintf('Processing save path: %s',strrep(gui_data.save_path,filesep,repmat(filesep,1,2)));

% Check for present files
% (images)
if raw_images_present
    n_raw_images = length(dir(fullfile(gui_data.image_path,'*.tif')));
else
    n_raw_images = 0;
end
n_raw_images_text = sprintf('Raw images: %d',n_raw_images);

if processed_images_present
    n_processed_images = length(dir(fullfile(gui_data.save_path,'*.tif')));
else
    n_processed_images = 0;
end
n_processed_images_text = sprintf('Processed images: %d',n_processed_images);

% (alignment)
if atlas_slices_present
    histology_atlas_text = 'Histology atlas slices: YES';
else
    histology_atlas_text = 'Histology atlas slices: NO';
end

if atlas_alignment_present
    alignment_text = 'Histology-atlas alignment: YES';
else
    alignment_text = 'Histology-atlas alignment: NO';
end

% (annotations)
if neuropixels_annotations_present
    annotations_text = 'Neuropixels probe annotations: YES';
else
    annotations_text = '';
end
gui_text = { ...
    '\bf --File paths \rm', ...
    image_path_text,save_path_text,'\newline', ...
    '\bf --Images \rm', ...
    n_raw_images_text,n_processed_images_text, '\newline', ...
    '\bf --Atlas alignment \rm', ...
    histology_atlas_text,alignment_text, '\newline', ...
    '\bf --Annotations \rm', ...
    annotations_text};

set(gui_data.gui_text,'String',gui_text(cellfun(@(x) ~isempty(x),gui_text)))

end



















