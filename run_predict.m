clear; clc; close all;

%target_dir = 'D:\Workspace\git\matlab\normal modes\at_helper\data\220503_range3-5km';
%target_dir = 'D:\Workspace\git\matlab\normal modes\at_helper\data\220503_range5-7km';
%target_dir = 'D:\Workspace\git\matlab\normal modes\at_helper\data\220503_range7-9km';
%target_dir = 'D:\Workspace\git\matlab\normal modes\at_helper\data\220503_range9-11km';
%target_dir = 'D:\Workspace\git\matlab\normal modes\at_helper\data\220503_range11-13km';


is_binary = true;
scaling = 1;

val_dir = fullfile(target_dir, '\val\');
pred_dir = fullfile(val_dir, 'val_pred');
indices = sort(readmatrix([val_dir 'indices.csv']));

overlayed_dir = fullfile(target_dir, '\overlay_mode');
val_overalyed_dir = fullfile(val_dir, '\overlay');
final_dir =  fullfile(val_dir, '\final');
lists = dir(fullfile(overlayed_dir, '*.jpg'));

val_indices = zeros(length(lists), 1);
for idx = 1:length(lists)
    name = lists(idx).name;
    for jdx = 1:length(indices)
        val_name = [num2str(indices(jdx)) '.jpg'];
        if strcmp(val_name, name)
            val_indices(idx) = 1;
            break;
        end
        
    end
end
val_indices = find(val_indices);
lists = lists(val_indices);
for idx = 1:length(val_indices)

    file = fullfile(lists(idx).folder, lists(idx).name);
    dst_path = [val_overalyed_dir '\' lists(idx).name];
    copyfile(file, dst_path);
    
    ov_data = imresize(imread(file), scaling);
    imshow(ov_data, 'InitialMagnification', 'fit')

%     pred_name = fullfile(pred_dir, [lists(idx).name(1:end-4) '.bmp']);
%     pred_data = imread(pred_name);
%     hold on;
%     h = imshow(pred_data, 'InitialMagnification', 'fit');
%     hold off;
%     set(h, 'AlphaData', pred_data);


end

