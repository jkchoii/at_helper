clear; clc; close all;

%target_dir = 'D:\Workspace\git\matlab\normal modes\at_helper\data\220502_range3-5km';
%target_dir = 'D:\Workspace\git\matlab\normal modes\at_helper\data\220502_range5-7km';
%target_dir = 'D:\Workspace\git\matlab\normal modes\at_helper\data\220502_range7-9km';
%target_dir = 'D:\Workspace\git\matlab\normal modes\at_helper\data\220502_range9-11km';

root_dir = 'D:\Workspace\git\matlab\normal modes\at_helper\data\';

%current_date = datestr(now, 'yymmdd-HHMMSS');
current_date = '';
current_dir = fullfile(root_dir, ['test_result/' current_date]);
mkdir(current_dir);
test_dir = fullfile(current_dir, 'train');

train_cand = ["3-5", "5-7", "7-9", "9-11", "11-13"];

target_size = [256,256];

fig = figure('Renderer', 'painters', 'Position', [-900 474 572 572]);


for loop = 1:length(train_cand)
    test_dir_target = [test_dir char(train_cand(loop)) 'km\'];
    mkdir(test_dir_target)
    for idx = 1:length(train_cand)
        target_dir = [root_dir [datestr(now, 'yymmdd') '_range'] char(train_cand(idx)) 'km'];
        pred_dir = [test_dir_target char(train_cand(idx))];
        mkdir(pred_dir);

        overlay_dir = fullfile(target_dir, 'val/overlay');
        pred_imag_lists = dir([pred_dir '\*.bmp']);
        overlay_image_lists = dir([overlay_dir '/*.jpg']);
        for jdx = 1:length(pred_imag_lists)
            pred_img = imread(fullfile(pred_imag_lists(jdx).folder, pred_imag_lists(jdx).name));
            pred_img = cat(3, pred_img, pred_img,pred_img);
            pred_img = imresize(pred_img, target_size);
            tmp = pred_imag_lists(jdx).name;
            raw_name = [tmp(1:end-4) 'comp' '.jpg'];
            name = [tmp(1:end-4) '.jpg'];
            overlay_img = imread(fullfile(overlay_image_lists(jdx).folder, name));
            overlay_img = imresize(overlay_img, target_size);
            
            img = [pred_img overlay_img];
            
            concat_dir = fullfile(pred_imag_lists(jdx).folder, 'concat');
            mkdir(concat_dir);
            imwrite(img, fullfile(concat_dir, raw_name));
        end

    end
    


end
