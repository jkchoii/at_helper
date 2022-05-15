classdef savePackage
    %SAVEPACKAGE 이 클래스의 요약 설명 위치
    %   자세한 설명 위치
    
    properties
        data_path
        cmap
        data_idx
        
    end
    
    methods
        function obj = savePackage(path)
            obj.data_path = fullfile(pwd, path);
            obj.cmap = flipud(colormap('gray'));
            close();
        end

        function [g_cells] = assignValidGroupTimeIndexes(obj, groups, target_freq)
            max_size = size(groups(1).vg, 1);
            num_ranges = size(groups(1).t_group, 3);
            num_modes = length(groups);
            g_cells = cell(length(groups), 1);
            

            for idx = 1:num_modes % num modes
                [~, max_freq_idx] = min(abs(groups(idx).f_group - max(target_freq)));
                g_cells{idx} = repmat(struct('start', 0, ...
                                             'end', 0, ...
                                             'valid_idx', zeros(max_size, 1)), ...
                                             [num_ranges, 1]);
                
                for jdx = 1:num_ranges % num ranges
                    group_t = groups(idx).t_group(:,:,jdx);
                    g_cells{idx}(jdx).start = group_t(max_freq_idx);
                    g_cells{idx}(jdx).end = max(group_t);
                    g_cells{idx}(jdx).valid_idx = group_t >= g_cells{idx}(jdx).start & group_t <= g_cells{idx}(jdx).end;
                    %groups(idx).t_group(:,:,jdx) = groups(idx).t_group(t_valid_idx,:,jdx);
                    %groups(idx).f_group = groups(idx).f_group(t_valid_idx);
                    %groups(idx).vg = groups(idx).vg(t_valid_idx);
                end

            end

        end

        function [g_grids] = get2DGrid(obj, spectro, g_cells, groups, target_freq)
            
            max_size = max([length(spectro(1).f), length(spectro(1).t)]);
            num_ranges = length(spectro);
            spectro_cells = cell(num_ranges, 1);

            for idx = 1:num_ranges
                spectro_cells{idx} = struct('f_idx', zeros(max_size, 1), ...
                                            't_idx', zeros(max_size, 1), ...
                                            't_start_point', 0, ...
                                            't_len', 0);

                t_start = realmax;
                t_end = -realmin;
                for jdx = 1:length(g_cells)
                    if g_cells{jdx}(idx).start < t_start
                        t_start = g_cells{jdx}(idx).start;
                    end
                    if g_cells{jdx}(idx).end > t_end
                        t_end = g_cells{jdx}(idx).end;
                    end
                end

                spectro_cells{idx}.t_idx = spectro(idx).t >= t_start & spectro(idx).t <= t_end;
                [~, spectro_cells{idx}.t_start_point] = min(abs(spectro(idx).t - t_start));
                %spectro_cells{idx}.t_idx = spectro_cells{idx}.t_idx | t_idx;

                spectro_cells{idx}.t_len = length(find(spectro_cells{idx}.t_idx));
                %f_idx = spectro(idx).f >= min(target_freq) & spectro(idx).f <= max(target_freq);
                f_idx = spectro(idx).f >= 0 & spectro(idx).f <= max(target_freq);
                spectro_cells{idx}.f_idx = f_idx;

            end

            max_cols = 0;
            max_range_idx = 0;
            for idx = 1:length(spectro_cells) % num_range
                if max_cols < spectro_cells{idx}.t_len
                    max_cols = spectro_cells{idx}.t_len;
                    max_range_idx = idx;
                end
            end
            
            t_vec = spectro(max_range_idx).t(spectro_cells{max_range_idx}.t_idx);
            f_vec = spectro(max_range_idx).f(spectro_cells{max_range_idx}.f_idx);
            g_grids = repmat(struct('grid', zeros(length(f_vec), length(t_vec)), 't', t_vec, 'f', f_vec), [num_ranges, 1]);
            for idx = 1:num_ranges
                slice_start = spectro_cells{idx}.t_start_point;
                slice_end = (spectro_cells{idx}.t_start_point+max_cols-1);
                if slice_end > size(spectro(idx).t, 1)
                    slice_end = size(spectro(idx).t, 1);
                    slice_start = slice_end - max_cols + 1;
                end

                g_grids(idx).t = spectro(idx).t(slice_start:slice_end);
                g_grids(idx).f = spectro(idx).f(spectro_cells{idx}.f_idx);

            end

        end

        function [imag] = get2DGroupVelocityImages(obj, g_grids, g_cells, groups)
            num_ranges = length(g_grids);
            num_modes = length(groups);
            imag = zeros(size(g_grids(1).grid, 1), size(g_grids(1).grid, 2), num_ranges);
            for idx = 1:num_modes
                for jdx = 1:num_ranges
                    % ground truth
                    t = groups(idx).t_group(:,:,jdx);
                    
                    % downsampled grids
                    [~, time_start_idx] = min(abs(g_cells{idx}(jdx).start - g_grids(jdx).t));
                    [~, time_end_idx] = min(abs(g_cells{idx}(jdx).end - g_grids(jdx).t));

                    [~, freq_start_idx] = min(abs(min(groups(idx).f_group) - g_grids(jdx).f));
                    [~, freq_end_idx] = min(abs(max(groups(idx).f_group) - g_grids(jdx).f));
                    
                    

                    imag(freq_start_idx:freq_end_idx, time_start_idx:time_end_idx, jdx) = 1;
                    %+ imag(freq_start_idx:freq_end_idx, time_start_idx:time_end_idx);
                    

                end
            end

            

            
        end

        function [groups] = getTruncatedGroupVelocities(obj, spectro, groups, target_freq)
            
            g_cells = assignValidGroupTimeIndexes(obj, groups, target_freq);
            g_grids = get2DGrid(obj, spectro, g_cells, groups, target_freq);
            imag = get2DGroupVelocityImages(obj, g_grids, g_cells, groups);
            
        end

        function saveSpectrogram(obj, sig_t, spectro, groups, target_freq)
            
            groups = getTruncatedGroupVelocities(obj, spectro, groups, target_freq);

            for idx = 1:length(spectro)
                f_idx = (spectro(idx).f >= min(target_freq)) & spectro(idx).f <= max(target_freq);
                
            end
            
        end

        function [obj, target_path] = makeDirectories(obj, dir_name)
            head_file_name = datestr(now, 'yymmdd');
            target_path = fullfile(obj.data_path, [head_file_name dir_name]);
            [~, ~, matlab_msg] = mkdir(target_path);

            if ~strcmp('MATLAB:MKDIR:DirectoryExists', matlab_msg)
                mkdir(fullfile(target_path, 'real_img'));
                mkdir(fullfile(target_path, 'imag_img'));
                mkdir(fullfile(target_path, 'abs_img'));
                mkdir(fullfile(target_path, 'real_bin'));
                mkdir(fullfile(target_path, 'imag_bin'));
                mkdir(fullfile(target_path, 'abs_bin'));
                mkdir(fullfile(target_path, 'abs_bin_log'));


                mkdir(fullfile(target_path, 'test'));
                mkdir(fullfile(target_path, 'train'));
                mkdir(fullfile([target_path, '\train'], 'overlay'));
                mkdir(fullfile(target_path, 'val'));
                mkdir(fullfile([target_path, '\val'], 'val_pred'));
                mkdir(fullfile([target_path, '\val'], 'overlay'));
                
                mkdir(fullfile(target_path, 'label_img'));
                mkdir(fullfile(target_path, 'label_bin'));
                mkdir(fullfile(target_path, 'label_bin_log'));
                mkdir(fullfile(target_path, 'overlay_mode'));
                mkdir(fullfile(target_path, 'overlay_full'));
                mkdir(fullfile(target_path, 'checkpoints'));
                obj.data_idx = 0;

            end

        end

        function [obj] = writeDataInformation(obj, target_path, info_name, data_idx, range_idx, pos, mode_number, type, data_size)
            
            sd = pos.sd;
            rd = pos.rd;
            range = pos.range_src2recv;

            
            fid = fopen(fullfile(target_path, info_name), 'a');
            
            msg = sprintf("filename: %s, data_size: (%d, %d), sd: %4f, rd: %4f, range: %4f, mode_number: %d, %s", num2str(data_idx), data_size(1), data_size(2), sd, rd, range, mode_number, type );
            fprintf(fid, "%s\n", msg);
            fclose(fid);

        end
        
        function [obj] = SetStartDataIndex(obj, num)
            obj.data_idx = num;

        end

        function [obj] = saveSpectrogramFixedTime(obj, dir_name, sig, spectro, groups, stft_pack, roi_freq, pos, full_modes, target_size)
            
            info_name = datestr(now, 'yymmdd_HHMMSS');
            if ~isempty(dir_name)
                dir_name = ['_' char(dir_name)];
            end

            info_name = [info_name dir_name '_info.txt'];
            

            [obj, target_path] = makeDirectories(obj, dir_name);
            
            
            num_range = length(spectro);

            border_pad = 2;
            %figure('Renderer', 'painters', 'Position', [-900 474 572 + border_pad 572 + border_pad]);
            
            init = false;

            colors = [];
            
            total_iter = num_range * size(spectro(1).s, 1);
            count = 0;
            for idx = 1:num_range
                num_modes = size(spectro(idx).s, 3);
                for jdx = 1:num_modes
                    count = count + 1;
                    
                    if mod(fix(count / total_iter * 100), 2 )
                        proc_ratio =  count / total_iter * 100;
                        disp(proc_ratio);
                    end

                    [~, min_freq_idx] = min(abs(min(roi_freq) - spectro(idx).f));
                    [~, max_freq_idx] = min(abs(max(roi_freq) - spectro(idx).f));
                    t = spectro(idx).t;
                    f = stft_pack.fftfreq(min_freq_idx:max_freq_idx);

                    %s = flipud(spectro(idx).s(min_freq_idx:max_freq_idx,:,jdx));
                    s = spectro(idx).s(min_freq_idx:max_freq_idx,:,jdx);

                    if init == false
                        %target_size = [512, 512];
                        %figure('Renderer', 'painters', 'Position', [-900 474 size(s, 2) + border_pad size(s, 1) + border_pad]);
                        figure('Renderer', 'painters', 'Position', [-900 474 target_size(2) + border_pad target_size(1) + border_pad]);
                        h = gca;
                        imagesc(zeros(10));
                        hold on; plot([1 2 3]);
                        hold off;
                        h.Children(1).XData = [];
                        h.Children(1).YData = [];
                        h.Children(2).CData = [];
                        init = true;

                        figure('Renderer', 'painters', 'Position', [-600 474 target_size(2) + border_pad target_size(1) + border_pad]);
                        h2 = gca;
                        imagesc(zeros(10));
                        hold on; 
                        for plot_idx = 1:length(groups{1})
                            plot(plot_idx:(plot_idx+1));
                        end
                        hold off;
                        for plot_idx = 1:length(groups{1})
                            h2.Children(plot_idx).XData = [];
                            h2.Children(plot_idx).YData = [];
                        end
                        h2.Children(end).CData = [];
                        
                    end

                    % normalization [0 1]
                    [s_real, s_imag, s_abs] = getNormalizedSpectrogram(obj, s, t, f, target_size);
                    s_abs = s_abs.^2;


                    % standardization [-1 1]
                    %[s_real, s_imag, s_abs] = getStandardizedSpectrogram(obj, s);

                    %file_name = sprintf("%s\\real\\%d.bmp", target_path, obj.data_idx);

%                     if num_modes ~= 1
%                         if groups{idx}(jdx).valid_flag == true
%                             need_save = true;
%                         else 
%                             need_save = false;
%                         end
%                     else
%                         need_save = true;
%                     end
                    need_save = true;
                    if need_save == true
                        writeDataInformation(obj, target_path, info_name, obj.data_idx, idx, pos(idx), jdx, 'real', size(s_abs));
                        saveJpgImage(obj, h, t, f, s_real, target_path, 'real');
                        
                        %file_name = sprintf("%s\\imag\\%d.bmp", target_path, obj.data_idx);
                        writeDataInformation(obj, target_path, info_name, obj.data_idx, idx, pos(idx), jdx, 'imag', size(s_abs));
                        saveJpgImage(obj, h, t, f, s_imag, target_path, 'imag');
    
                        %file_name = sprintf("%s\\abs\\%d.bmp", target_path, obj.data_idx);
                        writeDataInformation(obj, target_path, info_name, obj.data_idx, idx, pos(idx), jdx, 'abs', size(s_abs));
                        saveJpgImage(obj, h, t, f, s_abs, target_path, 'abs');
                    


                        s_empty = zeros(size(s_abs));
    %                     colors.image = obj.cmap;
    %                     colors.plot = [0 0 0];
    %                     colors.plot_marker = 'none';
    %                     colors.plot_style = '-';
                        writeDataInformation(obj, target_path, info_name, obj.data_idx, idx, pos(idx), jdx, 'label', size(s_abs));
    
                        if full_modes == true
                            plot_num_modes = 1:length(groups{idx});
                        else
                            plot_num_modes = jdx;
                        end
                        %file_name = sprintf("%s\\label\\%d_mask.gif",target_path, obj.data_idx);
                        %num_classes = saveLabelImage(obj, h, t, f, s_empty, idx, plot_num_modes, groups, colors, file_name);
                        num_classes = saveLabel(obj, h, t, f, s_empty, idx, plot_num_modes, groups{idx}, colors, target_path, 'label');
    
                        
                        %%
                        colors.image = 'parula';
                        colors.plot = [1 0 0];
                        colors.plot_marker = 'none';
                        colors.plot_style = '-';
                        
                        file_name = sprintf("%s\\overlay_mode\\%d.jpg", target_path, obj.data_idx);
                        saveLabelFullImage(obj, h2, t, f, s_abs, idx, plot_num_modes, groups{idx}, colors, file_name);
                        %saveLabelImage(obj, h, t, f, s_abs, idx, plot_num_modes, groups, colors, file_name);
    
    
    
                        obj.data_idx = obj.data_idx + 1;

                    end

                end
                file_name = sprintf("%s\\overlay_full\\range%f data_name%d_to_%d.jpg",target_path, pos(idx).range_src2recv, (obj.data_idx-num_modes), obj.data_idx-1);
                s_full = sum(spectro(idx).s(min_freq_idx:max_freq_idx,:,:), 3);
                [~, ~, s_abs] = getNormalizedSpectrogram(obj, s_full, t, f, target_size);
                colors.image = 'parula';
                colors.plot = [1 0 0];
                colors.plot_marker = '*';
                colors.plot_style = 'none';
                
                %saveLabelImage(obj, h, t, f, s_abs, idx, 1:length(groups), groups, colors, file_name);
                saveLabelFullImage(obj, h2, t, f, s_abs, idx, 1:length(groups{idx}), groups{idx}, colors, file_name);
            end 

            close;
        end

        function saveLabelFullImage(obj, h, t, f, s, idx, jdx, groups, colors, file_name)
            
            %fig = figure;
            
            h.Children(end).CData = s;
            h.Children(end).XData = t;
            h.Children(end).YData = f;
            
            for loop = 1:length(jdx)
                if groups(loop).valid_flag == true
                    h.Children(loop).Visible = true;
                    h.Children(loop).LineWidth = 1;
                    %h.Children(1).MarkerSize = 1;
                    h.Children(loop).Color = [1 0 0]; % black
                    h.Children(loop).XData = groups(jdx(loop)).t_group;
                    %h.Children(loop).XData = groups(jdx(loop)).t_group(:,1,idx);
                    h.Children(loop).YData = groups(jdx(loop)).f_group;
                end
                
            end
            setImageSetting(obj, h, t, f);
            colormap('parula');
            im = getframe(h);

            imwrite(im.cdata, file_name);
            %close(fig);

            for loop = 1:length(groups)
                h.Children(loop).Visible = false;
            end

        end

        function saveLabelBinary(obj, s, file_name)

            fid = fopen(file_name, 'w');
            
            s = s.';
            s = flipud(s);
            %s = flipud(s);

            fwrite(fid, s(:), "uint8");
            fclose(fid);


        end

        function [ x_map ] = getBackwardMapping(obj, x_groups, x)
            dx = diff(x);
            dx = dx(1);

            x_valid = false(length(x_groups), 1);
            x_hit = zeros(length(x_groups), 1);

            for idx = 1:length(x_groups)

                [cal_val, hit_idx] = min(abs(x - x_groups(idx)));
                if cal_val <= dx
                    x_hit(idx) = hit_idx;
                    x_valid(idx) = true;
                    %break;
                else
                    x_valid(idx) = false;
                end
            end
            x_map = [];
            x_map.valid = x_valid;
            x_map.hit = x_hit;

        end

        function [ label_mat ] = getBackwardMappingMatrix(obj, t_groups, f_groups, t, f, s)
            
            target_size = size(s);
            t = linspace(t(1), t(end), target_size(1));
            f = linspace(f(1), f(end), target_size(2)); 

            label_mat = zeros(target_size);

            pair = [];
            pair.t = 0;
            pair.f = 0;
            pair.valid = false;
            pair = repmat(pair, [length(t_groups), 1]);

            [ t_map ] = getBackwardMapping(obj, t_groups, t);
            [ f_map ] = getBackwardMapping(obj, f_groups, f);

            for idx = 1:length(pair)
                if t_map.valid(idx) == true && f_map.valid(idx) == true
                    pair(idx).t = t_map.hit(idx);
                    pair(idx).f = f_map.hit(idx);
                    pair(idx).valid = true;
                    label_mat(pair(idx).t, pair(idx).f) = 1;
                end
            end
            
            %label_mat = imgaussfilt(flipud(rot90(label_mat)), 0.5);
            label_mat = imgaussfilt(label_mat, 0.1);
            
            label_mat(label_mat ~= 0) = 1;

        end


        function [num_classes, label_mat] = saveLabelImage(obj, h, t, f, s, idx, jdx, groups, colors, file_name)

            ext = char(file_name);
            ext = ext(end-2:end);

            child_idx = 0;
            for kdx = 1:length(h.Children)
                if strcmp(h.Children(kdx).Type, 'image')
                    child_idx = kdx;
                    break;
                end
            end
            h.Children(child_idx).CData = s;
            h.Children(child_idx).XData = t;
            h.Children(child_idx).YData = f;
            
            setImageSetting(obj, h, t, f);
            %colormap(h, colors.image);


            child_idx = 0;
            for kdx = 1:length(h.Children)
                if strcmp(h.Children(kdx).Type, 'line')
                    child_idx = kdx;
                    break;
                end
            end
            h.Children(child_idx).Visible = true;
            h.Children(child_idx).LineWidth = 1;
            h.Children(child_idx).MarkerSize = 1;
            h.Children(child_idx).Color = [0 0 0]; % black
            %h.Children(child_idx).Color = colors.plot;
            
            tm_map = uint8(zeros([size(s), length(jdx)]));
            tm_color_map = uint8(zeros([size(s), length(jdx)]));
            for loop=1:length(jdx)
                if groups(jdx(loop)).valid_flag == true || length(jdx) == 1
                    %t_groups = groups(jdx(loop)).t_group(:,:,idx);
                    t_groups = groups(jdx(loop)).t_group;
                    f_groups = groups(jdx(loop)).f_group;
                    h.Children(child_idx).XData = t_groups;
                    h.Children(child_idx).YData = f_groups;
                else
                    t_groups = zeros(length(groups(jdx(loop)).t_group), 1);
                    f_groups = zeros(length(groups(jdx(loop)).f_group), 1);
                    h.Children(child_idx).XData = t_groups;
                    h.Children(child_idx).YData = f_groups;
                end
                %h.Children(child_idx).Marker = colors.plot_marker;
                %h.Children(child_idx).LineStyle = colors.plot_style;

                %% backward mapping
                label_mat2 = zeros(size(s));
                for group_idx = 1:length(jdx)
                    t_groups = groups(jdx(group_idx)).t_group;
                    f_groups = groups(jdx(group_idx)).f_group;
                    [ tmp ] = getBackwardMappingMatrix(obj, t_groups, f_groups, t, f, s);
                    label_mat2 = label_mat2 + tmp;
                end
                


                %%

                fig_img = getframe(h);

                label = fig_img.cdata(2:end-1,2:end-1,:);
                %if strcmp(ext, 'gif')
                    label = rgb2gray(label);
                    tm_color_map(:,:,loop) = label;
                    label(label~=255) = 1;
                    label(label==255) = 0;
                    label = uint8(label);
                    tm_map(:,:,loop) = label;
                %end
                

            end
            label_mat = uint8(sum(tm_map, 3));
            num_classes = max(label_mat(:)) + 1;
            if num_classes == 2
                label_mat = logical(label_mat);
            end
        
            is_help_img = false;
            tokens = strsplit(file_name, '\');
            for loop = 1:length(tokens)
                if strcmp(tokens(loop), "overlay_mode") || strcmp(tokens(loop), "overlay_full")
                    is_help_img = true;
                    break;

                end
            end

            if is_help_img == true
                fig_img = getframe(h);
                imwrite(fig_img.cdata, file_name);

            else
                imwrite(label_mat, file_name);
            end
            
            h.Children(child_idx).Visible = false;

        end

        function [num_classes] = saveLabel(obj, h, t, f, s, idx, jdx, groups, colors, target_path, dir_name)
            
            data_index = num2str(obj.data_idx);
            [num_classes, label] = saveLabelImage(obj, h, t, f, s, idx, jdx, groups, colors, [target_path, '\', [dir_name '_img'], '\', data_index, '_mask.gif']);
            saveLabelBinary(obj, label, [target_path, '\', [dir_name '_bin'], '\', data_index, '_mask.bin']);


        end

        function setImageSetting(obj, h, t, f)
            %axis xy;
            set(h, 'YDir', 'normal');
            set(h, 'XTick', []);
            set(h, 'YTick', []);
            set(h,'Position',[0 0 1 1]) % Make the axes occupy the hole figure
            set(h, 'XLim', [min(t), max(t)]);
            set(h, 'YLim', [min(f), max(f)]);
            colormap(h, obj.cmap);
            caxis(h, [0 1]);
            
        end

        function saveImage(obj, h, t, f, s, file_name)
%             child_idx = 0;
%             for idx = 1:length(h.Children)
%                 if strcmp(h.Children(idx).Type, 'image')
%                     child_idx = idx;
%                     break;
%                 end
%             end
%             h.Children(child_idx).CData = s;
%             h.Children(child_idx).XData = t;
%             h.Children(child_idx).YData = f;
% 
%             setImageSetting(obj, h, t, f);
%             fig_img = getframe(h);
%             img = fig_img.cdata(2:end-1,2:end-1,:);
%             imwrite(rgb2gray(img), file_name);

            img = flipud(uint8(mat2gray(s)*255));
            imwrite(img, file_name);
            
            

        end

        function saveBinaryData(obj, s, t, f, file_name)

            fid = fopen(file_name, 'w');
            s = s.';
            s = rot90(s,2);
            fwrite(fid, s(:), "single");
            fclose(fid);


        end

        function saveJpgImage(obj, h, t, f, s, file_name, dir_name)
            
            
            data_index = num2str(obj.data_idx);
            saveImage(obj, h, t, f, s, [file_name, '\', [dir_name '_img'], '\', data_index, '.bmp']);
            saveBinaryData(obj, s, t, f, [file_name, '\', [dir_name '_bin'], '\', data_index, '.bin']);

            if strcmp(dir_name, 'abs')
                p = 0.00001;
                s_logp = abs(log(s + p));
                s_logp = (s_logp - min(s_logp(:))) /  ( max(s_logp(:)) - min(s_logp(:)));
                %[~,~,s_logp] = getStandardizedSpectrogram(obj, s_logp);

                saveBinaryData(obj, s_logp, t, f, [file_name, '\', [dir_name '_bin_log'], '\', data_index, '.bin']);
            end

            %img_now = imread(file_name);
            %img_now = img_now(2:end-1, 2:end-1,:);
            %imwrite(img_now, file_name);

        end

        function [s] = getInterpolatedMatrix(obj, s, t, f, target_size)
            [T, F] = meshgrid(t, f);
            [Tq, Fq] = meshgrid(linspace(t(1), t(end), target_size(1)), linspace(f(1), f(end), target_size(2)));
            s = interp2(T,F, s, Tq, Fq, 'spline');
        end

        function [s_norm] = getNormalizedMatrix(obj, s)
            s_norm = (s - min(s(:))) ./ (max(s(:)) - min(s(:)));

        end

        function [s_real, s_imag, s_abs] = getNormalizedSpectrogram(obj, s, t, f, target_size)

            s_real = getInterpolatedMatrix(obj, real(s), t, f, target_size);
            s_imag = getInterpolatedMatrix(obj, imag(s), t, f, target_size);
            s_abs = getInterpolatedMatrix(obj, abs(s), t, f, target_size);

            s_real = getNormalizedMatrix(obj, s_real);
            s_imag = getNormalizedMatrix(obj, s_imag);
            s_abs = getNormalizedMatrix(obj, s_abs);

            
            
        end

        function [s_real, s_imag, s_abs] = getStandardizedSpectrogram(obj, s)
            s_real = (real(s) - mean(real(s(:)))) / std(real(s(:)));
            s_imag = (imag(s) - mean(imag(s(:)))) / std(imag(s(:)));
            s_abs = (abs(s) - mean(abs(s(:)))) / std(abs(s(:)));

        end

        function saveSpectrogramFocusedTime(obj, mode_t, mode_spectro, groups, stft_pack, roi_freq)

        end


    end
end

