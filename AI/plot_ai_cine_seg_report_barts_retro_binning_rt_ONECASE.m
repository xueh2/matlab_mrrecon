function plot_ai_cine_seg_report_barts_retro_binning_rt_ONECASE(dst_dir, aiDir, pt_id, roi_size, cut_thres)
% plot_ai_cine_seg_report_barts_retro_binning_rt_ONECASE(dstDir, aiDir)

linewidth = 4.0;

I_2ch = readNPY(fullfile(dst_dir, 'I_2ch_1mm.npy'));
I_4ch = readNPY(fullfile(dst_dir, 'I_4ch_1mm.npy'));
I_sax = readNPY(fullfile(dst_dir, 'I_sax_1mm.npy'));

size(I_2ch)
size(I_4ch)
size(I_sax)

metrics = analyze75read(fullfile(aiDir, 'metrics_1mm.hdr'));
endo_mask = analyze75read(fullfile(aiDir, 'sax_endo_mask_1mm.hdr'));
epi_mask = analyze75read(fullfile(aiDir, 'sax_epi_mask_1mm.hdr'));
endo_C = analyze75read(fullfile(aiDir, 'sax_endo_contour_1mm.hdr'));
epi_C = analyze75read(fullfile(aiDir, 'sax_epi_contour_1mm.hdr'));
ch2_pts = analyze75read(fullfile(aiDir, 'ch2_pts_1mm.hdr'));
ch4_pts = analyze75read(fullfile(aiDir, 'ch4_pts_1mm.hdr'));

ED = metrics(1)+1;
ES = metrics(2)+1;        
EF = metrics(3);
SV = metrics(4);
EDV = metrics(5);
ESV = metrics(6);
MASS = metrics(7);

RO = size(I_sax, 1);
E1 = size(I_sax, 2);
PHS = size(I_sax, 3);
SLC = size(I_sax, 4);

if(roi_size>E1)
    roi_size = min([RO, E1])
end
if(roi_size>RO)
    roi_size = min([RO, E1]);
end

endo_mask_sum = sum(sum(endo_mask, 3), 4);
if(sum(endo_mask_sum(:))<1)
    return
end

binaryImage = true(size(endo_mask_sum));
labeledImage = bwlabel(binaryImage);
measurements = regionprops(labeledImage, endo_mask_sum, 'WeightedCentroid')
centerOfMass = measurements.WeightedCentroid;

c_e1 = centerOfMass(1);
c_ro = centerOfMass(2);

s_ro = round(c_ro - roi_size/2);
s_e1 = round(c_e1 - roi_size/2);

if(s_ro<1)
    s_ro = 1;
end
if(s_e1<1)
    s_e1 = 1;
end

e_ro = s_ro+roi_size-1;
e_e1 = s_e1+roi_size-1;

if(e_ro>=RO)
    e_ro = RO;
    s_ro = e_ro-roi_size+1;
end
if(e_e1>=E1)
    e_e1 = E1;
    s_e1 = e_e1-roi_size+1;
end

I_sax_heart = I_sax(s_ro:e_ro, s_e1:e_e1, :,:);
figure; imagescn(I_sax_heart, [], [], [], 3);

endo_heart = endo_mask(s_ro:e_ro, s_e1:e_e1, :,:);
figure; imagescn(endo_heart, [], [], [], 3);
epi_heart = epi_mask(s_ro:e_ro, s_e1:e_e1, :,:);
figure; imagescn(epi_heart, [], [], [], 3);

% plot ED
h_ED = figure('Name', [pt_id ' ED ' num2str(ED)], 'visible', 'on');
imm = I_sax_heart(:,:,ED,:);
imagescn(I_sax_heart(:,:,ED,:), [0 0.45*max(imm(:))], [3 ceil(SLC/3)], 16);

h_axes=flipud(findobj(h_ED,'type','axes'));

need_cut_ED = zeros(SLC, 1);
wrong_valve_cut_ED = zeros(SLC, 1);

for slc=1:SLC
    axes(h_axes(slc))

    endo_slc_mask = endo_heart(:,:,1,slc);
    epi_slc_mask = epi_heart(:,:,1,slc);

    RO = size(epi_slc_mask, 1);
    E1 = size(epi_slc_mask, 2);
    
    ind = find(endo_slc_mask>0);
    if(numel(ind)<10)
        need_cut_ED(slc) = 1;
        continue;
    end

    ind = find(epi_slc_mask==0);
    endo_slc_mask(ind) = 0;

    cut_endo = 0;

    need_cut = 0;
    max_epi_slc_mask = max(epi_slc_mask(:));
    ind_all = find(epi_slc_mask>0);
    ind = find(epi_slc_mask>0.8*max_epi_slc_mask);
    if(max_epi_slc_mask<0.5 | numel(ind)<0.1*numel(ind_all))
        need_cut_ED(slc) = 1;
        continue;
    end
    
    if(numel(ind)>0.98*numel(ind_all))
        need_cut = 0;
    else
        need_cut = 1;
    end
    
    need_cut_ED(slc) = need_cut;
    
    if(need_cut)
        % if cut portion is on the wrong side of valve, do not cut
        ind_all = find(epi_slc_mask>0);
        
        max_prob = max(epi_slc_mask(:));
        [ind_high_y, ind_high_x] = find(epi_slc_mask>0 & epi_slc_mask>0.5*max_prob);
        [ind_low_y, ind_low_x] = find(epi_slc_mask>0 & epi_slc_mask<0.5*max_prob);
        
        if(numel(ind_low_y)<0.1*numel(ind_all))
            need_cut = 0;
            wrong_valve_cut_ED(slc) = 1;
        else
            c_high = [mean(ind_high_x), mean(ind_high_y)];
            c_low = [mean(ind_low_x), mean(ind_low_y)];

            hold on  
            plot(c_low(1), c_low(2), 'r+', 'MarkerSize', 12);
            plot(c_high(1), c_high(2), 'b+', 'MarkerSize', 12);
            hold off

            if(c_low(1)>c_high(1))
                disp(['wrong valve cuting detected ... ']);
                need_cut = 0;
                wrong_valve_cut_ED(slc) = 1;
            end
        end
    end
    
    % plot endo
    endo = endo_C(:,:,1,slc);
    ind = find(endo<0);
    if(numel(ind)==0)

        ind_endo = find(endo_slc_mask>=cut_thres(1) & endo_slc_mask<=cut_thres(2));
        if(need_cut & numel(ind_endo)>20)
            endo = mask2contour(endo_slc_mask, cut_thres, 400, 16);
            if(~isempty(endo))
                endo = [endo; endo(1,:)];
                cut_endo = 1;
                hold on  
                plot(endo(:,1)+1, endo(:,2)+1, 'r', 'LineWidth', linewidth);
                hold off
            end
        else
            endo = [endo; endo(1,:)];
            hold on  
            plot(endo(:,1)+1-s_e1, endo(:,2)+1-s_ro, 'r', 'LineWidth', linewidth);
            hold off
        end                                
    end

    epi = epi_C(:,:,1,slc);
    ind = find(epi<0);
    if(numel(ind)==0)
        if(need_cut & numel(ind_endo)>20)
            epi = mask2contour(epi_slc_mask, cut_thres, 400, 16);
            if(~isempty(epi))
                epi = [epi; epi(1,:)];
                hold on  
                plot(epi(:,1)+1, epi(:,2)+1, 'g', 'LineWidth', linewidth);
                hold off
            end
        else
            epi = [epi; epi(1,:)];
            hold on  
            plot(epi(:,1)+1-s_e1, epi(:,2)+1-s_ro, 'g', 'LineWidth', linewidth);
            hold off
        end 
    end
end

% plot ES
h_ES = figure('Name', [pt_id ' ES ' num2str(ES)], 'visible', 'on');
imm = I_sax_heart(:,:,ES,:);
imagescn(I_sax_heart(:,:,ES,:), [0 0.45*max(imm(:))], [3 ceil(SLC/3)], 16);

h_axes=flipud(findobj(h_ES,'type','axes'));

cut_thres_ES = cut_thres;

for slc=1:SLC
    axes(h_axes(slc))

    epi_slc_mask = epi_heart(:,:,1,slc);
    endo_slc_mask = endo_heart(:,:,2,slc);

    ind = find(epi_slc_mask>0);
    if(numel(ind)<10)
        continue;
    end
    
    ind = find(endo_slc_mask>0);
    if(numel(ind)<10)
        continue;
    end

    need_cut = 0;
    max_epi_slc_mask = max(epi_slc_mask(:));
    ind_all = find(epi_slc_mask>0);
    ind = find(epi_slc_mask>0.8*max_epi_slc_mask);
    if(max_epi_slc_mask<0.5 | numel(ind)<0.1*numel(ind_all))
        need_cut_ED(slc) = 1;
        continue;
    end
    
    if(numel(ind)>0.98*numel(ind_all))
        need_cut = 0;
    else
        need_cut = 1;
    end
    
    if(need_cut==0)
        ind_all = find(endo_slc_mask>0);
        ind = find(endo_slc_mask==1);
        if(numel(ind)<0.5*numel(ind_all))
            need_cut = 1;
        end
    end
    
    if(need_cut_ED(slc))
        need_cut = 1;
    end
    
    if(wrong_valve_cut_ED(slc))
        need_cut = 0;
    elseif (need_cut)
        ind_all = find(endo_slc_mask>0);
        
        max_prob = max(endo_slc_mask(:));
        [ind_high_y, ind_high_x] = find(endo_slc_mask>0 & endo_slc_mask>0.5*max_prob);
        [ind_low_y, ind_low_x] = find(endo_slc_mask>0 & endo_slc_mask<0.5*max_prob);
        
        c_high = [mean(ind_high_x), mean(ind_high_y)];
        c_low = [mean(ind_low_x), mean(ind_low_y)];
        
        hold on  
        plot(c_low(1), c_low(2), 'r+', 'MarkerSize', 12);
        plot(c_high(1), c_high(2), 'b+', 'MarkerSize', 12);
        hold off
        
        if(c_low(1)>c_high(1) | (abs(c_low(1)-c_high(1))<5) )
            disp(['wrong valve cuting detected ... ']);
            need_cut = 0;
        end
    end
    
    % plot endo
    endo = endo_C(:,:,2,slc);
    ind = find(endo<0);
    ind_endo = find(endo_slc_mask>=cut_thres_ES(1) & endo_slc_mask<=cut_thres_ES(2));
    if(numel(ind)==0)
        if(need_cut & numel(ind_endo)>20)
            endo = mask2contour(endo_slc_mask, cut_thres_ES, 400, 16);
            if(isempty(endo))
                continue;
            end
            
            ind = find(endo_slc_mask>0);
            ind_cut = find(endo_slc_mask>=cut_thres_ES(1) & endo_slc_mask<=cut_thres_ES(2));
            if(numel(ind_cut) < 0.5*numel(ind))
                continue;
            end
            
            endo = [endo; endo(1,:)];
            hold on  
            plot(endo(:,1)+1, endo(:,2)+1, 'r', 'LineWidth', linewidth);
            hold off
        else
            endo = [endo; endo(1,:)];
            hold on  
            plot(endo(:,1)+1-s_e1, endo(:,2)+1-s_ro, 'r', 'LineWidth', linewidth);
            hold off
        end 
    end

%             epi = epi_C(:,:,2,slc);
%             ind = find(epi<0);
%             if(numel(ind)==0)
%                 epi = [epi; epi(1,:)];
%                 hold on  
%                 plot(epi(:,1)+1-s_e1, epi(:,2)+1-s_ro, 'g', 'LineWidth', linewidth);
%                 hold off
%             end
end

% 2CH, 4CH

roi_size = 256;

ch2_pts_mean = squeeze(mean(mean(ch2_pts, 1), 2));
c_e1 = ch2_pts_mean(1);
c_ro = ch2_pts_mean(2);

s_ro = round(c_ro - roi_size/2);
s_e1 = round(c_e1 - roi_size/2);

if(s_ro<1)
    s_ro = 1;
end
if(s_e1<1)
    s_e1 = 1;
end

e_ro = s_ro+roi_size-1;
e_e1 = s_e1+roi_size-1;

if(e_ro>size(I_2ch, 1))
    e_ro = size(I_2ch, 1);
    s_ro = e_ro-roi_size+1;
    
    if(s_ro<1)
        s_ro = 1;
        e_ro = s_ro+roi_size-1;
    end
end

if(e_e1>size(I_2ch, 2))
    e_e1 = size(I_2ch, 2);
    s_e1 = e_e1-roi_size+1;
    
    if(s_e1<1)
        s_e1 = 1;
        e_e1 = s_e1+roi_size-1;
    end
end

if(e_e1>size(I_2ch,2))
    e_e1 = size(I_2ch,2);
end

I_2ch_heart = I_2ch(s_ro:e_ro, s_e1:e_e1, :);

h_2ch = figure('Name', [pt_id], 'visible', 'on');
imagescn(I_2ch_heart(:,:,[ED ES]), [], [1 2], 12);

h_axes=flipud(findobj(h_2ch,'type','axes'));

for slc=1:2
    axes(h_axes(slc))

    hold on  
    plot(ch2_pts(slc, 1, 1)+1-s_e1, ch2_pts(slc, 1, 2)+1-s_ro, 'r+', 'MarkerSize', 20, 'LineWidth', 2);
    plot(ch2_pts(slc, 2, 1)+1-s_e1, ch2_pts(slc, 2, 2)+1-s_ro, 'r+', 'MarkerSize', 20, 'LineWidth', 2);
    hold off
end

% 4ch
ch4_pts_mean = squeeze(mean(mean(ch4_pts, 1), 2));
c_e1 = ch4_pts_mean(1);
c_ro = ch4_pts_mean(2);

s_ro = round(c_ro - roi_size/2);
s_e1 = round(c_e1 - roi_size/2);

if(s_ro<1)
    s_ro = 1;
end
if(s_e1<1)
    s_e1 = 1;
end

e_ro = s_ro+roi_size-1;
e_e1 = s_e1+roi_size-1;

if(e_ro>size(I_4ch, 1))
    e_ro = size(I_4ch, 1);
    s_ro = e_ro-roi_size+1;
    
    if(s_ro<1)
        s_ro = 1;
        e_ro = s_ro+roi_size-1;
    end
end

if(e_ro>size(I_4ch,1))
    e_ro = size(I_4ch,1);
end

if(e_e1>size(I_4ch, 2))
    e_e1= size(I_4ch, 2);
    s_e1 = e_e1-roi_size+1;
    
    if(s_e1<1)
        s_e1 = 1;
        e_e1 = s_e1+roi_size-1;
    end
end

if(e_e1>size(I_4ch,2))
    e_e1 = size(I_4ch,2);
end

I_4ch_heart = I_4ch(s_ro:e_ro, s_e1:e_e1, :);        

h_4ch = figure('Name', [pt_id], 'visible', 'on');
imagescn(I_4ch_heart(:,:,[ED ES]), [], [1 2], 12);

h_axes=flipud(findobj(h_4ch,'type','axes'));

for slc=1:2
    axes(h_axes(slc))

    hold on  
    plot(ch4_pts(slc, 1, 1)+1-s_e1, ch4_pts(slc, 1, 2)+1-s_ro, 'r+', 'MarkerSize', 20, 'LineWidth', 2);
    plot(ch4_pts(slc, 2, 1)+1-s_e1, ch4_pts(slc, 2, 2)+1-s_ro, 'r+', 'MarkerSize', 20, 'LineWidth', 2);
    hold off
end

% generate report
RowNames = {'EF', 'SV', 'EDV', 'ESV', 'Mass'};
ColumnName = {'Value', 'Unit'}
T = table(EF, SV, EDV, ESV, MASS);
Data = {EF*100, '    %'; SV, '    ml'; EDV, '    ml'; ESV, '    ml'; MASS, '    g'}
h_metric = figure; uitable('Data',Data,'ColumnName',ColumnName,'RowName',RowNames,'fontSize',12);

saveas(h_ES, fullfile(dst_dir, 'ES'), 'jpg');
saveas(h_ES, fullfile(dst_dir, 'ES'), 'fig');
saveas(h_ED, fullfile(dst_dir, 'ED'), 'jpg');
saveas(h_ED, fullfile(dst_dir, 'ED'), 'fig');
saveas(h_4ch, fullfile(dst_dir, 'Ch4'), 'jpg');
saveas(h_4ch, fullfile(dst_dir, 'Ch4'), 'fig');
saveas(h_2ch, fullfile(dst_dir, 'Ch2'), 'jpg');
saveas(h_2ch, fullfile(dst_dir, 'Ch2'), 'fig');
saveas(h_metric, fullfile(dst_dir, 'Measure'), 'jpg');
saveas(h_metric, fullfile(dst_dir, 'Measure'), 'fig');
save(fullfile(dst_dir, 'Measure'), 'T', 'metrics');

