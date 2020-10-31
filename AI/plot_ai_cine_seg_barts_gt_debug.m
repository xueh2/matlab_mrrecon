function plot_ai_cine_seg_barts_gt_debug(I_sax, endo_mask, epi_mask, endo_C, epi_C, I_2ch, I_4ch, ch2_pts, ch4_pts, metrics)
% plot_ai_cine_seg_barts_gt_debug(I_sax, endo_mask, epi_mask, endo_C, epi_C, I_2ch, I_4ch, ch2_pts, ch4_pts, metrics)

linewidth = 2.0;
cut_thres = [0.8 1.0];

size(I_2ch)
size(I_4ch)
size(I_sax)

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

I_sax_heart = I_sax;
figure; imagescn(I_sax_heart, [], [], [], 3);

endo_heart = endo_mask;
figure('Name', 'Endo mask'); imagescn(endo_heart, [], [], [], 3);
epi_heart = epi_mask;
figure('Name', 'EPI mask'); imagescn(epi_heart, [], [], [], 3);

% plot ED
h_ED = figure('Name', [' ED ' num2str(ED)], 'visible', 'on');
imm = I_sax_heart(:,:,ED,:);
imagescn(I_sax_heart(:,:,ED,:), [0 0.45*max(imm(:))], [3 ceil(SLC/3)], 12);

h_axes=flipud(findobj(h_ED,'type','axes'));

need_cut_ED = zeros(SLC, 1);
wrong_valve_cut_ED = zeros(SLC, 1);

if(size(endo_C, 3)==PHS)
    endo_C = endo_C(:,:,[ED, ES],:);
end

if(size(epi_C, 3)==PHS)
    epi_C = epi_C(:,:,[ED, ES],:);
end

s_e1 = 0;
s_ro = 0;

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
                plot(endo(:,2)+1, endo(:,1)+1, 'r', 'LineWidth', linewidth);
                hold off
            end
        else
            endo = [endo; endo(1,:)];
            hold on  
            plot(endo(:,2)+1-s_e1, endo(:,1)+1-s_ro, 'r', 'LineWidth', linewidth);
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
                plot(epi(:,2)+1, epi(:,1)+1, 'g', 'LineWidth', linewidth);
                hold off
            end
        else
            epi = [epi; epi(1,:)];
            hold on  
            plot(epi(:,2)+1-s_e1, epi(:,1)+1-s_ro, 'g', 'LineWidth', linewidth);
            hold off
        end 
    end
end

% plot ES
h_ES = figure('Name', [' ES ' num2str(ES)], 'visible', 'on');
imm = I_sax_heart(:,:,ES,:);
imagescn(I_sax_heart(:,:,ES,:), [0 0.45*max(imm(:))], [3 ceil(SLC/3)], 12);

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
            plot(endo(:,2)+1, endo(:,1)+1, 'r', 'LineWidth', linewidth);
            hold off
        else
            endo = [endo; endo(1,:)];
            hold on  
            plot(endo(:,2)+1-s_e1, endo(:,1)+1-s_ro, 'r', 'LineWidth', linewidth);
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

I_2ch_heart = I_2ch;

h_2ch = figure('Name', ['ch 2'], 'visible', 'on');
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
I_4ch_heart = I_4ch;

h_4ch = figure('Name', ['ch 4'], 'visible', 'on');
imagescn(I_4ch_heart(:,:,[ED ES]), [], [1 2], 12);

h_axes=flipud(findobj(h_4ch,'type','axes'));

for slc=1:2
    axes(h_axes(slc))

    hold on  
    plot(ch4_pts(slc, 1, 1)+1-s_e1, ch4_pts(slc, 1, 2)+1-s_ro, 'r+', 'MarkerSize', 20, 'LineWidth', 2);
    plot(ch4_pts(slc, 2, 1)+1-s_e1, ch4_pts(slc, 2, 2)+1-s_ro, 'r+', 'MarkerSize', 20, 'LineWidth', 2);
    hold off
end
