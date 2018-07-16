
function [res, endo_res, epi_res] = SplitEndoEpiContourForAHAModel(endo, epi, lv_center, rv_insertion, num_seg, plot_flag)
% [res, endo, epi] = SplitEndoEpiContourForAHAModel(endo, epi, lv_center, rv_insertion, num_seg)
%
% res are the contours for num_seg
%
% endo and epi are the contours split transmurally into endo and epi

if(nargin<6)
    plot_flag = 0;
end

res_endo = SplitContourForAHAModel(endo, lv_center, rv_insertion, num_seg);
res_epi = SplitContourForAHAModel(epi, lv_center, rv_insertion, num_seg);

if(numel(res_endo)~=numel(res_epi))
    res = [];
    return;
end

res = cell(numel(res_endo), 1);

for ii=1:numel(res_endo)
    c1 = res_endo{ii};
    c2 = res_epi{ii};
    
    c3 = [c1; c2(end:-1:1, :); c1(1,:)];
    res{ii} = c3;
end

if nargout > 1 % compute endo and epi contours
    % interpolate c1 and c2 to have same number of samples and compute mid
    % points
    for ii=1:numel(res_endo)
        c1 = res_endo{ii};
        c2 = res_epi{ii};
        c1x = c1(:,1); c1y = c1(:,2);
        c1x_interp = interp1(linspace(0,1,length(c1x)),c1x, linspace(0,1,100));
        c1y_interp = interp1(linspace(0,1,length(c1y)),c1y, linspace(0,1,100));
        c1_interp = [c1x_interp(:), c1y_interp(:)];
        c2x = c2(:,1); c2y = c2(:,2);
        c2x_interp = interp1(linspace(0,1,length(c2x)),c2x, linspace(0,1,100));
        c2y_interp = interp1(linspace(0,1,length(c2y)),c2y, linspace(0,1,100));
        c2_interp = [c2x_interp(:), c2y_interp(:)];
        mid = (c1_interp + c2_interp)/2;
        
        endo_res{ii} = [c1_interp; mid(end:-1:1, :); c1_interp(1,:)];
        epi_res{ii} = [c2_interp; mid(end:-1:1, :); c2_interp(1,:)];
    end
end

if(plot_flag)
    cmap = hsv(numel(res));
    
    figure;
    hold on
    for ii=1:numel(res)
        C3 = res{ii};
        plot(C3(:,1), C3(:,2), '-', 'Color',cmap(ii,:));
    end
    if nargout> 1
        for ii=1:numel(res)
            C3 = endo_res{ii};
            plot(C3(:,1), C3(:,2), '-', 'Color',cmap(ii,:));
            C3 = epi_res{ii};
            plot(C3(:,1), C3(:,2), '-', 'Color',cmap(ii,:));
        end        
    end
    hold off
end

