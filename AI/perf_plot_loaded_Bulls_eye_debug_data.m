
function [h_fmap, h_Gd, h_mask] = perf_plot_loaded_Bulls_eye_debug_data(res, permute_data)
% plot loaded bulls eye data
% [h_fmap, h_mask] = perf_plot_loaded_Bulls_eye_debug_data(res)

if(nargin<2)
    permute_data = 1;
end

if(permute_data)
    fmap = permute(res.fmap, [2 1 3]);
    perf_for_seg = permute(res.perf_for_seg, [2 1 3 4]);
    endo_epi_rv_rvi_mask = permute(res.endo_epi_rv_rvi_mask, [2 1 3]);
else
    fmap = res.fmap;
    perf_for_seg = res.perf_for_seg;
    endo_epi_rv_rvi_mask = res.endo_epi_rv_rvi_mask;
end

h_fmap = figure; imagescn(fmap, [0 8], [1 size(res.fmap,3)], 12); PerfColorMap;

h_axes=flipud(findobj(h_fmap,'type','axes'));
plot_sectors(res, h_axes, permute_data);

h_Gd = figure; imagescn(perf_for_seg, [0 1.5], [1 size(res.perf_for_seg, 4)], 12, 3);
h_axes=flipud(findobj(h_Gd,'type','axes'));
plot_sectors(res, h_axes, permute_data);

h_mask = figure; imagescn(endo_epi_rv_rvi_mask, [], [1 size(res.perf_for_seg, 4)], 12);

end

function plot_sectors(res, h_axes, permute_data)
    SLC = numel(h_axes);
    for slc=1:SLC
        axes(h_axes(slc))

        if(slc==1)
            sectors = [1:6];
        end
        if(slc==2)
            sectors = [7:12];
        end
        if(slc==3)
            sectors = [13:16];
        end

        hold on    
        for k=1:numel(sectors)
            C = res.sectors_contours(:,:,sectors(k));
            lineW = 3;
            colorC = 'r';
%             if(k==1 | k==2)
%                 lineW = 4;
%                 colorC = 'y';
%                 if(k==2)
%                     colorC = 'g';
%                 end
%             end
                            
            % LAD
            if(sectors(k)==1 | sectors(k)==2 | sectors(k)==7 | sectors(k)==8 | sectors(k)==13  | sectors(k)==14)
                lineW = 3;
                colorC = 'y';
            end

            % RCA
            if(sectors(k)==3 | sectors(k)==4 | sectors(k)==9 | sectors(k)==10 | sectors(k)==15)
                lineW = 3;
                colorC = 'g';
            end
                        
            if(permute_data)
                plot(C(:,1), C(:,2), colorC, 'LineWidth', lineW);    
            else
                plot(C(:,2), C(:,1), colorC, 'LineWidth', lineW);    
            end
        end
        hold off
    end
end