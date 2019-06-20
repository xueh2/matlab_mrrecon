
function [h_slc, h_obj_slc] = plotGTContoursOnImages(data, endo_pt, epi_pt, linewidth, visibility)
% [h_slc, h_obj_slc] = plotGTContoursOnImages(data, endo_pt, epi_pt, linewidth)

if (nargin<4)
    linewidth = 1;
end

if (nargin<5)
    visibility = 'off';
end

S = size(data);

scale_factor = 20;

RO = S(1);
E1 = S(2);
PHS = S(3);
if(numel(S)==3)
    SLC = 1;
else
    SLC = S(4);
end

max_d = max(data(:));

h_slc = zeros(SLC,1);
for slc = 1:SLC
    h_slc(slc) = figure('visible', visibility, 'Name', ['SLC ' num2str(slc)]);
    imagescn(data(:,:,:,slc), [0 0.3*max_d], [3 ceil(PHS/3)], scale_factor);

    h_axes=flipud(findobj(h_slc(slc),'type','axes'));

    for phs=1:PHS
        axes(h_axes(phs))

        if(~isempty(endo_pt))
            endo = find_contour_entry(endo_pt, phs, slc);
            if(~isempty(endo))
                endo = [endo; endo(1,:)];
                hold on  
                plot(endo(:,1)+1, endo(:,2)+1, 'r', 'LineWidth', linewidth);
                hold off
            end
        end

        if(~isempty(epi_pt))
            epi = find_contour_entry(epi_pt, phs, slc);
            if(~isempty(epi))
                epi = [epi; epi(1,:)];
                hold on    
                plot(epi(:,1)+1, epi(:,2)+1, 'g', 'LineWidth', linewidth);
                hold off
            end
        end
    end
end

% plot with objects
h_obj_slc = zeros(SLC,1);
for slc = 1:SLC
    h_obj_slc(slc) = figure('Name', ['SLC ' num2str(slc)], 'visible', visibility);
    imagescn(data(:,:,:,slc), [0 0.3*max_d], [], 6, 3);

    h_axes=flipud(findobj(h_obj_slc(slc),'type','axes'));

    clear ObjStruct
    for time = 1:PHS;

        if(~isempty(endo_pt))
            endo = find_contour_entry(endo_pt, time, slc);
            if(isempty(endo))
                endo = [0 0];
            else
                endo = [endo; endo(1,:)];
            end
        else
            endo = [0 0];
        end

        if(~isempty(epi_pt))
            epi = find_contour_entry(epi_pt, time, slc);
            if(isempty(epi))
                epi = [0 0];
            else
                epi = [epi; epi(1,:)];
            end
        else
            epi = [0 0];
        end

        ObjStruct(1,time).name = 'Endo';
        ObjStruct(1,time).type = 'Line';
        ObjStruct(1,time).xdata = endo(:,1);
        ObjStruct(1,time).ydata = endo(:,2);
        ObjStruct(1,time).color = [1 0 0]; % red
        ObjStruct(1,time).marker = '.';
        ObjStruct(1,time).markersize=2;
        ObjStruct(1,time).markerfacecolor=[0 1 0];

        ObjStruct(2,time).name = 'Epi';
        ObjStruct(2,time).type = 'Line';
        ObjStruct(2,time).xdata = epi(:,1);
        ObjStruct(2,time).ydata = epi(:,2);
        ObjStruct(2,time).color = [0 1 0]; % green
        ObjStruct(2,time).marker = '.';
        ObjStruct(2,time).markersize=2;
        ObjStruct(2,time).markerfacecolor=[0 1 0];
    end

    setappdata(h_axes,'Objects',ObjStruct);
end
end

function endo = find_contour_entry(endo_pt, phs, slc)

N = size(endo_pt, 1);

endo = [];
for ii=1:N
    if(endo_pt{ii, 1} == phs-1 & endo_pt{ii, 2} == slc-1)
        endo = endo_pt{ii, 3};
        break;
    end
end

return
end


