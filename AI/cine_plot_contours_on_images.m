
function [h_slc, h_ED, h_ES, h_slc_obj] = cine_plot_contours_on_images(Cine, endoC, epiC, linewidth, ED, ES)

if (nargin<4)
    linewidth = 1;
end

if (nargin<5)
    ED = 1;
end

if (nargin<6)
    ES = 13;
end

S = size(Cine);

scale_factor = 20;

RO = S(1);
E1 = S(2);
PHS = S(3);
SLC = S(4);

% plot ED, ES

h_ED = plot_ED_ES(Cine, endoC, epiC, ED, 'ED', scale_factor, linewidth);
h_ES = plot_ED_ES(Cine, endoC, epiC, ES, 'ES', scale_factor, linewidth);

h_slc = zeros(SLC,1);
for slc=1:SLC

    h = figure('Name', ['SLC ' num2str(slc)], 'visible', 'off');
    imagescn(Cine(:,:,:,slc), [], [3 ceil(PHS/3)], scale_factor);

    h_axes=flipud(findobj(h,'type','axes'));
    
    for phs=1:PHS
        axes(h_axes(phs))
           
        if(~isempty(endoC))
            endo = endoC{slc, phs};            
            if(~isempty(endo))
                endo = [endo; endo(1,:)];
                hold on  
%                 if(isunix())
%                     plot(endo(:,2)+1, endo(:,1)+1, 'r', 'LineWidth', linewidth);
%                 else
                    plot(endo(:,1)+1, endo(:,2)+1, 'r', 'LineWidth', linewidth);
%                 end
                hold off
            end
        end
        
        if(~isempty(epiC))
            epi = epiC{slc, phs};            
            if(~isempty(epi))
                epi = [epi; epi(1,:)];
                hold on    
%                 if(isunix())
%                     plot(epi(:,2)+1, epi(:,1)+1, 'g', 'LineWidth', linewidth);
%                 else
                    plot(epi(:,1)+1, epi(:,2)+1, 'g', 'LineWidth', linewidth);
%                 end
                hold off
            end
        end
    end
    
    h_slc(slc) = h;
end

% plot with objects
h_slc_obj = zeros(SLC,1);
for slc=1:SLC

    % h = figure('Name', ['SLC ' num2str(slc)], 'visible', 'off');
    h = figure('Name', ['SLC ' num2str(slc)]);
    imagescn(Cine(:,:,:,slc), [], [], 6, 3);

    h_axes=flipud(findobj(h,'type','axes'));
    
    clear ObjStruct
    for time = 1:PHS;

        if(~isempty(endoC))
            endo = endoC{slc, time};
            if(isempty(endo))
                endo = [0 0];
            else
                endo = [endo; endo(1,:)];
            end
        else
            endo = [0 0];
        end
        
        if(~isempty(epiC))
            epi = epiC{slc, time};
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
    
    h_slc_obj(slc) = h;
end

end

function h = plot_ED_ES(Cine, endoC, epiC, ED, phs_name, scale_factor, linewidth)

    PHS = size(Cine, 3);
    SLC = size(Cine, 4);

    h = figure('Name', [phs_name ' Plot']);
    
    ind = [ED:PHS 1:ED-1];
    
    CineSelected = Cine(:,:, ind, :);
    imagescn(CineSelected, [], [3 ceil(SLC/3)], scale_factor, 3);
    % imagescn(CineSelected, [], [1 SLC], scale_factor, 3);
    h_axes=flipud(findobj(h,'type','axes'));

    for slc=1:SLC    

        axes(h_axes(slc))

        if(~isempty(endoC))
            endo = endoC{slc, ED};
            if(~isempty(endo))
                endo = [endo; endo(1,:)] + 1;
                hold on    
                plot(endo(:,1), endo(:,2), 'r', 'LineWidth', linewidth);
                hold off
            end
        end

        if(~isempty(epiC))
            epi = epiC{slc, ED};
            if(~isempty(epi))
                epi = [epi; epi(1,:)] + 1;
                hold on    
                plot(epi(:,1), epi(:,2), 'g', 'LineWidth', linewidth);
                hold off
            end
        end
    end
end


