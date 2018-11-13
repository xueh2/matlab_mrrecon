function [im, endo, epi, rv, rv_insertion] = get_rois_from_manual_labelled_sax(h_axes, RO, E1)
% [im, endo, epi, rv, rv_insertion] = get_rois_from_manual_labelled_sax(h_axes, RO, E1)
% endo/epi are 1-based, not 0-based

h_image=findobj(h_axes(1),'type','image');
im = get(h_image(1),'cdata');

h_lines=findobj(h_axes(1),'type','line');

rv_insertion = [];
endo = [];
epi = [];
rv = [];
C1 = endo;
C2 = epi;
C3 = rv;

fmap0 = zeros(RO, E1);

for k=1:numel(h_lines)
    X = get(h_lines(k),'XData');
    Y = get(h_lines(k),'YData');

%     numel(X)
    
    if(numel(X)==2)
        rv_insertion = [X' Y'];
    end

    if(numel(X)>20)
        if (isempty(C1))
            C1 = [X' Y'];
        else
            if (isempty(C2)) 
                C2 = [X' Y'];
            else
                if (isempty(C3)) 
                    C3 = [X' Y'];
                end
            end
        end
    end
end

if(~isempty(C1) & isempty(C2) & isempty(C3))
    BW1=roipoly(fmap0,C1(:,1), C1(:,2));

    rvi = round(rv_insertion);

    inBW1 = 0;        
    if( BW1(rvi(1,2), rvi(1,1)) | BW1(rvi(2,2), rvi(2,1)) )
        inBW1 = 1;
    end

    if(inBW1)
        epi = C1;
        rv = [];
        endo = [];
    end

    if(~inBW1)
        epi = [];
        rv = C1;
        endo = [];
    end
end

if(~isempty(C1) & ~isempty(C2) & isempty(C3))

    BW1=roipoly(fmap0,C1(:,1), C1(:,2));
    BW2=roipoly(fmap0,C2(:,1), C2(:,2));

    if(isempty(rv_insertion))
        if(size(C1,1)>size(C2, 1))
            endo = C2;
            epi = C1;
        else
            endo = C1;
            epi = C2;
        end
    else        
        rvi = round(rv_insertion);

        inBW1 = 0;        
        if( BW1(rvi(1,2), rvi(1,1)) | BW1(rvi(2,2), rvi(2,1)) )
            inBW1 = 1;
        end
        inBW2 = 0;        
        if( BW2(rvi(1,2), rvi(1,1)) | BW2(rvi(2,2), rvi(2,1)) )
            inBW2 = 1;
        end

        if(inBW1 & inBW2)
            rv = [];
            if(size(C1,1)>size(C2, 1))
                endo = C2;
                epi = C1;
            else
                endo = C1;
                epi = C2;
            end                
        end

        if(inBW1 & ~inBW2)
            rv = C2;
            endo = [];
            epi = C1;
        end
        if(~inBW1 & inBW2)
            rv = C1;
            endo = [];
            epi = C2;
        end
    end
end

if(~isempty(C1) & ~isempty(C2) & ~isempty(C3))       
    % check endo/epi
    BW1=roipoly(fmap0,C1(:,1), C1(:,2));
    BW2=roipoly(fmap0,C2(:,1), C2(:,2));
    BW3=roipoly(fmap0,C3(:,1), C3(:,2));

    rvi = round(rv_insertion);

    inBW1 = 0;        
    if( BW1(rvi(1,2), rvi(1,1)) | BW1(rvi(2,2), rvi(2,1)) )
        inBW1 = 1;
    end
    inBW2 = 0;        
    if( BW2(rvi(1,2), rvi(1,1)) | BW2(rvi(2,2), rvi(2,1)) )
        inBW2 = 1;
    end
    inBW3 = 0;        
    if( BW3(rvi(1,2), rvi(1,1)) | BW3(rvi(2,2), rvi(2,1)) )
        inBW3 = 1;
    end
%     disp(['inBW1 = ' num2str(inBW1) ' inBW2 = ' num2str(inBW2) ' inBW3 = ' num2str(inBW3)]);

    if(~inBW1 & inBW2 & inBW3)
        rv = C1;

        if(size(C2,1)>size(C3, 1))
            endo = C3;
            epi = C2;
        else
            endo = C2;
            epi = C3;
        end
    end

    if(inBW1 & ~inBW2 & inBW3)
        rv = C2;

        if(size(C1,1)>size(C3, 1))
            endo = C3;
            epi = C1;
        else
            endo = C1;
            epi = C3;
        end
    end

    if(inBW1 & inBW2 & ~inBW3)
        rv = C3;

        if(size(C1,1)>size(C2, 1))
            endo = C2;
            epi = C1;
        else
            endo = C1;
            epi = C2;
        end
    end            
end

