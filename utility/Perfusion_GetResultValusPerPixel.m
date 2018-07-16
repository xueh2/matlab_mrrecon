function pt = Perfusion_GetResultValusPerPixel(sf, unused_list)
% take the ROI pixels and put then into one array

N = size(sf, 1);
M = size(sf, 2);

pt = [];

if(nargin<2)
    unused_list = [];
end

for n=1:N
    
    if(~isempty(unused_list))
        ul = unused_list{n};
    end
    
    for m=1:M                
        if(~isempty(unused_list))
            if(~isempty(ul))        
                pind = find(ul==m);
                if(~isempty(pind))
                    continue;
                end
            end
        end
        
        pt_curr = sf{n, m};        
        if(pt_curr(1)~=-1)
            pt = [pt; pt_curr(:)];
        end        
    end
end