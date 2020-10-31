function [ro_s, ro_e, e1_s, e1_e] = detect_heart_region_LGE_roi(data, NN_RO, NN_E1)
% detect heart region
% [ro_s, ro_e, e1_s, e1_e] = detect_heart_region_LGE_roi(data, NN_RO, NN_E1)

RO = size(data,1);
E1 = size(data,2);
SLC = size(data,3);
   
if(SLC>1)
    figure; imagescn(data(:,:, round(SLC/2)+2), [4096-800 4096+800], [], [8]);
else
    figure; imagescn(data, [4096-800 4096+800], [], [8]);
end

[x,y] = getpts;

c_ro = y;
c_e1 = x;

ro_s = c_ro - NN_RO/2;
if(ro_s<1)
    ro_s = 1;
end
ro_e = ro_s + NN_RO -1;
if(ro_e>RO)
    ro_e = RO;
    ro_s = ro_e - NN_RO + 1;
end

e1_s = c_e1 - NN_E1/2;
if(e1_s<1)
    e1_s = 1;
end
e1_e = e1_s + NN_E1 -1;
if(e1_e>E1)
    e1_e = E1;
    e1_s = e1_e - NN_E1 + 1;
end
