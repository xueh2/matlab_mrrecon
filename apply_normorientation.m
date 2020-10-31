function [I, gt_h_norm] = apply_normorientation(I, Rotate90, Flip, gt_h);
% function [I, gt_h_norm] = apply_normorientation(I, Rotate90, Flip, gt_h);
% 

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  NIH NHLBI                          *
%     ***************************************
% 1 ---- 2  E1, phase_dir
% |      |
% |      |
% 4 ---- 3 
% RO, read_dir
    
if(nargin<4)
    gt_h = [];
end

Isize=size(I);
tmp=I(:,:,:);

gt_h_norm = gt_h;

% rotate then flip
if ndims(tmp)>2
    for i=1:size(tmp,3)
        tmp2(:,:,i)=rot90(tmp(:,:,i), Rotate90);
    end   
    tmp=tmp2;
else
    tmp=rot90(tmp,Rotate90);
end

RO = size(I, 1);
E1 = size(I, 2);

if(~isempty(gt_h_norm))
    if(Rotate90==1) % counter clockwise 90 degree

        [p1, p2, p3, p4] = gt_2_world(gt_h_norm, RO, E1);
        gt_h_norm = world_2_gt(p2, p3, p4, p1, gt_h_norm);

    %     read_dir = gt_h_norm.read_dir;
    %     phase_dir = gt_h_norm.phase_dir;
        fov = gt_h_norm.FOV;
    %     
    %     gt_h_norm.read_dir = phase_dir;
    %     gt_h_norm.phase_dir = -read_dir;
        gt_h_norm.FOV(1) = fov(2);
        gt_h_norm.FOV(2) = fov(1);
    end

    if(Rotate90==-1) % clockwise 90 degree
    %     read_dir = gt_h_norm.read_dir;
    %     phase_dir = gt_h_norm.phase_dir;
        fov = gt_h_norm.FOV;
    %     
    %     gt_h_norm.read_dir = -phase_dir;
    %     gt_h_norm.phase_dir = read_dir;
        gt_h_norm.FOV(1) = fov(2);
        gt_h_norm.FOV(2) = fov(1);

        [p1, p2, p3, p4] = gt_2_world(gt_h_norm, RO, E1);
        gt_h_norm = world_2_gt(p4, p1, p2, p3, gt_h_norm);
    end

    if(Rotate90==2) % counter clockwise 180 degree    
        [p1, p2, p3, p4] = gt_2_world(gt_h_norm, RO, E1);
        gt_h_norm = world_2_gt(p3, p4, p1, p2, gt_h_norm);    
    end
end

RO = size(tmp, 1);
E1 = size(tmp, 2);

switch Flip
    case 0 % do nothing
    case 1
        tmp=tmp(end:-1:1,:,:);
    case 2
        tmp=tmp(:,end:-1:1,:);
end

if(~isempty(gt_h_norm))
    if(Flip==1)
    %     read_dir = gt_h_norm.read_dir;
    %     phase_dir = gt_h_norm.phase_dir;
    %     
    %     gt_h_norm.phase_dir = -phase_dir;

        [p1, p2, p3, p4] = gt_2_world(gt_h_norm, RO, E1);
        gt_h_norm = world_2_gt(p4, p3, p2, p1, gt_h_norm);
    end

    if(Flip==2)
    %     read_dir = gt_h_norm.read_dir;
    %     phase_dir = gt_h_norm.phase_dir;
    %     
    %     gt_h_norm.read_dir = -read_dir;

        [p1, p2, p3, p4] = gt_2_world(gt_h_norm, RO, E1);
        gt_h_norm = world_2_gt(p2, p1, p4, p3, gt_h_norm);
    end
end

I = reshape(tmp,[size(tmp,1) size(tmp,2) Isize(3:end)]);

end

function [p1, p2, p3, p4] = gt_2_world(gt_h, RO, E1)

    % 1 ---- 2  E1, phase_dir
    % |      |
    % |      |
    % 4 ---- 3 
    % RO, read_dir
    
    c_ro = floor(RO/2);
    c_e1 = floor(E1/2);
    
    ps = max(gt_h.FOV)/max([RO, E1]);
    
    p1 = gt_h.PatientPosition + (-0.5-c_ro) * ps * gt_h.read_dir + (-0.5-c_e1) * ps * gt_h.phase_dir;
    p2 = gt_h.PatientPosition + (-0.5-c_ro) * ps * gt_h.read_dir + (E1-1+0.5-c_e1) * ps * gt_h.phase_dir;
    p3 = gt_h.PatientPosition + (RO-1+0.5-c_ro) * ps * gt_h.read_dir + (E1-1+0.5-c_e1) * ps * gt_h.phase_dir;
    p4 = gt_h.PatientPosition + (RO-1+0.5-c_ro) * ps * gt_h.read_dir + (-0.5-c_e1) * ps * gt_h.phase_dir;
end

function gt_h = world_2_gt(p1, p2, p3, p4, gt_h)

    % 1 ---- 2
    % |      |
    % |      |
    % 4 ---- 3 

    gt_h.read_dir = p4 - p1;
    gt_h.read_dir = gt_h.read_dir / norm(gt_h.read_dir);
    
    gt_h.phase_dir = p2 - p1;
    gt_h.phase_dir = gt_h.phase_dir / norm(gt_h.phase_dir);
    
    gt_h.PatientPosition = (p1+p2+p3+p4)/4;
end

   
