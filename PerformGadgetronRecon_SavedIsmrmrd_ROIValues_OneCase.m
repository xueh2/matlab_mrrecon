
function res = PerformGadgetronRecon_SavedIsmrmrd_ROIValues_OneCase(s1, s2, s3, a)
% res = PerformGadgetronRecon_SavedIsmrmrd_ROIValues_OneCase(s1, s2, s3, stress)
% given the ROIs s1/s2/s3, get the flow and other values
% res has [flow, Ki, E, Visf, Vp, PS, SD, Ki_MF, Ki_Fermi, Ki_TwoCompExp, Ki_BTEX]

if(isfield(a, 'flow_stress'))
    %% stress
    
    stress = a;
    
    [f1, f2, f3] = get_roi_values(stress.flow_stress, s1, s2, s3);
    res.flow = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(stress.E_stress, s1, s2, s3);
    res.E = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(stress.Visf_stress, s1, s2, s3);
    res.Visf = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(stress.Vp_stress, s1, s2, s3);
    res.Vp = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(stress.PS_stress, s1, s2, s3);
    res.PS = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(stress.SDMap_stress, s1, s2, s3);
    res.SD = [f1.m f2.m f3.m];

    num_maps = size(stress.Ki_stress, 3);
    
    [f1, f2, f3] = get_roi_values( squeeze(stress.Ki_stress(:,:,1,:)), s1, s2, s3);
    res.Ki_MF = [f1.m f2.m f3.m];

    if(num_maps==4)
        [f1, f2, f3] = get_roi_values( squeeze(stress.Ki_stress(:,:,2,:)), s1, s2, s3);
        res.Ki_Fermi = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_roi_values( squeeze(stress.Ki_stress(:,:,3,:)), s1, s2, s3);
        res.Ki_TwoCompExp = [f1.m f2.m f3.m];
    else
        res.Ki_Fermi = [-1 -1 -1];
        res.Ki_TwoCompExp = [-1 -1 -1];
    end
    
    [f1, f2, f3] = get_roi_values( squeeze(stress.Ki_stress(:,:,end,:)), s1, s2, s3);
    res.Ki_BTEX = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(stress.Delay_stress, s1, s2, s3);
    res.delay = [f1.m f2.m f3.m];

    res.flow_i = [-1 -1 -1];
    res.E_i = [-1 -1 -1];
    res.Visf_i = [-1 -1 -1];
    res.Vp_i = [-1 -1 -1];
    res.PS_i = [-1 -1 -1];
    res.SD_i = [-1 -1 -1];
    res.Ki_MF_i = [-1 -1 -1];
    res.Ki_Fermi_i = [-1 -1 -1];
    res.Ki_TwoCompExp_i = [-1 -1 -1];
    res.Ki_BTEX_i = [-1 -1 -1];
    
    %% if having second roi
    if( (~isempty(s1) & numel(s1.ROI_info_table)==2) ... 
        | (~isempty(s2) & numel(s2.ROI_info_table)==2) ...
        | (~isempty(s3) & numel(s3.ROI_info_table)==2) )

        [f1, f2, f3] = get_2nd_roi_values(stress.flow_stress, s1, s2, s3);
        res.flow_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(stress.E_stress, s1, s2, s3);
        res.E_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(stress.Visf_stress, s1, s2, s3);
        res.Visf_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(stress.Vp_stress, s1, s2, s3);
        res.Vp_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(stress.PS_stress, s1, s2, s3);
        res.PS_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(stress.SDMap_stress, s1, s2, s3);
        res.SD_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(stress.Ki_stress(:,:,1,:)), s1, s2, s3);
        res.Ki_MF_i = [f1.m f2.m f3.m];

        if(num_maps==4)
            [f1, f2, f3] = get_2nd_roi_values( squeeze(stress.Ki_stress(:,:,2,:)), s1, s2, s3);
            res.Ki_Fermi_i = [f1.m f2.m f3.m];

            [f1, f2, f3] = get_2nd_roi_values( squeeze(stress.Ki_stress(:,:,3,:)), s1, s2, s3);
            res.Ki_TwoCompExp_i = [f1.m f2.m f3.m];
        else
            res.Ki_Fermi_i = [-1 -1 -1];
            res.Ki_TwoCompExp_i = [-1 -1 -1];
        end
        
        [f1, f2, f3] = get_2nd_roi_values( squeeze(stress.Ki_stress(:,:,end,:)), s1, s2, s3);
        res.Ki_BTEX_i = [f1.m f2.m f3.m];
        
        [f1, f2, f3] = get_2nd_roi_values(stress.Delay_stress, s1, s2, s3);
        res.delay_i = [f1.m f2.m f3.m];
    end
else
    %% rest
    
    rest = a;
    
    [f1, f2, f3] = get_roi_values(rest.flow_rest, s1, s2, s3);
    res.flow = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(rest.E_rest, s1, s2, s3);
    res.E = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(rest.Visf_rest, s1, s2, s3);
    res.Visf = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(rest.Vp_rest, s1, s2, s3);
    res.Vp = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(rest.PS_rest, s1, s2, s3);
    res.PS = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(rest.SDMap_rest, s1, s2, s3);
    res.SD = [f1.m f2.m f3.m];

    num_maps = size(rest.Ki_rest, 3);

    [f1, f2, f3] = get_roi_values( squeeze(rest.Ki_rest(:,:,1,:)), s1, s2, s3);
    res.Ki_MF = [f1.m f2.m f3.m];

    if(num_maps==4)
        [f1, f2, f3] = get_roi_values( squeeze(rest.Ki_rest(:,:,2,:)), s1, s2, s3);
        res.Ki_Fermi = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_roi_values( squeeze(rest.Ki_rest(:,:,3,:)), s1, s2, s3);
        res.Ki_TwoCompExp = [f1.m f2.m f3.m];
    else
        res.Ki_Fermi = [-1 -1 -1];
        res.Ki_TwoCompExp = [-1 -1 -1];
    end
    
    [f1, f2, f3] = get_roi_values( squeeze(rest.Ki_rest(:,:,end,:)), s1, s2, s3);
    res.Ki_BTEX = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(rest.Delay_rest, s1, s2, s3);
    res.delay = [f1.m f2.m f3.m];
    
    res.flow_i = [-1 -1 -1];
    res.E_i = [-1 -1 -1];
    res.Visf_i = [-1 -1 -1];
    res.Vp_i = [-1 -1 -1];
    res.PS_i = [-1 -1 -1];
    res.SD_i = [-1 -1 -1];
    res.Ki_MF_i = [-1 -1 -1];
    res.Ki_Fermi_i = [-1 -1 -1];
    res.Ki_TwoCompExp_i = [-1 -1 -1];
    res.Ki_BTEX_i = [-1 -1 -1];
        
    %% if having second roi
    if( (~isempty(s1) & numel(s1.ROI_info_table)==2) ... 
        | (~isempty(s2) & numel(s2.ROI_info_table)==2) ...
        | (~isempty(s3) & numel(s3.ROI_info_table)==2) )

        [f1, f2, f3] = get_2nd_roi_values(rest.flow_rest, s1, s2, s3);
        res.flow_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(rest.E_rest, s1, s2, s3);
        res.E_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(rest.Visf_rest, s1, s2, s3);
        res.Visf_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(rest.Vp_rest, s1, s2, s3);
        res.Vp_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(rest.PS_rest, s1, s2, s3);
        res.PS_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(rest.SDMap_rest, s1, s2, s3);
        res.SD_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(rest.Ki_rest(:,:,1,:)), s1, s2, s3);
        res.Ki_MF_i = [f1.m f2.m f3.m];

        if(num_maps==4)
            [f1, f2, f3] = get_2nd_roi_values( squeeze(rest.Ki_rest(:,:,2,:)), s1, s2, s3);
            res.Ki_Fermi_i = [f1.m f2.m f3.m];

            [f1, f2, f3] = get_2nd_roi_values( squeeze(rest.Ki_rest(:,:,3,:)), s1, s2, s3);
            res.Ki_TwoCompExp_i = [f1.m f2.m f3.m];
        else
            res.Ki_Fermi_i = [-1 -1 -1];
            res.Ki_TwoCompExp_i = [-1 -1 -1];
        end
        
        [f1, f2, f3] = get_2nd_roi_values( squeeze(rest.Ki_rest(:,:,end,:)), s1, s2, s3);
        res.Ki_BTEX_i = [f1.m f2.m f3.m];
    end
end

end

function [f1, f2, f3] = get_roi_values(a, s1, s2, s3)
    if(~isempty(s1))
        f1 = roi_statistics(a(:,:,1), s1.ROI_info_table(1,1));
    else
        f1.m = -1;
    end
    
    if(~isempty(s2))
        f2 = roi_statistics(a(:,:,2), s2.ROI_info_table(1,1));
    else
        f2.m = -1;
    end
    
    if(~isempty(s3))
        f3 = roi_statistics(a(:,:,3), s3.ROI_info_table(1,1));
    else
        f3.m = -1;
    end
end

function [f1, f2, f3] = get_2nd_roi_values(a, s1, s2, s3)

    if(~isempty(s1) & numel(s1.ROI_info_table)==2)
        f1 = roi_statistics(a(:,:,1), s1.ROI_info_table(1,1));
    else
        f1.m = -1;
    end
    
    if(~isempty(s2) & numel(s2.ROI_info_table)==2)
        f2 = roi_statistics(a(:,:,2), s2.ROI_info_table(1,1));
    else
        f2.m = -1;
    end
    
    if(~isempty(s3) & numel(s3.ROI_info_table)==2)
        f3 = roi_statistics(a(:,:,3), s3.ROI_info_table(1,1));
    else
        f3.m = -1;
    end
end


