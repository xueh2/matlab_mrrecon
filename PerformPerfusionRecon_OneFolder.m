function [t_linear, t_nl] = PerformPerfusionRecon_OneFolder(UT_Dir, home, gt_host, runLinear, runNL, linear_res_dir, nl_res_dir, deleteh5, xml, xmlNL, copy_debug)
% PerformPerfusionRecon_OneFolder('E', 'DualBolus\NewData\KAROLINSKA', 'palau', runLinear, runNL, 'grappa_flow_res_BTEX20_TwoCompWithShifts', 'slep_flow_res_BTEX20_TwoCompWithShifts', 0)
% PerformPerfusionRecon_OneFolder('F', 'DualBolus\FFR', 'palau', runLinear, runNL, 'grappa_flow_res_BTEX20_TwoCompWithShifts', 'slep_flow_res_BTEX20_TwoCompWithShifts', 0)

if(nargin<6)
    linear_res_dir = 'grappa_flow_res_BTEX20_TwoCompWithShifts';
end

if(nargin<7)
    linear_res_dir = 'slep_flow_res_BTEX20_TwoCompWithShifts';
end

if(nargin<8)
    deleteh5 = 0;
end

if(nargin<9)
    xml = 'GTPrep_2DT_Perfusion_AIF_TwoEchoes_Interleaved_R2_QuantitativeFlow_Mapping.xml';
end

if(nargin<10)
    xmlNL = 'GTPrep_2DT_Perfusion_AIF_TwoEchoes_Interleaved_R2_QuantitativeFlow_Mapping_NonLinear.xml';
end

if(nargin<11)
    copy_debug = 0;
end
%% 
set_UT_Dir(UT_Dir)
UTDir = getenv('GTPLUS_UT_DIR')

home_all = fullfile(UTDir, home)

exDir = {'ICERecon', 'grappa', 'TXMapping', 'moco', 'mocoSyn', 'mocoPS' , 'mocoPSSyn', 'seg', 'magFitting', 'slep_res', 'slep_cloud_res', 'grappa_res', 'ICE', 'DebugOutput', 'PhantomT1maps', 'slep_flow_res', 'grappa_flow_res', 'slep_cloud_flow_res', 'slep_cloud_flow_res2', 'T2Map', 'test', 'Dotarem_r1_r2_reduced', 'grappa_flow_res_BTEX20_TwoCompWithShifts', 'slep_flow_res_BTEX20_TwoCompWithShifts'};
[subdir, num] = FindAllEndDirectoryExclusive(home_all, exDir)

startI = 1;

t_linear = [];
t_nl = [];

for ii=startI:num  

    if ( ~isempty(strfind( lower(subdir{ii}), 'mini')) | ~isempty(strfind(subdir{ii}, 'MINI')) | ~isempty(strfind(subdir{ii}, 'Mini')) )
        continue;
    end

    if ( ~isempty(strfind( lower(subdir{ii}), 't1map')) )
        continue;
    end

    [tl, tnl] = PerformPerfusionRecon_OneCase(UT_Dir, home, subdir{ii}, gt_host, runLinear, runNL, linear_res_dir, nl_res_dir, deleteh5, xml, xmlNL, copy_debug); 
    
    t_linear = [t_linear; {subdir{ii} tl}];
    t_nl = [t_nl; {subdir{ii} tnl}];
end
