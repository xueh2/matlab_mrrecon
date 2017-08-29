
function PerformGadgetronRecon_Encrypt_GtPrep_config(src_dir, dst_dir)
% PerformGadgetronRecon_Encrypt_GtPrep_config(src_dir, dst_dir)
% PerformGadgetronRecon_Encrypt_GtPrep_config('D:\gtuser\mrprogs\install\share\gadgetron\config', 'D:\gtuser\mrprogs\gt_scanner_setup\VE11\config')
% PerformGadgetronRecon_Encrypt_GtPrep_config('D:\gtuser\mrprogs\install\share\gadgetron\config', 'D:\gtuser\mrprogs\gt_scanner_setup\VE11C\config')
% PerformGadgetronRecon_Encrypt_GtPrep_config('D:\gtuser\mrprogs\install\share\gadgetron\config', 'D:\gtuser\mrprogs\gt_scanner_setup\VE11B\config')

[names, num] = findFILE(src_dir, '*.xml');

key = 'GadgetronGoodForMRI';

delete(fullfile(dst_dir, '*.xml'));

for ii=1:num    
    disp(['Encrypting ' names{ii}])
    [path, name, ext] = fileparts(names{ii});    
    
    if( ~isempty(strfind(name, 'localhost')) )
        continue;
    end
    if( ~isempty(strfind(name, 'Generic_Cartesian')) )
        continue;
    end
    
    encrypted_config = Matlab_gt_config_xml_encryption(names{ii}, fullfile(dst_dir, [name '.xml']), key);
end