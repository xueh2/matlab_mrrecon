
function timeUsed = BuildGadgetronMarsInstallers(target_dir, filename, remove_gtprep_xml)
%% build and mars installer
% remove_gtprep_xml = 0
% timeUsed = BuildGadgetronMarsInstallers(target_dir, filename, remove_gtprep_xml)
% timeUsed = BuildGadgetronMarsInstallers('\\hl-share\RawMRI\Lab-Kellman\Share\Installers\20220208_gt4px', 'gadgetron-20220208-1829-a541c4d6-gtprep-42cd70b1', remove_gtprep_xml)

if(nargin<3)
    remove_gtprep_xml = 0;
end

tic;

mkdir(target_dir);
       
if(remove_gtprep_xml)
    PerformGadgetronRecon_Encrypt_GtPrep_config('D:\gtuser\mrprogs\install\share\gadgetron\config', 'D:\gadgetron_scanner_setup\VE11\config');
    PerformGadgetronRecon_Encrypt_GtPrep_config('D:\gtuser\mrprogs\install\share\gadgetron\config', 'D:\gadgetron_scanner_setup\VE11C\config');    
    PerformGadgetronRecon_Encrypt_GtPrep_config('D:\gtuser\mrprogs\install\share\gadgetron\config', 'D:\gadgetron_scanner_setup\VE11E\config');    

    command = [getenv('GADGETRON_SCRIPTS_FOLDER') '\compile_gadgetron_package_MARS_All_Versions.bat ' filename ' 4 2 0'];
else
    delete('D:\gtuser\mrprogs\gt_scanner_setup\VE11\config\*.xml')
    delete('D:\gtuser\mrprogs\gt_scanner_setup\VE11C\config\*.xml')
    delete('D:\gtuser\mrprogs\gt_scanner_setup\VE11E\config\*.xml')

    command = [getenv('GADGETRON_SCRIPTS_FOLDER') '\compile_gadgetron_package_Sites_All_Versions.bat ' filename ' 4 2 0'];
end

command
dos(command, '-echo');

disp(['target_dir = ' target_dir]);
disp(['filename = ' filename]);
CopyGadgetronInstallersAll(target_dir, filename, '4.2.0');

timeUsed = toc;
