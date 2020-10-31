function set_env_apply_ai(fid)
% set_env_apply_ai(fid)

GT_HOME = getenv('GADGETRON_HOME');
GADGETRON_SCRIPTS_FOLDER = getenv('GADGETRON_SCRIPTS_FOLDER');
GT_CMR_ML_UNITTEST_DIRECTORY = getenv('GT_CMR_ML_UNITTEST_DIRECTORY');

if(isunix())
    fprintf(fid, '%s\n', ['set GADGETRON_HOME=' GT_HOME]);
    fprintf(fid, '%s\n', ['set GT_CMR_ML_UNITTEST_DIRECTORY=' GT_CMR_ML_UNITTEST_DIRECTORY]);    
    fprintf(fid, '%s\n', ['cd ' GT_HOME '/share/gadgetron/python/cmr_ml']);
else
    fprintf(fid, '%s\n', 'chcp 65001');
    fprintf(fid, '%s\n', ['set GADGETRON_HOME=' GT_HOME]);
    fprintf(fid, '%s\n', ['set GT_CMR_ML_UNITTEST_DIRECTORY=' GT_CMR_ML_UNITTEST_DIRECTORY]);
    fprintf(fid, '%s\n', ['cd /D ' GT_HOME '/share/gadgetron/python/cmr_ml']);
end