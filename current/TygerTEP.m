closeall
clear all

host = 'localhost'
port = '9002'

% ---------------------------------

data_dir = '/data/raw_data/tyger/rtcine'

findAndMoveMeasDat(data_dir)
[names, num] = FindSubDirs(data_dir)

UTCases = set_up_UT_cases_RTCine;
UTCases{1, 1} = data_dir;

i=1;
UTCases{1, 2} = names{i};
UTCases{1, 4} = 'GT_RTCine_AI.xml';
UTCases{1, 5} = 'res_GT_RTCine_AI';

UTCases{1, 4} = 'GT_RTCine_AI_NormOrientation_OR.xml';
UTCases{1, 5} = 'res_GT_RTCine_AI_NormOrientation_OR';

names{i}

performUTValidation(UTCases(1,:), 0, 0, host, port, 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')

data = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 1);
res = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 2101);

size(data)

h = figure('Name', names{i}); imagescn(cat(4, data, res), [0 2*abs(mean(data(:)))], [1 2], [18], 3);
saveas(h, fullfile('/data/raw_data/tyger', [names{i} '.fig']));

cd /data/raw_data/tyger/rtcine/RT_Cine_LIN_00000_237143837_237143846_525_00000000-000000/res_GT_RTCine_AI_NormOrientation_OR/


% ---------------------------------

data_dir = '/data/raw_data/tyger/retrocine'

findAndMoveMeasDat(data_dir)
[names, num] = FindSubDirs(data_dir)

UTCases = set_up_UT_cases_RTCine;
UTCases{1, 1} = data_dir;

i=1;
UTCases{1, 2} = names{i}
UTCases{1, 4} = 'default_measurement_dependencies.xml'
performUTValidation(UTCases(1,:), 0, 0, host, port, 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')
names{i}

i=2;
UTCases{1, 2} = names{i};
UTCases{1, 4} = 'GT_RetroCine_AI.xml';
UTCases{1, 5} = 'res_GT_RetroCine_AI';
names{i}

performUTValidation(UTCases(1,:), 0, 0, host, port, 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')

data = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 1);
res = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 2101);

size(data)

h = figure('Name', names{i}); imagescn(cat(4, data, res), [0 4*abs(mean(data(:)))], [1 2], [18], 4);
saveas(h, fullfile('/data/raw_data/tyger', [names{i} '.fig']));


% ---------------------------------

data_dir = '/data/raw_data/tyger/spine'

findAndMoveMeasDat(data_dir)
[names, num] = FindSubDirs(data_dir)

UTCases = set_up_UT_cases_3D;
UTCases{1, 1} = data_dir;

i=1;
UTCases{1, 2} = names{i};
UTCases{1, 4} = 'GT_2D_AI.xml';
UTCases{1, 5} = 'res_GT_2D_AI';
names{i}

performUTValidation(UTCases(1,:), 0, 0, host, port, 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')

data = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 1);
res = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 2101);

size(data)

h = figure('Name', names{i}); imagescn(cat(4, data, res), [0 4*abs(mean(data(:)))], [2 2], [12], 3);
saveas(h, fullfile('/data/raw_data/tyger', [names{i} '.fig']));

% ---------------------------------

data_dir = '/data/raw_data/tyger/neuro'

findAndMoveMeasDat(data_dir)
[names, num] = FindSubDirs(data_dir)

UTCases = set_up_UT_cases_3D;
UTCases{1, 1} = data_dir;

i=1;
UTCases{1, 2} = names{i};
UTCases{1, 4} = 'GT_3D_AI.xml';
UTCases{1, 5} = 'res_GT_3D_AI';
names{i}

performUTValidation(UTCases(1,:), 0, 0, host, port, 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')

data = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 1);
res = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 2101);

data = analyze75read(fullfile(data_dir, names{i}, UTCases{1, 5}, 'GT_3D_AI.xml_SLC0_CON0_PHS0_REP0_SET0_AVE0_1_1.hdr'));
res = analyze75read(fullfile(data_dir, names{i}, UTCases{1, 5}, 'GT_3D_AI.xml_SLC0_CON0_PHS0_REP0_SET0_AVE0_1_2101.hdr'));

size(data)

h = figure('Name', names{i}); imagescn(cat(4, data, res), [0 4*abs(mean(data(:)))], [1 2], [12], 3);
saveas(h, fullfile('/data/raw_data/tyger', [names{i} '.fig']));

