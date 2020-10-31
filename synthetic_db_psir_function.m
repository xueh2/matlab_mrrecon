function [PSIR,DB_PSIR] = synthetic_db_psir_function(t1map,t2map,T1blood, T1myo, T2blood, T2myo);
% function [PSIR,DB_PSIR] = synthetic_db_psir_function(t1map,t2map,T1blood, T1myo, T2blood, T2myo);

% MI Case
% % case 1
% idir = 'Z:\ReconResults\Barts\20200306\T1SR_Mapping_SASHA_HC_T1T2_41837_1082381512_1082381521_216_20200306-115104_dicom';
% % DB_LGE_MOCO_AVE_OnTheFly_41837_1082381512_1082381521_217_20200306-115149_figure
% 
% 
% % case 2
% % (pre-contrast) idir = 'Z:\ReconResults\Barts\20200306\T1SR_Mapping_SASHA_HC_T1T2_42110_1142134903_1142134912_201_20200306-114402_dicom';
% % idir = 'Z:\ReconResults\Barts\20200306\T1SR_Mapping_SASHA_HC_T1T2_42110_1142134903_1142134912_220_20200306-121539_dicom';
% % LGE_MOCO_AVE_OnTheFly_42110_1142134903_1142134912_217_20200306-121335_figure
% 
% t1 = dicomread([idir,filesep,'T1_SLC0_CON0_PHS0_REP0_SET0_AVE0_1.dcm']);
% t2 = dicomread([idir,filesep,'T2_SLC0_CON0_PHS0_REP0_SET0_AVE0_1.dcm']);
% % figure; imagescn(t1,[0 2000]); siemens_map
% % figure; imagescn(t2,[0 120]); siemens_map

t1 = double(t1map);
t2 = double(t2map);

% optimize synthetic DB PSIR.

% % case 2
% T1blood = 520;
% T2blood = 157;
% T1myo = 680;
% T2myo = 43;
% % case 1
% T1blood = 420;
% T2blood = 157;
% T1myo = 650;
% T2myo = 39;

% null myocardium
te = [10:5:70];

clear td1_* td2_*
for j = 1:length(te)
    TE = te(j);

    TD2 = 0; % fixed
    % binary search to find TD1 that nulls myo
    TD1 = 400; % range 0 to 800
    delta = 200;
    for i = 1:10;
        M = db_psir_function(T1myo, T2myo, TD1, TD2, TE);
        if M > 0
            TD1 = TD1 - delta;
        else
            TD1 = TD1 + delta;
        end
        delta = delta/2;
    end
    td1_myo(j) = TD1;

    TD1 = 0; % fixed
    % binary search to find TD2 that nulls myo
    TD2 = 400; % range 0 to 800
    delta = 200;
    for i = 1:10;
        M = db_psir_function(T1myo, T2myo, TD1, TD2, TE);
        if M > 0
            TD2 = TD2 - delta;
        else
            TD2 = TD2 + delta;
        end
        delta = delta/2;
    end
    td2_myo(j) = TD2;


    TD2 = 0; % fixed
    % binary search to find TD1 that nulls blood
    TD1 = 400; % range 0 to 800
    delta = 200;
    for i = 1:10;
        M = db_psir_function(T1blood, T2blood, TD1, TD2, TE);
        if M > 0
            TD1 = TD1 - delta;
        else
            TD1 = TD1 + delta;
        end
        delta = delta/2;
    end
    td1_blood(j) = TD1;

    TD1 = 0; % fixed
    % binary search to find TD2 that nulls blood
    TD2 = 400; % range 0 to 800
    delta = 200;
    for i = 1:10;
        M = db_psir_function(T1blood, T2blood, TD1, TD2, TE);
        if M > 0
            TD2 = TD2 - delta;
        else
            TD2 = TD2 + delta;
        end
        delta = delta/2;
    end
    td2_blood(j) = TD2;
    
    % calculate line for myo null: TD1 = m_myo*TD2 + b_myo;
    b_myo = td1_myo(j); % point p2;    
    m_myo = -b_myo/td2_myo(j); 
    % calculate TD1 and TD2 that null myo
    % and that Mblood is desired level
    Mblood_desired = -0.15; % make this a config variable
    t = [1:800];
    for k = 1:length(t);
        TD2 = t(k);
        TD1 = m_myo*TD2 + b_myo;
        Mblood(k) = db_psir_function(T1blood, T2blood, TD1, TD2, TE);
    end
    ind = min(find(Mblood < Mblood_desired));
    if ~isempty(ind); calcTD2(j) = t(ind); else; calcTD2(j) = 0; end
    calcTD1(j) = m_myo*calcTD2(j) + b_myo;
        
%     % check for crossing of myo & blood where myo & blood both null      
%     % calculate line for myo null: TD1 = m_myo*TD2 + b_myo;
%     b_myo = td1_myo(j); % point p2;    
%     m_myo = -b_myo/td2_myo(j);    
%     b_blood = td1_blood(j); % point p2;    
%     m_blood = -b_blood/td2_blood(j);     
%     delta_td1_target1 = -50;
%     calcTD2(j) = (delta_td1_target1 - (b_myo - b_blood))/(m_myo - m_blood);
%     calcTD1(j) = m_myo * calcTD2(j) + b_myo;
%     calcTD1b(j) = m_blood * calcTD2(j) + b_blood;
%     dTD1(j) = calcTD1(j) - calcTD1b(j);
end % TE loop

for j = 1:length(calcTD1)
    if (calcTD1(j) > 0 & calcTD2(j) > 0)
        TD1 = calcTD1(j);
        TD2 = calcTD2(j);
        TE = te(j);
        break;
    end
end

% disp(['TD1 =  ',num2str(TD1)])
% disp(['TD2 =  ',num2str(TD2)])
% disp(['TE =  ',num2str(TE)])

DB_PSIR = db_psir_function(t1, t2, TD1, TD2, TE);
% figure; imagescn(cat(3,abs(DB_PSIR),DB_PSIR),[],[],8) % DB PSIR

% compute bright blood PSIR
    % binary search to find TD1 (TI) that nulls myo
    TD1 = 400; % range 0 to 800
    delta = 200;
    for i = 1:10;
        M = db_psir_function(T1myo, T2myo, TD1, 0, 0);
        if M > 0
            TD1 = TD1 - delta;
        else
            TD1 = TD1 + delta;
        end
        delta = delta/2;
    end
    td1_myo(j) = TD1;

PSIR = db_psir_function(t1, t2, TD1, 0, 0);
% figure; imagescn(cat(3,PSIR,DB_PSIR),[],[],10) % DB PSIR


return












