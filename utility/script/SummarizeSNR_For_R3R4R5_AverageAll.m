
R3Data = [];
R4Data = [];
R5Data = [];

R3SNR = [];
R4SNR = [];
R5SNR = [];

R3scoreMean = [];
R4scoreMean = [];
R5scoreMean = [];

R3scoreMax = [];
R4scoreMax = [];
R5scoreMax = [];

R3score10pMax = [];
R4score10pMax = [];
R5score10pMax = [];

R3scorePValue = [];
R4scorePValue = [];
R5scorePValue = [];

R3SNRDelta = [];
R4SNRDelta = [];
R5SNRDelta = [];

SNRRecord_Problem = [];
PathWrong = [];

for d=1:num
    
    currDir = fullfile(home, subdir{d})

    % only process stress
    if ( (isempty(strfind(currDir, 'STRESS'))==1) & (isempty(strfind(currDir, 'stress'))==1) )
        isStress = 0;
    else
        isStress = 1;
    end
    
    [R, isStress] = findDataInfoR3R4R5(subdir{d});
    if ( ~isStress | R==3 )
        continue;
    end
    
    cd(currDir);
    load(fullfile(currDir, 'totalSlices.mat'))
     
    for s=1:numOfSlice
        
        sDir = ['Slice' num2str(s)];
        
        cd(fullfile(currDir, sDir))
        
        load(SNRMat)
        
        if ( ~strcmp(currSNR{1}, subdir{d}) | ~strcmp(currSNR{2}, sDir) )
            disp('path wrong');
            PathWrong = [PathWrong; {subdir{d}} {sDir}];
        end
        
        SNR = currSNR{5};
        R_Ori = currSNR{3};

%         if ( isempty(strfind(currDir, '_R3_'))~=1 )        
%             R = 3;
%         end
% 
%         if ( isempty(strfind(currDir, 'R3_'))~=1 )
%             R = 3;
%         end
%         
%         if ( isempty(strfind(currDir, '_3_'))~=1 )        
%             R = 3;
%         end
% 
%         if ( isempty(strfind(currDir, '_R4_'))~=1 )        
%             R = 4;
%         end
% 
%         if ( isempty(strfind(currDir, '_r4_'))~=1 )        
%             R = 4;
%         end
% 
%         if ( isempty(strfind(currDir, '_4_'))~=1 )        
%             R = 4;
%         end
% 
%         if ( isempty(strfind(currDir, '_R5_'))~=1 )
%             R = 5;
%         end
% 
%         if ( isempty(strfind(currDir, '_r5_'))~=1 )
%             R = 5;
%         end
% 
%         if ( isempty(strfind(currDir, '_5_'))~=1 )        
%             R = 5;
%         end
%         R        

        R = -1;
        if ( isempty(strfind(currDir, 'R3_'))~=1 )        
            R = 3;
        end

        if ( isempty(strfind(currDir, '_R3_'))~=1 )        
            R = 3;
        end

        if ( isempty(strfind(currDir, '_3_'))~=1 )        
            R = 3;
        end

        if ( isempty(strfind(currDir, '_R4_'))~=1 )        
            R = 4;
        end

        if ( isempty(strfind(currDir, '_4_'))~=1 )        
            R = 4;
        end

        if ( isempty(strfind(currDir, '_R5_'))~=1 )
            R = 5;
        end

        if ( isempty(strfind(currDir, '_r5_'))~=1 )
            R = 5;
        end

        if ( isempty(strfind(currDir, '_5_'))~=1 )        
            R = 5;
        end
        R 
    
        if ( R ~= R_Ori )
            subdir{d}
        end
        
        % R = R_Ori;
        
        % isStress = currSNR{4};
        if ( ~isStress )
            continue;
        end

%         if ( R==3 )        
%             if ( isempty(R3SNR) )            
%                 R3SNR = [SNR(1:3)];
%                 R3SNRDelta = [SNR(1:3)./SNR(1)];
%                 if ( R3SNRDelta(3)<1 )
%                    d 
%                    SNRRecord_Problem = [SNRRecord_Problem; {subdir{d}} {sDir} ];
%                 end
%                 R3scoreMean = [SNR([4 6 9])];
%                 R3scoreMax = [SNR([5 7 10])];
%                 R3scorePValue = [SNR([8 11])];
%                 R3Data = {[subdir{d} '/' sDir]};
%             else
%                 R3SNR = [R3SNR; SNR(1:3)];
%                 R3SNRDelta = [R3SNRDelta; SNR(1:3)./SNR(1)];
%                 if ( R3SNRDelta(end,3)<1 )
%                     d
%                     SNRRecord_Problem = [SNRRecord_Problem; {subdir{d}} {sDir}];
%                 end
%                 R3scoreMean = [R3scoreMean; SNR([4 6 9])];
%                 R3scoreMax = [R3scoreMax; SNR([5 7 10])];
%                 R3scorePValue = [R3scorePValue; SNR([8 11])];
%                 R3Data = [R3Data {[subdir{d} '/' sDir]}];
%             end        
%         end
        
        asScore = currSNR{6};            
        numOfFrame = size(asScore, 1);
        asTGRAPPA = asScore(:,1);
        asTSPIRIT = asScore(:,2);
        asSPIRIT = asScore(:,3);

        asTGRAPPA10pMax = sort(asTGRAPPA);
        asTGRAPPA10pMax = mean(asTGRAPPA10pMax(floor(0.9*numOfFrame):end));
   
        asTSPIRIT10pMax = sort(asTSPIRIT);
        asTSPIRIT10pMax = mean(asTSPIRIT10pMax(floor(0.9*numOfFrame):end));

        asSPIRIT10pMax = sort(asSPIRIT);
        asSPIRIT10pMax = mean(asSPIRIT10pMax(floor(0.9*numOfFrame):end));
        
        if ( R==4 )                        
            if ( isempty(R4SNR) )            
                R4SNR = [SNR(1:3)];
                R4SNRDelta = [SNR(1:3)./SNR(1)];
                if ( R4SNRDelta(3)<1 )
                   d 
                   SNRRecord_Problem = [SNRRecord_Problem; {subdir{d}} {sDir} ];
                end
                R4scoreMean = [SNR([4 6 9])];
                R4scoreMax = [SNR([5 7 10])];
                R4scorePValue = [SNR([8 11])];
                R4Data = {[subdir{d} '/' sDir]};
                R4score10pMax = [asTGRAPPA10pMax asTSPIRIT10pMax asSPIRIT10pMax];                                
            else
                R4SNR = [R4SNR; SNR(1:3)];
                R4SNRDelta = [R4SNRDelta; SNR(1:3)./SNR(1)];
                if ( R4SNRDelta(end,3)<1 )
                    d
                    SNRRecord_Problem = [SNRRecord_Problem; {subdir{d}} {sDir}];
                end
                R4scoreMean = [R4scoreMean; SNR([4 6 9])];
                R4scoreMax = [R4scoreMax; SNR([5 7 10])];
                R4scorePValue = [R4scorePValue; SNR([8 11])];
                R4Data = [R4Data {[subdir{d} '/' sDir]}];
                R4score10pMax = [R4score10pMax; asTGRAPPA10pMax asTSPIRIT10pMax asSPIRIT10pMax];
            end
        end
    
        if ( R==5 )        
            if ( isempty(R5SNR) )            
                R5SNR = [SNR(1:3)];
                R5SNRDelta = [SNR(1:3)./SNR(1)];
                if ( R5SNRDelta(3)<1 )
                   d 
                   SNRRecord_Problem = [SNRRecord_Problem; {subdir{d}} {sDir} ];
                end
                R5scoreMean = [SNR([4 6 9])];
                R5scoreMax = [SNR([5 7 10])];
                R5scorePValue = [SNR([8 11])];
                R5Data = {[subdir{d} '/' sDir]};
                R5score10pMax = [asTGRAPPA10pMax asTSPIRIT10pMax asSPIRIT10pMax]; 
            else
                R5SNR = [R5SNR; SNR(1:3)];
                R5SNRDelta = [R5SNRDelta; SNR(1:3)./SNR(1)];
                if ( R5SNRDelta(end,3)<1 )
                    d
                    SNRRecord_Problem = [SNRRecord_Problem; {subdir{d}} {sDir}];
                end
                R5scoreMean = [R5scoreMean; SNR([4 6 9])];
                R5scoreMax = [R5scoreMax; SNR([5 7 10])];
                R5scorePValue = [R5scorePValue; SNR([8 11])];
                R5Data = [R5Data {[subdir{d} '/' sDir]}];
                R5score10pMax = [R5score10pMax; asTGRAPPA10pMax asTSPIRIT10pMax asSPIRIT10pMax];
            end        
        end
    end
end

SNRRecord_Problem
PathWrong
