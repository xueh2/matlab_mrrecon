R3SNR = [];
R4SNR = [];
R5SNR = [];

R3scoreMean = [];
R4scoreMean = [];
R5scoreMean = [];

R3scoreMax = [];
R4scoreMax = [];
R5scoreMax = [];

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
    if ( isempty(strfind(currDir, 'STRESS'))==1 )
        isStress = 0;
    else
        isStress = 1;
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
        R = currSNR{3};

        isStress = currSNR{4};
        if ( ~isStress )
            continue;
        end

        if ( R==3 )        
            if ( isempty(R3SNR) )            
                R3SNR = [SNR(1:2)];
                R3SNRDelta = [SNR(1:2)./SNR(1)];
                if ( R3SNRDelta(2)<1 )
                   d 
                   SNRRecord_Problem = [SNRRecord_Problem; {subdir{d}} {sDir} ];
                end
                R3scoreMean = [SNR([3 5])];
                R3scoreMax = [SNR([4 6])];
                R3scorePValue = [SNR(7)];
            else
                R3SNR = [R3SNR; SNR(1:2)];
                R3SNRDelta = [R3SNRDelta; SNR(1:2)./SNR(1)];
                if ( R3SNRDelta(end,2)<1 )
                    d
                    SNRRecord_Problem = [SNRRecord_Problem; {subdir{d}} {sDir}];
                end
                R3scoreMean = [R3scoreMean; SNR([3 5])];
                R3scoreMax = [R3scoreMax; SNR([4 6])];
                R3scorePValue = [R3scorePValue; SNR(7)];
            end        
        end

        if ( R==4 )
            if ( isempty(R4SNR) )            
                R4SNR = [SNR(1:2)];
                R4SNRDelta = [SNR(1:2)./SNR(1)];
                if ( R4SNRDelta(2)<1.1 )
                    d
                    SNRRecord_Problem = [SNRRecord_Problem; {subdir{d}} {sDir}];
                end
                R4scoreMean = [SNR([3 5])];
                R4scoreMax = [SNR([4 6])];
                R4scorePValue = [SNR(7)];
            else
                R4SNR = [R4SNR; SNR(1:2)];
                R4SNRDelta = [R4SNRDelta; SNR(1:2)./SNR(1)];
                if ( R4SNRDelta(end,2)<1.1 )
                    d
                    SNRRecord_Problem = [SNRRecord_Problem; {subdir{d}} {sDir}];
                end            
                R4scoreMean = [R4scoreMean; SNR([3 5])];
                R4scoreMax = [R4scoreMax; SNR([4 6])];
                R4scorePValue = [R4scorePValue; SNR(7)];
            end
        end
    end
    
    if ( R==5 )        
        if ( isempty(R5SNR) )            
            R5SNR = [SNR(1:2)];
            R5SNRDelta = [SNR(1:2)./SNR(1)];
            if ( R5SNRDelta(2)<1.1 )
                d
                SNRRecord_Problem = [SNRRecord_Problem; {subdir{d}} {sDir}];
            end            
            R5scoreMean = [SNR([3 5])];
            R5scoreMax = [SNR([4 6])];
            R5scorePValue = [SNR(7)];
        else
            R5SNR = [R5SNR; SNR(1:2)];
            R5SNRDelta = [R5SNRDelta; SNR(1:2)./SNR(1)];
            if ( R5SNRDelta(end,2)<1.1 )
                d
                SNRRecord_Problem = [SNRRecord_Problem; {subdir{d}} {sDir}];
            end            
            R5scoreMean = [R5scoreMean; SNR([3 5])];
            R5scoreMax = [R5scoreMax; SNR([4 6])];
            R5scorePValue = [R5scorePValue; SNR(7)];
        end         
    end
end

SNRRecord_Problem
PathWrong
