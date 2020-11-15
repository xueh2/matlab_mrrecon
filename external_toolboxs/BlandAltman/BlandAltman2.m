function [means,diffs,meanDiff,CR,linFit] = BlandAltman2(var1, var2, flag, symbol, plot_precentage, markersize)
% function [means,diffs,meanDiff,CR,linFit] = BlandAltman2(var1, var2, flag, symbol, plot_precentage, markersize)
 
    %%%Plots a Bland-Altman Plot
    %%%INPUTS:
    %%% var1 and var2 - vectors of the measurements
    %%%flag - how much you want to plot
        %%% 0 = no plot
        %%% 1 = just the data
        %%% 2 = data and the difference and CR lines
        %%% 3 = above and a linear fit
    %%%
    %%%OUTPUTS:
    %%% means = the means of the data
    %%% diffs = the raw differences
    %%% meanDiff = the mean difference
    %%% CR = the 2SD confidence limits
    %%% linfit = the paramters for the linear fit
    
    
    if (nargin<1)
       %%%Use test data
       var1=[512,430,520,428,500,600,364,380,658,445,432,626,260,477,259,350,451];%,...
       var2=[525,415,508,444,500,625,460,390,642,432,420,605,227,467,268,370,443];
       flag = 3;
    end
    
    if nargin==2
        flag = 0;
    end
    
    if nargin<5
        plot_precentage = 0;
    end
    
    if nargin<6
        markersize = 12;
    end
    
%     var1 = var1(:);
%     var2 = var2(:);
    
    means = mean([var1;var2]);
    diffs = var1-var2;
    
    meanDiff = mean(diffs);
    sdDiff = std(diffs);
    CR = [meanDiff + 1.96 * sdDiff, meanDiff - 1.96 * sdDiff]; %%95% confidence range
    
    CR_percentage = 100* [meanDiff + 1.96 * sdDiff, meanDiff - 1.96 * sdDiff] / meanDiff; 
    
    linFit = polyfit(means,diffs,1); %%%work out the linear fit coefficients
    
    %%%plot results unless flag is 0
    if flag ~= 0
        if(plot_precentage)
            plot(means,100*diffs/means,symbol, 'MarkerSize', markersize)
            hold on
            if flag > 1
                plot([0  1.1* max(means)], [CR_percentage(1) CR_percentage(1)],'k--', 'LineWidth', 2.0, 'Color', 'r'); %%%plot the upper CR
                plot([0  1.1* max(means)], [CR_percentage(2) CR_percentage(2)],'k--', 'LineWidth', 2.0, 'Color', 'r'); %%%plot the upper CR
                plot([0  1.1* max(means)], [0 0],'k-', 'LineWidth', 1.0); %%%plot the upper CR
                plot([0  1.1* max(means)], (CR_percentage(1)+CR_percentage(2))/2*[1 1],'k:', 'LineWidth', 2.0, 'Color', 'b'); %%%plot the upper CR                     
                set(gca, 'FontSize', 16);
                set(gca, 'LineWidth', 1.0);
            end
            title('Bland-Altman Plot (%)','fontsize',14)
        else
            plot(means,diffs,symbol, 'MarkerSize', markersize)
            hold on
            if flag > 1
                plot([0  1.1* max(means)], [CR(1) CR(1)],'k--', 'LineWidth', 2.0, 'Color', 'r'); %%%plot the upper CR
                plot([0  1.1* max(means)], [CR(2) CR(2)],'k--', 'LineWidth', 2.0, 'Color', 'r'); %%%plot the upper CR
                plot([0  1.1* max(means)], [0 0],'k-', 'LineWidth', 1.0); %%%plot the upper CR
                plot([0  1.1* max(means)], (CR(1)+CR(2))/2*[1 1],'k:', 'LineWidth', 2.0, 'Color', 'b'); %%%plot the upper CR                     
                set(gca, 'FontSize', 16);
                set(gca, 'LineWidth', 1.0);
            end
            title('Bland-Altman Plot','fontsize',14)
        end
        
        if flag > 2
            plot(means, means.*linFit(1)+linFit(2),'k--'); %%%plot the linear fit
        end        
    end
    
return