function plotGMMMIX_Parzon_1D(mix, data, plot_a)
% x is the intensity values that you want to draw, not the original intensity data

pp = data(find(data>0));
n = hist(pp, max(data));

n = n ./ max(n);
figure;
stem(n, 'Marker', 'none', 'Color', [100 100 100]./255);
a = plot_a;
pa = a * (mix.priorx)';
% pa = gmmprob_drawGMM(mix, mix.priorx, x);
x = [min(data(:)):max(data(:))]';

hold on
for i = 1:mix.ncentres
    
    if ( mod(i,8) == 0 )
        MarkerEdgeColor = 'b';
    end
    if ( mod(i,8) == 1 )
        MarkerEdgeColor = [0.21568627450980   0.07843137254902   0.41960784313725];
    end
    if ( mod(i,8) == 2 )
        MarkerEdgeColor = 'g';
    end
    if ( mod(i,8) == 3 )
        MarkerEdgeColor = 'm';
    end
    if ( mod(i,8) == 4 )
        MarkerEdgeColor = 'c';
    end
    if ( mod(i,8) == 5 )
        MarkerEdgeColor = [139 124 211]./255;
    end
    if ( mod(i,8) == 6 )
        MarkerEdgeColor = 'k';
    end
    if ( mod(i,8) == 7 )
        MarkerEdgeColor = [0.5 0.5 0.5];
    end
    
    plot(x, a(:,i), 'Color', MarkerEdgeColor, 'LineWidth', 1);
    
end

plot(x, pa./max(pa), 'Color', 'r', 'LineWidth', 2);
% plot(x, pa, 'Color', 'r', 'LineWidth', 2);

hold off
box on
return