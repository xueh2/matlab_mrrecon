% The accuracy is characterized by an error between the true signal and its
% estimate. The error is characterized by the bias and variance of
% estimates. This example illustrates an asymptotic analysis for the bias
% and variance for h and delta (see Chapter 5 Discrete LPA Accuracy).
%
% Dmitriy Paliy. Tampere University of Technology. TICSP. 06-06-2005
% dmitriy.paliy@tut.fi
clear all

delta = 0.0001;
sigma_g = 1;
x = 0:delta:5; % arguments
gamma = 1/2; % the parameter gamma is a ratio of the upper bound bias to 
             % the standard deviation calculated for the ideal scale h
             % (5.25)-(5.26), (5.43).

% constats computed analytically
Ag = 0.5*sigma_g^2;
Bg = 1/sqrt(4*pi);

sigma = Ag/sqrt(delta*Bg*gamma^2);

% the estimated signal
y = x.^2.*cos(3*x);

derivative_kernel=[1 -1]/1;
% the first derivative of the estimated signal
y1 = conv(y,derivative_kernel)/delta;
y1 = y1(2:end);
y1(1)=y1(2);
y1(end)=y1(end-1);
% the second derivative of the estimated signal
y2 = conv(y1,derivative_kernel)/delta;
y2 = y2(2:end);
y2(1:2)=y2(3);
y2(end-1:end)=y2(end-2);


version -release; % get matlab release
matlab_R=str2num(ans);

figure,

    plot(x, y, 'LineWidth', 2); hold on;

    plot(x, y1, 'r-', 'LineWidth', 2); hold on;

    plot(x, y2, 'g-', 'LineWidth', 2); hold on;
    
    title('Estimated function and two derivatives'),

    if matlab_R>=14,
        xlabel('$x$','interpreter','latex');
    else
        xlabel('\itx');
    end;
    grid on
    if matlab_R>=14,
        legend('y(x)','y^{(1)}(x)','y^{(2)}(x)',...
            'Location','NorthWest');
    else
        legend('y(x)','y^{(1)}(x)','y^{(2)}(x)',...
            2);
    end;
% -------------------------------------------------------


h_const = gamma^2*delta*sigma^2*Bg/(Ag^2);
L_m = abs(y2);

figure,

    % the ideal varying scale h(x) (5.41)
    h = (h_const./(L_m.^2)).^(1/5);
    plot(x, h, 'LineWidth', 2); hold on; YLim([0 4]),
    title('Ideal varying scale h^*(x)')
    if matlab_R>=14,
        xlabel('$x$','interpreter','latex');
        ylabel('$h^*(x)$','interpreter','latex');
    else
        xlabel('\itx');
        ylabel('\ith^*(\itx)');
    end;
    grid on

    

figure,

    % FOR THE IDEAL VARYING SCALE
    % the squared bias
    bias = h.^2.*L_m.*Ag;
    plot(x, bias.^2, 'r-', 'LineWidth', 2); hold on;

    % the variance
    var = (sigma^2*Bg*delta)./h;
    plot(x, var, 'g-', 'LineWidth', 2); hold on;

    % the upper bound of the risk l(x) (5.42)
    l_risk = bias.^2 + var;
    plot(x, l_risk, '-', 'LineWidth', 2); hold on;

    
    % the risk for the ideal INVARIANT scale
    L_m_2 = sum(L_m.^2)./size(L_m,2);
    h1 = (h_const/L_m_2).^(1/5);

    bias = h1.^2.*L_m.*Ag;
    var = (sigma^2*Bg*delta)./h1;
    l_risk1 = bias.^2 + var;
    plot(x, l_risk1, 'k-', 'LineWidth', 2); hold on;
    if matlab_R>=14,
        xlabel('$x$','interpreter','latex');
    else
        xlabel('\itx');
    end;
    grid on
    
    
    if matlab_R>=14,
        legend('Squared bias','Variance',...
            'Upper bound of risk for ideal varying h',...
            'Upper bound of risk for ideal invariant h',...
            'Location','NorthWest')
    else
        legend('Squared bias','Variance',...
            'Upper bound of risk for ideal varying h',...
            'Upper bound of risk for ideal invariant h',...
            2)
    end;