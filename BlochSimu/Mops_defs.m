function M = Mops_defs;
% function M = Mops_defs
% Export function handles of the operators below
% Going through this procedure allows you to access operators in
%  both standard MATLAB scripts and functions. It also allows for
%  many different simulations to use the same basis set of
%  operators, and visualization functions

% Creat a cell-list of available functions
fs={'Mx','My','Mz','Mxy','Mmag','Mphase','RotX','RotY','RotZ','Renorm',...
    'RFHard','RFSoft','RFwaveform','Precess','Relax','Evolve','rad','deg',...
	'IdealCrush','gamma','gammabar'};

% Convert each name into a function handle available from structure M
for i=1:length(fs),
	M.(fs{i}) = str2func(fs{i});
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = gammabar
g = 42.577e3;            % gamma in MHz/T or Hz/mT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = gamma
g = gammabar*2*pi;       % gamma in rad/s/mT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Magnetization Operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = Mx(M);
M = squeeze(M(1,:,:,:))';

function M = My(M);
M = squeeze(M(2,:,:,:))';

function M = Mz(M);
M = squeeze(M(3,:,:,:))';

function M = Mxy(M)
M = squeeze(sqrt(sum(M(1:2,:,:,:).^2, 1)));

function M = Mmag(M);
M = sqrt(sum(M(1:3,:,:,:).^2,1));

function ph = Mphase(M)
ph = atan2(M(2,:,:,:),M(1,:,:,:)); 

function M = Renorm(M)
% Function to normalize a M vector
M = M./M(4);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Rotation Operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
function Mop = RotX(c,s)
% x rotation operator
if nargin ==1, s = sin(c); c = cos(c); end;
Mop = [...
    1 , 0 , 0 , 0;...
    0 , c ,-s , 0;...
    0 , s , c , 0;...
    0 , 0 , 0 , 1];


%%
function Mop = RotY(c,s)
% Y rotation operator
if nargin ==1, s = sin(c); c = cos(c); end;
Mop = [...
    c , 0 , s , 0;...
    0 , 1 , 0 , 0;...
    -s , 0 , c , 0;...
    0 , 0 , 0 , 1];

%%
function Mop = RotZ(c,s);
% Z rotation operator
c1 = c;
if nargin ==1, s=sin(c); c = cos(c); end;
% disp(num2str([c1, s,c]));
Mop = [...
    c ,-s , 0 , 0;...
    s , c , 0 , 0;...
    0 , 0 , 1 , 0;...
    0 , 0 , 0 , 1];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Standard MR Operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mop = RFHard(alpha, phase)
if nargin==1, phase = zeros(size(alpha)); end;
Mop = RotZ(phase) * RotX(alpha) * RotZ(-phase);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mop = RFSoft(theta, phi, psi)
if nargin == 2, psi = zeros(size(theta)); end;
if nargin == 1, phi = zeros(size(theta)); end;
% Needs to be expanded to include all the RF points (which it doesn't do
% now)
Mop = RotZ(-phi) * RotY(-psi) * RotZ(theta) * RotY(psi) * RotZ(phi);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mop = Precess(theta)
Mop = RotZ(theta);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mop = Relax(time, T1, T2)
E1 = exp(-time/T1);
E2 = exp(-time/T2);
Mop = [...
    E2 , 0 , 0 , 0;...
    0 , E2 , 0 , 0;...
    0 , 0  ,E1 , 1-E1;...
    0 , 0  , 0 , 1];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mop = IdealCrush;
Mop = [...
    0 , 0 , 0 , 0;...
    0 , 0 , 0 , 0;...
    0 , 0  ,1 , 0;...
    0 , 0  ,0 , 1];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mop = Evolve(delta_B0, B1_amp, B1_phase, time, T1, T2, gamma)
B_effective = sqrt(B1_amp^2 + delta_B0^2);
MRelax   = Relax(time/2, T1,T2);
MPrecess = ...
    RotZ( B1_phase) * ...
    RotY(delta_B0/B_effective,  B1_amp/B_effective) *...
    RotZ(2*pi*gamma* B_effective * time) * ...
    RotY(delta_B0/B_effective, -B1_amp/B_effective) *...
    RotZ(-B1_phase);
Mop = MRelax * MPrecess * MRelax;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Support functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mop = MakeOp4x4(Mop)
% Function to turn a 3x3 operator into a 4x4 op
if size(Mop,1) ==3
    Mop = [Mop, [0;0;0]; [0 0 0 1]];
end;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alpha = rad(alpha)
alpha = alpha *pi/180;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alpha = deg(alpha)
alpha = alpha /pi*180;









