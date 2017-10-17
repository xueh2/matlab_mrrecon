

load('d:/work/imagescn/New_MV_Tool_20070807/test_ims.mat');
imagescn(ims, [ -12 165], [], 6, 3);

% Create a default "object" structure
ObjStruct = struct(...
    'name', [], ...
    'type', [], ...  % Should be Point,Line or Patch
    'xdata',[], ...  % coordinates for points and lines, vertices for patches
    'ydata',[], ...
    'color',[], ...  % color to use for drawing, string or [1x3] vector
    'linestyle',[],... % type of lines - for lines only
    'markerstyle', []...
    );  % type of markers to use for points and lines only

% make up fake data    
[XX,YY ] = meshgrid(20:5:size(ims,1)/2-20,20:5:size(ims,2)/2-20);
th = linspace(0,2*pi,100);
circx = cos(th);
circy = sin(th);

c= [size(ims,1)/2, size(ims,2)/2];

base_trig = [ 0 0; 1 0; 1,1; 0 1];
th = linspace(0,2*pi,size(ims,3));
trig_cmap = [ones(size(ims,3),1), linspace(1,0,size(ims,3))', linspace(1,0,size(ims,3))'];

% Name will appear on the pulldown menu

% for each c-phase
for i=1:size(ims,3) 
    ObjStruct(1,i).name = 'Grid Points 1';
    ObjStruct(1,i).type = 'Points';
    ObjStruct(1,i).xdata = YY+i;
    ObjStruct(1,i).ydata = XX+i;
    ObjStruct(1,i).color = trig_cmap(i,:);
    ObjStruct(1,i).marker = '.';

    
    ObjStruct(2,i).name = 'ROI1 ';
    ObjStruct(2,i).type = 'Line';
    ObjStruct(2,i).xdata = circx*i*2 + c(2)+ c(2)/2;
    ObjStruct(2,i).ydata = circy*i*2 + c(1);
    ObjStruct(2,i).color = 'b';
    %ObjStruct(2,i).marker = '.';
    
    ObjStruct(3,i).name ='Square Patch';
    ObjStruct(3,i).type ='Patch';
    ObjStruct(3,i).xdata = 10*(base_trig(:,1)*cos(th(i))+base_trig(:,2)*sin(th(i)))+40;
    ObjStruct(3,i).ydata = 10*(base_trig(:,2)*cos(th(i))+base_trig(:,1)*sin(th(i)))+100;
    ObjStruct(3,i).color = shiftdim(trig_cmap(i,:),-1);
    ObjStruct(3,i).facealpha = (i-1)/(size(ims,3)-1);
    
    ObjStruct(4,i).name = 'Grid Points 2';
    ObjStruct(4,i).type = 'Points';
    ObjStruct(4,i).xdata = YY+i+40;
    ObjStruct(4,i).ydata = XX+i+40;
    ObjStruct(4,i).color = trig_cmap(i,:);
    ObjStruct(4,i).marker = '.';
end

h_ax = findobj(gcf, 'type', 'axes');
for i = 1:length(h_ax)
    setappdata(h_ax(i),'Objects',ObjStruct);
end;

% Now MV_tool will use this data for display during movie