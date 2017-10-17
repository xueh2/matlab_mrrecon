function [out] = roi_statistics(I, ROI_info_table);

% function [out]=roi_statistics(I, ROI_info_table);
%
% function to take stack of images, I(row,col,frame), and an roimatfile
% created by imagescn.m containing 1 or more ROIs and output the ROI data
% and ROI mean in matlab structure (out) for all frames. User specified
% which ROI number, and if the ROI is drawn on a multi-axis imagescn figure
% then the use specifies which image_no (axis)

%      ***************************************
%      *  Peter Kellman  (kellman@nih.gov)   *
%      *  Laboratory for Cardiac Energetics  *
%      *  NIH NHLBI                          *
%      ***************************************


BW=zeros(size(I(:,:,1)));

BW=roipoly(I(:,:,1),ROI_info_table.ROI_x_coordinates,ROI_info_table.ROI_y_coordinates);

index=find(BW >0);

x=linspace(-50, 50,20);
for j=1:size(I,4)
    for i=1:size(I,3)
        tmp=I(:,:,i,j);
        data(i,j,:)=tmp(index);
        m(i,j)=mean(tmp(index));
        s(i,j)=std(tmp(index));
        snr(i,j) = m(i,j)/s(i,j);
    %     h(i,:)=hist(tmp(index),x);
    end
end

% for i=1:size(I,3)
%     tmp=I(:,:,i);
%     tmp(find(tmp)==0)=nan;
%     m(i)=nanmean(tmp(index));
%     s(i)=nanstd(tmp(index));
% end

out.data=data;
out.m=m;
out.s=s;
out.snr=snr;
out.BW=BW;
% out.h=h;

return
