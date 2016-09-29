clear;clc;

%% add path
addpath(genpath('FISTA'));
addpath(genpath('operators'));

%% set up global variables
global sizeImage;
% it contains five elements  [1d,2d,c,t]

global coilSensitivity;
% it contains four elements  [1d,2d,c]

global samplingLocations;
% it contains four elements [1d,2d,t]

%% load data and CSm

cd D:\gtuser\mrprogs\install\gadgetron\DebugOutput

kspace = readGTPlusExportData('incomingKSpace');
kspace = squeeze(kspace);
norm(kspace(:))

cd D:\gtuser\mrprogs\gadgetron\toolboxes\gtplus\ut\result
coilmap = readGTPlusExportData('coilMap');
norm(coilmap(:))

cd D:\software\example

% header = CreateGtImageHeader(kspace);
% kspace = single(kspace);
% Matlab_gt_write_analyze(real(kspace), header, 'kspace_REAL');
% Matlab_gt_write_analyze(imag(kspace), header, 'kspace_IMAG');
% norm(kspace(:))
% 
% meanKspace = sum(kspace(:,:,:,1:15), 4);
% norm(meanKspace(:))
% plotKSpace(meanKspace);
% 
% meanIm = ifft2c(meanKspace);
% 
% [combinedIm, senMap_SPIRIT, eigD] = Spirit2DForSensitivity(meanIm, [5 5], 0.01);
% 
% y = kspace(:,:,:,1);
% RO = size(kspace, 1);
% E1 = size(kspace, 2);
% CHA = size(kspace, 3);
% senMap = zeros(RO, E1, CHA, 2);
% senMap(:,:,:,1) = coilmap;
% senMap(:,:,:,2) = 0;
% 
% im1 = SensitivityCoilCombination(ifft2c(y), senMap(:,:,:,1));
% im2 = SensitivityCoilCombination(ifft2c(y), senMap(:,:,:,2));
% im = cat(3, im1, im2);
% [res, RESVEC] = cgSENSEWithoutFOV(y, senMap, 100, 1e-6, im);
% [res, RESVEC] = cgSENSEWithoutFOV_NotLSQR(y, senMap, 100, 1e-6, im);
% [res, RESVEC] = cgSENSE(y, senMap, 100, 1e-6, im);
% 
% [ImCombined, Im, eigD, senMap, x0] = senseWithoutFOV(y, meanKspace, 100, 1e-6, senMap);
% 
% header = CreateGtImageHeader(meanIm);
% meanIm = single(meanIm);
% Matlab_gt_write_analyze(real(meanIm), header, 'ref_REAL');
% Matlab_gt_write_analyze(imag(meanIm), header, 'ref_IMAG');
% 
% meanIm = ifft2c(meanIm);
% plotComplexImage(meanIm);
% 
% senMap = CoilSensitivityMapEstimation(meanIm, 'Jun', [7 7], [5 5], 95, [7 7]);
% 
% senMap2 = CoilSensitivityMapEstimation(meanIm, 'Souheil', [7 7], [5 5], 95, [7 7]);
% 
% imagescn(abs(senMap));
% imagescn(abs(senMap2));
% 
% opts=[];
% opts.choice=1; % use opts.tau
% opts.percentage=96;
% kSize=[6,6];
% 
% opts.thd=0.98;
% 
% ACS = fft2c(meanIm);
% tic;
% [senMap, eigD, tau]=ED_eigen_2D_parallel_fov(ACS, kSize, size(meanIm), opts);
% toc;
% 
% imagescn(eigD);
% 
% im1 = SensitivityCoilCombination(meanIm, senMap(:,:,:,1));
% im2 = SensitivityCoilCombination(meanIm, senMap(:,:,:,2));
% figure; imagescn(abs(im1));
% figure; imagescn(abs(im2));
% 
% imagescn(abs(im1+im2));
% imagescn(abs(senMap(:,:,:,1)))

for id=[18, 26, 35, 43]
    
    for mask_num= 1:2
        
        fname='TPAT_05_07_2012\';
        
        
        for csm_approach=3:3
            switch csm_approach
                case 1
                    load(strcat(fname,'MID_',num2str(id),'_mask_',num2str(mask_num), '_csm_Siemens.mat'));
                    csm=csm_Siemens;
                case 2
                    load(strcat(fname,'MID_',num2str(id),'_mask_',num2str(mask_num), '_csm_b1Map.mat'));
                    csm=b1map;
                case 3
                    load(strcat(fname,'MID_',num2str(id),'_mask_',num2str(mask_num), '_csm_eigen_vector.mat'));
                    csm=sensitivityMap;
            end
            % load csm, size: 1d x 2d x c
            
            
            %% set up the input size, csm, k-space etc.
            sizeImage=size(kspace);
            sizeReconImage=[sizeImage(1:2), sizeImage(4)];
            coilSensitivity=[];
            samplingLocations=[];
            vCSimage=[];
            
            % set csm
            for i=1:sizeImage(3)
                coilSensitivity(:,:,i)=fftshift( csm(:,:,i) );
            end
            % clear csm;
            
            % set sampling locations
            sLocations=squeeze(kspace(:,:,1,:)~=0);
            for j=1:sizeImage(4)    %t
                samplingLocations(:,:,j)=fftshift(sLocations(:,:,j));
            end
            
            % set k-space data
            for i=1:sizeImage(3)      %c
                vCSimage_temp=[];
                for j=1:sizeImage(4)    %t
                    
                    akspace=fftshift(kspace(:,:,i,j));
                    
                    temp=akspace;
                    temp=temp( samplingLocations(:,:,j)==1 ) ;
                    vCSimage_temp=[vCSimage_temp; temp];
                end
                vCSimage(:,i)=vCSimage_temp;
            end
            
            
            WW=redundantHaarOperator(1)*redundantHaarOperator(2)*redundantHaarOperator(3);
            W = @(x) reshape(WW * reshape(x, [sizeImage(1), sizeImage(2), sizeImage(4)]), prod([2*sizeImage(1), 2*sizeImage(2), 2*sizeImage(4)]),1);
            WT= @(x) reshape(WW'* reshape(x, [2*sizeImage(1), 2*sizeImage(2), 2*sizeImage(4)]), prod([sizeImage(1), sizeImage(2), sizeImage(4)]),1);
            
            
            opts.n=prod(sizeReconImage);
            opts.maxIter=25;
            opts.init=0;
            lambda=1e-3;
            opts.rFlag=1;
            opts.maxFlag=0;
            opts.quickWarmStart=0;
            opts.scale=5;
            
            opts.pFlag=1; % 1: use the one for wavelet
            % 0: use the one for L1
            opts.W=W;
            opts.WT=WT;
            
            %     for lambda=[1e-4, 2e4, 5e-4, 1e-3]
            for scale= 5 %[5, 10, 20]
                
                opts.scale=scale;
                
                tic;
                [x, funVal,  ValueL]=fista_slep_2D_csm_wavelet(@cs_2D_csm_cartesian, @csCt_2D_csm_cartesian,...
                    vCSimage,...
                    lambda, opts);
                toc;
                
                X=reshape(x,sizeReconImage);
                
                XX=[];
                for i=1:sizeReconImage(3)
                    XX(:,:,i)=((fftshift(X(:,:,i))));
                end
                
                switch csm_approach
                    case 1
                        str=strcat(fnameStore,'MID_',num2str(id),'_mask_',num2str(mask_num),'_csm_Siemens','_scale_',num2str(opts.scale),'_reg_',num2str(lambda));
                    case 2
                        str=strcat(fnameStore,'MID_',num2str(id),'_mask_',num2str(mask_num),'_csm_b1Map','_scale_',num2str(opts.scale),'_reg_',num2str(lambda));
                    case 3
                        str=strcat(fnameStore,'MID_',num2str(id),'_mask_',num2str(mask_num),'_csm_eigen_vector','_scale_',num2str(opts.scale),'_reg_',num2str(lambda),'_40');
                        % for some reason, 40 is used in the one of May 15,
                        % 2012.
                        % here _40 does not have any meaning
                end
                
                
                save( strcat(str,'.mat'), 'XX');
                
            end
            %     end
        end
    end
end