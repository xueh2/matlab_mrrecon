% finds radial thickening by iterating over each pixel in the 
% edge of epi mask


endo = retro_cine_endo_mask;
epi = retro_cine_epi_mask;
mask = epi - endo;
imsz = size(mask);
phase = 1;

PHS = imsz(3);
[Rmesh, Emesh] = meshgrid(1:imsz(2), 1:imsz(1));
for i = 1:(PHS)
    % find centroid of mask
    centroidR(i) = sum(sum(Rmesh(mask(:, :, i) == 1)))/sum(sum(mask(:, :, i) == 1));
    centroidE(i) = sum(sum(Emesh(mask(:, :, i) == 1)))/sum(sum(mask(:, :, i) == 1));
end
dists = sqrt((repmat(Rmesh, 1, 1, PHS) - permute(repmat(centroidR, imsz(1), 1, imsz(2)), [1, 3, 2])).^2 + ...
    (repmat(Emesh, 1, 1, PHS) - permute(repmat(centroidE, imsz(1), 1, imsz(2)), [1, 3, 2])).^2);

X = repmat(Rmesh, 1, 1, PHS) - permute(repmat(centroidR, imsz(1), 1, imsz(2)), [1, 3, 2]);
Y = -(repmat(Emesh, 1, 1, PHS) - permute(repmat(centroidE, imsz(1), 1, imsz(2)), [1, 3, 2]));
theta = atan(Y./(X+1e-8)) + pi*(X < 0) + pi*2*(X>=0).*(Y <0);
imagescn(theta)


X_in = Rmesh + 1*cos(theta);
X_out = Rmesh - 1*cos(theta);
Y_in = Emesh - 1*sin(theta);
Y_out = Emesh + 1*sin(theta);


edge_endo = zeros(imsz);
edge_epi = zeros(imsz);
for i = 1:PHS
    edge_endo(:, :, i) = abs(interp2(endo(:, :, i), X_out(:,:,i), Y_out(:,:,i)) - interp2(endo(:, :, i), X_in(:,:,i), Y_in(:,:,i), 'linear'));
    edge_epi(:, :, i) = abs(interp2(epi(:, :, i), X_out(:,:,i), Y_out(:,:,i)) - interp2(epi(:, :, i), X_in(:,:,i), Y_in(:,:,i), 'linear'));
end
edge_endo = (edge_endo > 0.8);
edge_epi = (edge_epi > 0.8);
% figure;imagescn(edge_endo, [], [], [], 3)
% figure;imagescn(edge_epi, [], [], [], 3)

ref_epi_edge = edge_epi(:, :, phase);
ref_endo_edge = edge_endo(:, :, phase);
upsample = 1;

mask = imresize(mask, upsample);
theta = imresize(theta, upsample);
imsz_upsampled = size(mask);
% dx_res = zeros(imsz_upsampled);
% dy_res = zeros(imsz_upsampled);
strain_res = zeros(imsz_upsampled);
samples = max(imsz_upsampled);
edge_endo = imresize(edge_endo, upsample);
edge_epi = imresize(edge_epi, upsample);
ref_phase = 1;


tic
for e = 1:imsz_upsampled(1)
    for r = 1:imsz_upsampled(2)
        if (edge_epi(e, r, ref_phase) == 1)
            theta_pt = theta(e, r, ref_phase);
            
            check_e = e - samples/8*sin(theta_pt);
            check_r = r + samples/8*cos(theta_pt);
%             if sum(abs([pi/4, 3*pi/4, 5*pi/4, 7*pi/4] - theta_pt ) < 0.1) == 1
%                 [check_r, check_e, theta_pt]
%             end
%           
%             check_e - centroidE(ref_phase);
%             check_r - centroidR(ref_phase);
            test_e_ref = linspace(centroidE(ref_phase)*upsample, check_e, samples);
            test_r_ref = linspace(centroidR(ref_phase)*upsample, check_r, samples);
            
            [~, epi_ind_ref] = max(interp2(single(edge_epi(:, :, ref_phase)),test_r_ref, test_e_ref, 'linear'));
            [~, endo_ind_ref] = max(interp2(single(edge_endo(:, :, ref_phase)),test_r_ref, test_e_ref, 'linear'));
            epi_e_ref = (test_e_ref(epi_ind_ref));
            epi_r_ref = (test_r_ref(epi_ind_ref));
            endo_e_ref = (test_e_ref(endo_ind_ref));
            endo_r_ref = (test_r_ref(endo_ind_ref));
            myo_dist_ref = sqrt((epi_e_ref-endo_e_ref)^2 + (epi_r_ref-endo_r_ref)^2) ;

            for p = 1:imsz_upsampled(3)
                check_e_new = centroidE(p) - centroidE(ref_phase) + check_e;
                check_r_new = centroidR(p) - centroidR(ref_phase) + check_r;
%                 check_e_new - centroidE(p);
%                 check_r_new - centroidR(p);
                test_e = linspace(centroidE(p)*upsample, check_e_new, samples);
                test_r = linspace(centroidR(p)*upsample, check_r_new, samples);
                [~, endo_ind] = max(interp2(single(edge_endo(:, :, p)), test_r, test_e, 'linear'));
                [~, epi_ind] = max(interp2(single(edge_epi(:, :, p)), test_r, test_e, 'linear'));
                endo_e = (test_e(endo_ind));
                endo_r = (test_r(endo_ind));
                epi_e = (test_e(epi_ind));
                epi_r = (test_r(epi_ind));
                myo_dist = sqrt((epi_e-endo_e)^2 + (epi_r-endo_r)^2) ;

                strain_res(e, r, p) = (myo_dist - myo_dist_ref)/myo_dist_ref;
            end
        end
    end
end
toc

imagescn(strain_res, [], [], [], 3)


%second try: per theta, find endo and epi to find compression
%tic
% for i = 1:numel(test_thetas)
%     theta_pt = test_thetas(i);
%     if (theta_pt >= 7*pi/4) || ((theta_pt >= 0) && (theta_pt < pi/4))
%         len_r = imsz_upsampled(2) - centroidR_mask*upsample;
%         len_e = -len_r*tan(theta_pt);
%     elseif (theta_pt >= pi/4) && (theta_pt < 3*pi/4)
%         len_e = -centroidE_mask*upsample;
%         len_r = abs(len_e)/tan(theta_pt);
%     elseif (theta_pt >= 3*pi/4) && (theta_pt < 5*pi/4)
%         len_r = -centroidR_mask*upsample;
%         len_e = -len_r*tan(theta_pt);
%     else
%         len_e = imsz_upsampled(1) - centroidE_mask*upsample;
%         len_r = -len_e/tan(theta_pt);
%     end
%     test_e = linspace(centroidE_mask*upsample, centroidE_mask*upsample + len_e, samples);
%     test_r = linspace(centroidR_mask*upsample, centroidR_mask*upsample + len_r, samples);
%     
%     [~, endo_ind_ref] = max(interp2(single(edge_endo(:, :, 1)),test_r, test_e, 'spline'));
%     [~, epi_ind_ref] = max(interp2(single(edge_epi(:, :, 1)), test_r, test_e, 'spline'));
%     endo_e_ref = round(test_e(endo_ind_ref));
%     endo_r_ref = round(test_r(endo_ind_ref));
%     epi_e_ref = round(test_e(epi_ind_ref));
%     epi_r_ref = round(test_r(epi_ind_ref));
%     myo_dist_ref = sqrt((epi_e_ref-endo_e_ref)^2 + (epi_r_ref-endo_r_ref)^2) ;
%     
%     for p = 1:imsz_upsampled(3)
%         
%         [~, endo_ind] = max(interp2(single(edge_endo(:, :, p)),test_r, test_e, 'spline'));
%         [~, epi_ind] = max(interp2(single(edge_epi(:, :, p)), test_r, test_e, 'spline'));
%         endo_e = round(test_e(endo_ind));
%         endo_r = round(test_r(endo_ind));
%         epi_e = round(test_e(epi_ind));
%         epi_r = round(test_r(epi_ind));
%         myo_dist = sqrt((epi_e-endo_e)^2 + (epi_r-endo_r)^2) ;
%         
%         strain_thetas(i, p) = (myo_dist - myo_dist_ref)/myo_dist_ref;
%     end
% end
% toc
% 
% 
% num_segments = 6;
% segments = linspace(0, 2*pi, num_segments+1);
% rad_strain = zeros(imsz_upsampled);
% seg_len = numel(test_thetas)/num_segments;
% for p = 1:PHS
%     for segs = 1:num_segments
%         seg_val = sum(strain_thetas(seg_len*(segs - 1)+1:seg_len*(segs), p))/seg_len;
%         temp = (theta(:, :) < segments(segs+1)) & (theta(:, :) > segments(segs));
%         rad_strain(:, :, p) = rad_strain(:, :, p) + temp.*mask(:, :, p).*seg_val;
%     end    
% end
% figure;imagescn(rad_strain, [-0.1, 0.6], [], [], 3)


% First try: per pixel in myocardium mask thickening calculation
% 
%
% tic
% for e = 1:imsz_upsampled(1)
%     for r = 1:imsz_upsampled(2)
%         for p = 1:1
%             if mask(e, r, p) == 1
%                 theta_pt = theta(e, r);
%                 if ((theta_pt >= 7*pi/4) && (theta_pt < 0)) || ((theta_pt >= 0) && (theta_pt < pi/4))
%                     len = imsz_upsampled(2) - centroidR_mask*upsample;
%                 elseif (theta_pt >= pi/4) && (theta_pt < 3*pi/4)
%                     len = centroidE_mask*upsample;
%                 elseif (theta_pt >= 3*pi/4) && (theta_pt < 5*pi/4)
%                     len = centroidR_mask*upsample;
%                 else
%                     len = imsz_upsampled(1) - centroidE_mask*upsample;
%                 end
%                 dist_e = -len*sin(theta_pt);
%                 dist_r = len*cos(theta_pt);
%                 test_e = linspace(centroidE_mask*upsample, centroidE_mask*upsample + dist_e, samples);
%                 test_r = linspace(centroidR_mask*upsample, centroidR_mask*upsample + dist_r, samples);
%                 [~, endo_ind] = max(interp2(single(warped_endo(:, :, p)),test_r, test_e, 'spline'));
%                 [~, epi_ind] = max(interp2(single(warped_epi(:, :, p)), test_r, test_e, 'spline'));
%                 [endoarr_i, endoarr_j] = ind2sub(imsz_upsampled, endo_ind);
%                 endo_e = round(test_e(endoarr_i));
%                 endo_r = round(test_r(endoarr_j));
%                 [epiarr_i, epiarr_j] = ind2sub(imsz_upsampled, epi_ind);
%                 epi_e = round(test_e(epiarr_i));
%                 epi_r = round(test_r(epiarr_j));
%                 dist_endo = sqrt((e-endo_e)^2 + (r-endo_r)^2) ;
%                 dist_epi = sqrt((e-epi_e)^2 + (r-epi_r)^2) ;
%                 
%                 dx_res(e, r, p) = (dist_epi*interp2(dx_endo(:, :, p), r/upsample, e/upsample) + dist_endo*interp2(dx_epi(:, :, p), r/upsample, e/upsample))/(dist_epi + dist_endo);
%                 dy_res(e, r, p) = (dist_epi*interp2(dy_endo(:, :, p), r/upsample, e/upsample) + dist_endo*interp2(dy_epi(:, :, p), r/upsample, e/upsample))/(dist_epi + dist_endo);
%             end
%         end
%     end
% end
% toc
% 
% figure;imagescn([dx_res(:, :, 1), dy_res(:, :, 1)])
% [rad, circ] = numericalStrain(mask(:, :, 1), dx_res(:, :, 1), dy_res(:, :, 1));
% figure; imagescn([rad, circ])
% 
% centroidE_endo = zeros(PHS, 1);
% centroidR_endo = zeros(PHS, 1);
% centroidE_epi = zeros(PHS, 1);
% centroidR_epi = zeros(PHS, 1);
% for i = 1:(PHS)
%     % find centroid of mask
%     centroidR_endo(i) = sum(sum(Rmesh(endo(:, :, i) == 1)))/sum(sum(endo(:, :, i) == 1));
%     centroidE_endo(i) = sum(sum(Emesh(endo(:, :, i) == 1)))/sum(sum(endo(:, :, i) == 1));
%     centroidR_epi(i) = sum(sum(Rmesh(epi(:, :, i) == 1)))/sum(sum(epi(:, :, i) == 1));
%     centroidE_epi(i) = sum(sum(Emesh(epi(:, :, i) == 1)))/sum(sum(epi(:, :, i) == 1));
% end
% 
% 
% 
% 
% dx = zeros(endo);
% dy = zeros(endo);
% ref_endo = endo(:, :, phase);
% ref_epi = epi(:, :, :, phase);
% 
% 
% 
% 
% return