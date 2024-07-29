function dkl = lms2dkl_d65( lms )
% Convert from LMS color space to DKL color space assuming adaptation to
% D65 white.

% lms_gray - the LMS coordinates of the white point
lms_gray = [0.739876529525622   0.320136241543338   0.020793708751515];

mc1 = lms_gray(1)/lms_gray(2);
mc2 = (lms_gray(1)+lms_gray(2))/lms_gray(3);

M_lms_dkl = [ 1  1 0;
               1 -mc1 0;
              -1 -1 mc2 ];

dkl = cm_colorspace_transform( lms, M_lms_dkl );          
          
end
