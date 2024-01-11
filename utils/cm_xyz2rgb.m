function rgb = cm_xyz2rgb(xyz, rgb_space)
% CM_RGB2XYZ  Convert (selected) RGB to XYZ color space
%
% rgb = cm_xyz2rgb(xyz, rgb_space)
%
% "xyz" could be an image [width x height x 3] or [rows x 3] matrix of colour
%     vectors, or [width x height x 3 x frames] colour video

if ~exist( 'rgb_space', 'var' )
    rgb_space = 'rec709';
end

switch rgb_space
    case 'Adobe'
        M_xyz2rgb = [
            2.04148   -0.969258  0.0134455
            -0.564977 1.87599    -0.118373
            -0.344713 0.0415557  1.01527 ];
    case 'NTSC'
        M_xyz2rgb = [ 1.9099961 -0.5324542 -0.2882091
             -0.9846663  1.9991710 -0.0283082
             0.0583056  -0.1183781  0.8975535]; 
    case { 'sRGB', 'rec709' }
        M_xyz2rgb = [ 3.2406 -1.5372 -0.4986;
             -0.9689  1.8758  0.0415;
              0.0557 -0.2040  1.0570];
    case 'rec2020'
        % Source: https://colour.readthedocs.io/en/v0.3.7/colour.models.dataset.rec_2020.html
        M_xyz2rgb = inv( [0.636953507, 0.144619185, 0.168855854;
               0.262698339, 0.678008766, 0.0592928953;
               4.99407097e-17, 0.0280731358, 1.06082723] );
    otherwise
        error( 'Unknown RGB colour space' );
end

rgb = cm_colorspace_transform( xyz, M_xyz2rgb );

end
    