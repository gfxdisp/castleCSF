function XYZ = lms2006_2xyz( LMS )
% Transform from CIE 1931 XYZ trichormatic colour values to CIE 2006 LMS cone responses on an LED LCD display. 
%
% XYZ = lms2006_2xyz( LMS )
%
% LMS can be an image (h x w x 3) or a colour vectors (n x 3).
%
% Note that because CIE XYZ 1931 and CIE 2006 LMS are based on different
% colour matching functions, this transformation depends on the colour
% spectra of a display. This transformation was derived for the spectra of LED LCD. 
%
% The CIE 2006 LMS cone responses can be found at http://www.cvrl.org/.

M_lms2006_xyz = [ 
   2.629129278399650  -3.780202391780134  10.294956387893450;
   0.865649062438827   1.215555811642301  -0.984175688105352;
  -0.008886561474676   0.081612628990755  51.371024830897888 ];
  
XYZ = cm_colorspace_transform( LMS, M_lms2006_xyz );

end
