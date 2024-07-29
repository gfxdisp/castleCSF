% Plot CSF functions

if ~contains(path, '../')
    addpath(genpath('../'));
end

csf_model = CSF_castleCSF();

% 3D plot - as the function of spatial and temporal frequency
figure;

%%%%%%%%%%%%%%%%%%%% Inputs to the model %%%%%%%%%%%%%%%%%%%%%%%%%%

t_frequency = logspace( 0.1, log10(64), 30 );           % Temporal frequency in Hz
s_frequency = logspace( log10(0.5), log10(64), 30 );    % Spatial frequency in cycles per degree


orientation = 0;    % Orientation of grating in degrees
area = 1;           % Area of stimulus in visual sq. degrees
eccentricity = 0;   % Retinal eccentricity in degrees

luminance = 100;                    % Mean luminance of background in cd/m^2
xy_background = [0.3127, 0.3290];   % xy chromaticity coordinates of stimulus background

luminance_modulation = 200;         % Peak luminance of modulation; Same as background luminance in case of isoluminant stimuli
xy_modulation = [0.3127, 0.3290];   % Direction ( in xy chromaticity coordinates ) of stimulus chromatic modulation

% example chromaticity coordinates: 
%   white/black: [0.3127, 0.3290]
%   red/green:  [0.5027, 0.2322]
%   yellow/violet: [0.3313, 0.3738]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lms_background = xyz2lms2006(Yxy2XYZ([luminance, xy_background]));
lms_modulation_peak = xyz2lms2006(Yxy2XYZ([luminance_modulation, xy_modulation]));
lms_modulation_delta = lms_modulation_peak-lms_background;
lms_delta_norm = lms_modulation_delta./norm(lms_modulation_delta);


[ss, tt] = meshgrid( s_frequency, t_frequency );

csf_pars = struct( 's_frequency', ss(:), 't_frequency', tt(:),... 
    'orientation', orientation, 'area', area, 'eccentricity', eccentricity,...
    'lms_bkg', lms_background,...
    'lms_delta', lms_delta_norm);     

S = csf_model.sensitivity( csf_pars );        
S = reshape( S, size(ss) );

surf( s_frequency, t_frequency, S, 'FaceColor', 'interp', 'FaceLighting', 'phong' );
set( gca, 'XScale', 'log' );
set( gca, 'YScale', 'log' );
set( gca, 'ZScale', 'log' );
zlim( [1 1000] );
xlabel( 'Spatial frequency [cpd]')
ylabel( 'Temporal frequency [Hz]')
zlabel( 'Sensitivity')
title( 'castleCSF');

if ~contains(path, '../')
    rmpath(genpath('../'));
end