% Plot CSF functions

addpath('../matlab/');
addpath('../utils/');

csf_model = CSF_castleCSF();

% 2D plot - as the function of temporal frequency
figure;

%%%%%%%%%%%%%%%%%%%% Inputs to the model %%%%%%%%%%%%%%%%%%%%%%%%%%
s_frequency = 0.5;    % Spatial frequency in cycles per degree
orientation = 0;    % Orientation of grating in degrees
area = 1;           % Area of stimulus in visual sq. degrees
eccentricity = 0;   % Retinal eccentricity in degrees
t_frequency = linspace( 0, 60 )'; %Hz, must be a column vector

luminance = 100;                    % Mean luminance of background in cd/m^2
xy_background = [0.3127, 0.3290];   % xy chromaticity coordinates of stimulus background

luminance_modulation = 100;         % Peak luminance of modulation; Same as background luminance in case of isoluminant stimuli
xy_modulation = [0.5027, 0.2322];   % Direction ( in xy chromaticity coordinates ) of stimulus chromatic modulation

% example chromaticity coordinates: 
%   white/black: [0.3127, 0.3290]
%   red/green:  [0.5027, 0.2322]
%   yellow/violet: [0.3313, 0.3738]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lms_background = xyz2lms2006(Yxy2XYZ([luminance, xy_background]));
lms_modulation_peak = xyz2lms2006(Yxy2XYZ([luminance_modulation, xy_modulation]));
lms_modulation_delta = lms_modulation_peak-lms_background;
lms_delta_norm = lms_modulation_delta./norm(lms_modulation_delta);

csf_pars = struct( 's_frequency',s_frequency, 't_frequency', t_frequency,... 
    'orientation', orientation, 'area', area, 'eccentricity', eccentricity,...
    'lms_bkg', lms_background,...
    'lms_delta', lms_delta_norm);     

S = csf_model.sensitivity( csf_pars );        

plot( t_frequency, S );
set( gca, 'YScale', 'log' );
xlabel( 'Temporal frequency [Hz]' );
ylabel( 'Sensitivity' );


rmpath('../matlab/');
rmpath('../utils/');
