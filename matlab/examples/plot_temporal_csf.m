% Plot CSF functions

if ~exist('CSF_castleCSF', 'file')
    addpath('../../matlab/');
    addpath('../../utils/');
end

csf_model = CSF_castleCSF();

% 2D plot - as the function of temporal frequency
figure;

%%%%%%%%%%%%%%%%%%%% Inputs to the model %%%%%%%%%%%%%%%%%%%%%%%%%%
s_frequency = 1;    % Spatial frequency in cycles per degree
orientation = 0;    % Orientation of grating in degrees
area = pi*(1)^2;           % Area of stimulus in visual sq. degrees
eccentricity = 0;   % Retinal eccentricity in degrees
t_frequency = linspace( 0, 20 )'; %Hz, must be a column vector

luminance = 3;                    % Mean luminance of background in cd/m^2
xy_background = [0.3127, 0.3290];   % xy chromaticity coordinates of stimulus background

luminance_modulation = 5;         % Peak luminance of modulation; Same as background luminance in case of isoluminant stimuli
xy_modulation = [0.3127, 0.3290];   % Direction ( in xy chromaticity coordinates ) of stimulus chromatic modulation

% example chromaticity coordinates: 
%   white/black: [0.3127, 0.3290]
%   red/green:  [0.5027, 0.2322]
%   yellow/violet: [0.3313, 0.3738]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lms_background = xyz2lms2006(Yxy2XYZ([luminance, xy_background]));
lms_modulation_peak = xyz2lms2006(Yxy2XYZ([luminance_modulation, xy_modulation]));
lms_modulation_delta = lms_modulation_peak-lms_background;
lms_delta_norm = lms_modulation_delta./norm(lms_modulation_delta);

% csf_pars = struct( 's_frequency',s_frequency, 't_frequency', t_frequency,... 
%     'orientation', orientation, 'area', area, 'eccentricity', eccentricity,...
%     'lms_bkg', lms_background,...
%     'lms_delta', lms_delta_norm);     

csf_pars = struct( 's_frequency',s_frequency, 't_frequency', t_frequency,... 
    'orientation', orientation, 'area', area, 'eccentricity', eccentricity,...
    'luminance', luminance);  

S = csf_model.sensitivity( csf_pars );        

plot( t_frequency, S );
set( gca, 'YScale', 'log' );
xlabel( 'Temporal frequency [Hz]' );
ylabel( 'Sensitivity' );

grid on

rmpath('../../matlab/');
rmpath('../../utils/');
