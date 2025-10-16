% Plot CSF functions

if ~exist( 'CSF_castleCSF', 'file' )
    addpath( fullfile( fileparts(mfilename('fullpath')), '..' ) );
end
addpath('../utils/');

csf_model = CSF_castleCSF();

% 3D plot - as the function of spatial and temporal frequency
figure(1);
clf

%%%%%%%%%%%%%%%%%%%% Inputs to the model %%%%%%%%%%%%%%%%%%%%%%%%%%

t_frequency = 0;           % Temporal frequency in Hz

orientation = 0;    % Orientation of grating in degrees
radius = 0.25;      % Radius of disc in visual degrees
eccentricity = 0;   % Retinal eccentricity in degrees

luminance = logspace(-2, 3, 100);                    % Mean luminance of background in cd/m^2
lms_background_norm = xyz2lms2006(whitepoint('d65'));
lms_background = luminance'.*lms_background_norm;

col_modulation_matrix = eye(3)/((csf_model.get_lms2acc)');
col_dir = 1;        % 1: achormatic, 2: red/green, 3: yellow/violet
lms_delta_norm = col_modulation_matrix(col_dir, :);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% castleCSF fully support broadcasting:
% s_frequency has the dimensions [1 N]
% t_frequency has the dimensions [M 1]
% The resulting S will have the dimension [M N] and provide sensitivity
% for all combinations of spatial and temporal frequencies
% Note that the LMS colour coordinates must be provides as the last
% dimension (3rd dimension) and hence reshape(lms_background, [1 1 3]).

csf_pars = struct( 't_frequency', t_frequency',... 
    'orientation', orientation, 'ge_sigma', radius, 'eccentricity', eccentricity,...
    'lms_bkg', reshape(lms_background, 1, [], 3),...
    'lms_delta', reshape(lms_delta_norm, [1 1 3]) );

S = csf_model.sensitivity_edge( csf_pars );    

plot(luminance, S);

set( gca, 'XScale', 'log' );
set( gca, 'YScale', 'log' );

grid on

xlabel( 'Luminance [cd/m^2]')
ylabel( 'Sensitivity')
title( 'castleCSF');

