classdef CSF_castleCSF_chrom < CSF_base
    % The class predicts contrast sensitivity as the functioon of
    % Spatio-Temporal frequency, Eccentricity, Luminance and Area for
    % red-green and yellow-violet isoluminant stimuli

    properties( Constant )
        Y_min = 0.001;  % The minimum luminance
        Y_max = 10000;  % The maximum luminance
        rho_min = 2^-4  % The minimum spatial frequency
        rho_max = 64;   % The maximum spatial frequency
        ecc_max = 120;  % The maximum eccentricity
    end

    properties
        use_gpu = true;
        ps_beta = 1;
    end

    methods

        function obj = CSF_castleCSF_chrom( colour )

            obj.par = obj.get_default_par( colour );
        end

        function name = short_name( obj )
            % A short name that could be used as a part of a file name
            name = 'castle-csf-chrom';
        end

        function name = full_name( obj )
            name = 'castleCSF-chrom';
        end

        % Return contrast sensitivity for a given set of parameters. 
        % The sensitivity is assumed to be either the inverse of luminance
        % cotrast (L/\Delta L) or the inverse of cone contrast.
        %
        % This is a general interface that takes as input a
        % structure with the parameters values. This allows to add/remove
        % parameters without changing the interface, more flexible use of
        % parameters (e.g. either LMS or luminance of the background) and
        % should be less error prone. 
        % 
        % pars is the structure with the field names that closely match
        % the column names in the data files:
        %
        % luminance - luminance in cd/m^2. D65 background is assumed. 
        % lms_bkg - specify background colour and luminance 
        % s_frequency - spatial frequency in cpd
        % t_frequency - temporal frequency in Hz (default is 0)
        % orientation - orientation in deg (default is 0)
        % lms_delta - modulation direction in the LMS colour space (default
        %             D65 luminance modulation)
        % area - area of the stimulus in deg^2
        % ge_sigma - radius of the gaussian envelope in deg 
        % eccentricity - eccentricity in deg (default 0)       
        % vis_field - orientation in the visual field. See README.md
        %             (default 0)
        %
        % 'lms_bkg' and 'lms_delta' should be either [1 3] or [N 3] matrix. 
        % Any other field should be a column vector of size N or a scalar. 
        % For best performance, pass vectors with the the required number
        % of parameters. Do not call the sensitibity() function repetitively.
        % 
        %
        % You must specify 'luminance' or 'lms_bkg' but not both. 
        % You must specify 'area' or 'ge_sigma' but not both. 
        %
        % Example: 
        % csf_model = CSF_stelaCSF();
        % csf_pars = struct( 's_frequency', 4, 't_frequency', 1, 'orientation', 0, 'lms_bkg', [0.7443 0.3054 0.0157], 'area', 1, 'eccentricity', 0 );          
        % S = csf_model.sensitivity( csf_pars );        
        function S = sensitivity( obj, csf_pars )

            csf_pars = obj.test_complete_params(csf_pars, { 'luminance', 'ge_sigma' } );

            ecc = csf_pars.eccentricity;
            sigma = csf_pars.ge_sigma;
            rho = csf_pars.s_frequency;
            omega = csf_pars.t_frequency;
            lum = csf_pars.luminance;

            [R_sust] = get_sust_trans_resp(obj, omega);

            A = pi*(sigma).^2; % Stimulus area

            S_sust = obj.csf_chrom( rho, A, lum, ecc, obj.par.ch_sust );

            S = R_sust.*S_sust;           
            
            % The drop of sensitivity with the eccentricity (the window of
            % visibiliy model + extension)
            alpha = min(1, abs(csf_pars.vis_field-180)/90 );
            ecc_drop = alpha .* obj.par.ecc_drop + (1-alpha) .* obj.par.ecc_drop_nasal;
            ecc_drop_f = alpha .* obj.par.ecc_drop_f + (1-alpha) .* obj.par.ecc_drop_f_nasal;
            a = ecc_drop + rho.*ecc_drop_f;
            S = S .* 10.^(-a.*ecc);

        end
        
        function S = sensitivity_edge(obj, csf_pars)
            csf_pars.s_frequency = logspace( log10(0.125), log10(16), 100 )';
            csf_pars = obj.test_complete_params(csf_pars, { 'luminance', 'ge_sigma' } );
            radius = csf_pars.ge_sigma; % Store ge_sigma for multiple receptor circumference model
            radius = radius(:);
            csf_pars = rmfield(csf_pars, 'ge_sigma');
            csf_pars.area = 3.09781; % Replace area parameter with fixed optimized area
            
            beta =4;
            S_gabor = sensitivity(obj, csf_pars);
            S = S_gabor.* (radius'.^(1/beta));
            S = max(S);
            
            S = permute(S, circshift(1:numel(size(S)), -1)); 
        end
        
        % Get the sustained temporal response functions
        % omega - temporal frequency in Hz
        function [R_sust] = get_sust_trans_resp(obj, omega)
            sigma_sust = obj.par.sigma_sust;
            beta_sust = obj.par.beta_sust;
            R_sust = exp( -omega.^beta_sust / (sigma_sust) );
        end

        % Chromatic CSF model
        function S = csf_chrom( obj, freq, area, lum, ecc, ch_pars )
            % Internal. Do not call from outside the object.
            % A nested CSF as a function of luminance

            S_max = obj.get_lum_dep( ch_pars.S_max, lum );
            f_max = obj.get_lum_dep( ch_pars.f_max, lum );
            bw = ch_pars.bw;

            % Truncated log-parabola for chromatic directions
            S_LP = 10.^( -abs(log10(freq) - log10(f_max)).^2./(2^bw) );
            ss = (freq<f_max);
            max_S_LP = 1;
            S_LP(ss) = max_S_LP;

            S_peak = S_max .* S_LP;


            % The stimulus size model from the paper:
            %
            % Rovamo, J., Luntinen, O., & N�s�nen, R. (1993).
            % Modelling the dependence of contrast sensitivity on grating area and spatial frequency.
            % Vision Research, 33(18), 2773�2788.
            %
            % Equation on the page 2784, one after (25)

            if isfield( obj.par, 'f_0' )
                f0 = obj.par.f_0;
            else
                f0 = 0.65;
            end
            if isfield( obj.par, 'A_0' )
                A0 = obj.par.A_0;
            else
                A0 = 270;
            end

            Ac = A0./(1+(freq/f0).^2);

            S = S_peak .* sqrt( Ac ./ (1+Ac./area)).*(freq);
        end

        function pd = get_plot_description( obj )
            pd = struct();
            pp = 1;
            pd(pp).title = 'Sustained and transient response';
            pd(pp).id = 'sust_trans';
            pp = pp+1;
            pd(pp).title = 'Peak sensitivity';
            pd(pp).id = 'peak_s';
            pp = pp+1;
        end

        function plot_mechanism( obj, plt_id )
            switch( plt_id )
                case 'sust_trans' % sust-trans-response
                    clf;
                    html_change_figure_print_size( gcf, 10, 10 );
                    omega = linspace( 0, 100 );
                    [R_sust, R_trans] = obj.get_sust_trans_resp(omega);
                    hh(1) = plot( omega, R_sust, 'DisplayName', 'Sustained');
                    hold on
                    hh(2) = plot( omega, R_trans, 'DisplayName', 'Transient');
                    hold off
                    xlabel( 'Temp. freq. [Hz]' );
                    ylabel( 'Response' );
                    legend( hh, 'Location', 'Best' );
                    grid on;

                case { 'peak_s' }

                    f = logspace( -2, log10(5), 1024 );
                    L = logspace( -2, 4 );
                    [LL, ff] = meshgrid( L, f );
                    OMEGAs = [0 5 16];
                    COLORs = lines(length(OMEGAs));

                    for pp=1:length(OMEGAs)
                        csfpar.luminance = LL(:);
                        csfpar.s_frequency = ff(:);
                        csfpar.t_frequency = OMEGAs(pp);
                        csfpar.area = pi*(1.5).^2;
                        csfpar.eccentricity = 0;

                        S = obj.sensitivity( csfpar );

                        S = reshape( S, size(ff) );

                        S_max = max(S);

                        hh(pp) = plot( L, S_max, 'Color', COLORs(pp,:), 'DisplayName', sprintf('%g Hz', OMEGAs(pp)) );
                        hold on
                    end

                    L_dvr = logspace( -1, 0 );
                    hh(pp+1) = plot( L_dvr, sqrt(L_dvr)*50, '--k', 'DisplayName', 'DeVries-Rose law' );

                    set_axis_tick_label( 'x', 'luminance', L );
                    set_axis_tick_label( 'y', 'sensitivity', [1 1000]);

                    legend( hh, 'Location', 'best' );
                    ylabel( 'Peak sensitivity' );

                    ylim( [1 1000] );
                    grid on


                case { 'sust_peak_s', 'trans_peak_s' }
                    clf;
                    html_change_figure_print_size( gcf, 10, 10 );
                    L = logspace( -2, 4 );
                    if strcmp( plt_id, 'sust_peak_s' )
                        S_max = obj.par.ch_sust.S_max;
                    else
                        S_max = obj.par.ch_trans.S_max;
                    end
                    plot( L,  obj.get_lum_dep( S_max, L ) );
                    hold on
                    L_dvr = logspace( -1, 1 );
                    plot( L_dvr, sqrt(L_dvr)*100, '--k' );

                    set_axis_tick_label( 'x', 'luminance', L );
                    set_axis_tick_label( 'y', 'sensitivity', [1 100000] );
                    grid on;
                case 'sust_peak_f'
                    clf;
                    html_change_figure_print_size( gcf, 10, 10 );
                    L = logspace( -2, 4 );
                    f_max = obj.par.ch_sust.f_max;

                    plot( L,  obj.get_lum_dep( f_max, L ) );
                    set_axis_tick_label( 'x', 'luminance', L );
                    set_axis_tick_label( 'y', 'frequency', [0.01 60] );
                    grid on;
                otherwise
                    error( 'Wrong plt_id' );
            end
        end
        
        function print( obj, fh )
            % Print the model parameters in a format ready to be pasted into
            % get_default_par()
            
            fn = fieldnames( obj.par.ch_sust );
            for ff=1:length(fn)
                fprintf( fh, '\t\t\t\tp.ch_sust.%s = ', fn{ff} );
                obj.print_vector( fh, obj.par.ch_sust.(fn{ff}) );
                fprintf( fh, ';\n' );
            end
            fprintf( 1, '\n' )
            
            fn = fieldnames( obj.par );
            for ff=1:length(fn)
                if ismember( fn{ff}, { 'ch_sust', 'ds' } )
                    continue;
                end
                fprintf( fh, '\t\t\t\tp.%s = ', fn{ff} );
                obj.print_vector( fh, obj.par.(fn{ff}) );
                fprintf( fh, ';\n' );
            end
            
        end

    end

    methods( Static )

        function p = get_default_par( colour )

            p = CSF_base.get_dataset_par();
      
            switch colour
                case 'rg'
                    p.ch_sust.S_max = [ 681.434 38.0038 0.480386 ];
                    p.ch_sust.f_max = 0.0178364;
                    p.ch_sust.bw = 2.42104;
                    p.A_0 = 2816.44;
                    p.f_0 = 0.0711058;
                    p.sigma_sust = 16.4325;
                    p.beta_sust = 1.15591;
                    p.ecc_drop = 0.0591402;
                    p.ecc_drop_nasal = 2.89615e-05;
                    p.ecc_drop_f = 2.04986e-69;
                    p.ecc_drop_f_nasal = 0.18108;

                case 'yv'
                    p.ch_sust.S_max = [ 166.683 62.8974 0.41193 ];
                    p.ch_sust.f_max = 0.00425753;
                    p.ch_sust.bw = 2.68197;
                    p.A_0 = 2.82789e+07;
                    p.f_0 = 0.000635093;
                    p.sigma_sust = 7.15012;
                    p.beta_sust = 0.969123;
                    p.ecc_drop = 0.00356865;
                    p.ecc_drop_nasal = 5.85804e-141;
                    p.ecc_drop_f = 0.00806631;
                    p.ecc_drop_f_nasal = 0.0110662;

                otherwise
                    error('Invalid colour direction supplied');
            end
            

        end


    end

end