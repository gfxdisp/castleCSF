classdef CSF_castleCSF < CSF_base
    % Colour, Area, Spatial frequency, Temporal frequency, Luminance,
    % Eccentricity dependent Contrast Sensitivity Function (CSF)
    %
    % Model from: 
    % Ashraf, M., Mantiuk, R. K., Chapiro, A., & Wuerger, S. (2024). 
    % castleCSF—A contrast sensitivity function of color, area, spatiotemporal 
    % frequency, luminance and eccentricity. Journal of Vision, 24(4), 5-5.
    %
    % Refer to CSF_base.m for the documentation of the main interface
    % functions. 
    %
    % Example: 
    % csf_model = CSF_castleCSF();
    % csf_pars = struct( 's_frequency', 4, 't_frequency', 1, 'orientation', 0,... 
    %   'lms_bkg', [0.7443 0.3054 0.0157], 'lms_delta', [0.9182, 0.3953, 0.0260],... 
    %   'area', 1, 'eccentricity', 0 );          
    % S = csf_model.sensitivity( csf_pars );   
    %
    % Subfunctions: CSF_stelaCSF_lum_peak, CSF_castleCSF_chrom

    properties( Constant )
       % which entries in the meachism matrix should be fixed to 1
        Mones = [ 1 1 0;
            1 0 0;
            1 1 0 ];
        
        % LMS to DKL matrix parameters used in Wuerger, 2020 (JoV)
        colmat = [2.3112 0 0 50.9875];
        chrom_ch_beta = 2;
        
        Y_min = 0.001;  % The minimum luminance
        Y_max = 10000;  % The maximum luminance
        rho_min = 2^-4  % The minimum spatial frequency
        rho_max = 64;   % The maximum spatial frequency
        ecc_max = 120;  % The maximum eccentricity        
        
    end
    
    properties  
        use_gpu = true;
        ps_beta = 1;
        
        castleCSF_ach = [];
        castleCSF_rg = [];
        castleCSF_yv = [];
    end
    
    methods
       
        function obj = CSF_castleCSF(  )
            
            obj.castleCSF_ach = CSF_stelaCSF_lum_peak();
            obj.castleCSF_rg = CSF_castleCSF_chrom('rg');
            obj.castleCSF_yv = CSF_castleCSF_chrom('yv');
            obj.par = obj.get_default_par();
            obj = obj.update_parameters(); % Update parameters in subclasses
            
        end
        
        function name = short_name( obj )
            % A short name that could be used as a part of a file name
            name = 'castle-csf';
        end

        function name = full_name( obj )
            name = 'castleCSF';
        end
        
         function M_lms2acc = get_lms2acc( obj )
            % Get the colour mechanism matrix
            
            M_lms2acc = ones(3,3);
            % Set only the cells in the array that should be ~= 1
            if isfield(obj.par, 'colmat')
                M_lms2acc(~obj.Mones(:)) = obj.par.colmat;
            else
                M_lms2acc(~obj.Mones(:)) = obj.colmat;
            end 
            % Get the right sign
            M_lms2acc =  M_lms2acc .* [ 1 1 1; 1 -1 1; -1 -1 1];
         end
        
        function S = sensitivity( obj, csf_pars )
            
            csf_pars = obj.test_complete_params(csf_pars, { 'lms_bkg', 'ge_sigma' }, true );

            k = obj.det_threshold( csf_pars );

            LMS_delta_thr = k .* csf_pars.lms_delta;
            last_dim = numel(size(LMS_delta_thr));  % the LMS channels are stored in the last dimension of the tensor
            S = 1./ (sqrt(sum((LMS_delta_thr ./ csf_pars.lms_bkg).^2, last_dim))/sqrt(3));

        end

        % Returns the detection threshold as k * lms_delta
        % The function retuns k
        function k_thr = det_threshold( obj, csf_pars )
            
            csf_pars = obj.test_complete_params(csf_pars, { 'lms_bkg', 'ge_sigma' }, true );
        
            lms_bkg = csf_pars.lms_bkg;
            lms_delta = csf_pars.lms_delta;
                        
            [C_A, C_R, C_Y] = csf_chrom_directions( obj, lms_bkg, lms_delta);
                        
            C_A_n = C_A.*obj.castleCSF_ach.sensitivity(csf_pars);
            C_R_n = C_R.*obj.castleCSF_rg.sensitivity(csf_pars);
            C_Y_n = C_Y.*obj.castleCSF_yv.sensitivity(csf_pars);
            
            C = (C_A_n.^obj.chrom_ch_beta + C_R_n.^obj.chrom_ch_beta + C_Y_n.^obj.chrom_ch_beta).^(1/obj.chrom_ch_beta);            

            k_thr = (C.^(-1));
        end
        
        function S = sensitivity_edge(obj, csf_pars)
            % Contrast sensitivity for discs, based on: 
            % "Modeling contrast sensitivity of discs", Maliha Ashraf,
            % Rafał K. Mantiuk and Alexandre Chapiro., http://dx.doi.org/10.2352/EI.2023.35.10.HVEI-246            
            
            t_freqs = csf_pars.t_frequency;
            tf_sel = t_freqs == 0;   % Check whether disk is static or temporally modulated
                        
            if isempty(tf_sel)
                S = [];
                return;
            end
            
            if isfield( csf_pars, 's_frequency' ) || isfield( csf_pars, 'area' )
                error( "Frequency or area cannot be specified when computing sensitivity for discs" );
            end
            csf_pars = obj.test_complete_params(csf_pars, { 'luminance', 'ge_sigma' } );

            % add extra dimension for spatial frequency
            fn = fieldnames(csf_pars);
            for kk=1:length(fn)
                if ~isscalar( csf_pars.(fn{kk}))
                    csf_pars.(fn{kk}) = reshape( csf_pars.(fn{kk}), [1 size(csf_pars.(fn{kk}))] );
                end
            end
            csf_pars.s_frequency = logspace( log10(0.125), log10(16), 64 )';

            radius = csf_pars.ge_sigma; % Store ge_sigma for multiple receptor circumference model
            csf_pars = rmfield(csf_pars, 'ge_sigma');
            
            if isfield(obj.par, 'disc_area')
                csf_pars.area = obj.par.disc_area;
            else
                % Optimized area and beta params from Gabor vs disc csf paper
                csf_pars.area = 2.42437; % Replace area parameter with fixed optimized area
            end
            if isfield(obj.par, 'disc_beta')
                beta = obj.par.disc_beta;
            else
                % Optimized area and beta params from Gabor vs disc csf paper
                beta = 3.01142;
            end
            
            S_gabor = sensitivity(obj, csf_pars);
            S1 = S_gabor.* (radius.^(1/beta));
            S1 = max(S1, [], 1);

            % Remove the 1st dimension, which we used for spatial
            % frequencies
            sz = size(S1);
            S = reshape( S1, sz(2:end) );
        end
        
        
        function [C_A, C_R, C_Y] = csf_chrom_directions( obj, LMS_mean, LMS_delta)
           
            M_lms2acc = obj.get_lms2acc();
            
            % if (numel(size(LMS_mean)) > 2)  && (numel(LMS_mean)~=3)
            %     dim3_size = size(LMS_mean);
            %     dim1_size = dim3_size(1:end-1);
            % else
            %     dim1_size = [max( size(LMS_mean, 1), size(LMS_delta, 1)), 1];
            % end

            target_size = max( [size(LMS_mean); size(LMS_delta)] );
            assert( target_size(end)==3 ); % The last dimension must be LMS
            target_size = target_size(1:(end-1)); % remove the last dim
 
            % Post-receptoral
            ACC_mean = reshape( abs(reshape(LMS_mean, numel(LMS_mean)/3, 3) * M_lms2acc'), size(LMS_mean) );
            ACC_delta = reshape( abs(reshape(LMS_delta, numel(LMS_delta)/3, 3) * M_lms2acc'), size(LMS_delta) );

            otherdims_delta = repmat({':'},1,ndims(ACC_delta)-1);
            otherdims_mean = repmat({':'},1,ndims(ACC_delta)-1);

            alpha = 0;
            C_A = obj.reshape_fix(abs(ACC_delta(otherdims_delta{:},1)./...
                ACC_mean(otherdims_mean{:},1)), target_size);
            C_R = obj.reshape_fix(abs(ACC_delta(otherdims_delta{:},2)./...
                (alpha*ACC_mean(otherdims_mean{:},2) + (1-alpha)*ACC_mean(otherdims_mean{:},1)) ), target_size);
            C_Y = obj.reshape_fix(abs(ACC_delta(otherdims_delta{:},3)./...
                (alpha*ACC_mean(otherdims_mean{:},3) + (1-alpha)*ACC_mean(otherdims_mean{:},1))), target_size);
        end
        
        % Handle reshape for cases when size = [1]
        function Y = reshape_fix( obj, X, sz )
            if isscalar(sz)
                Y = X;
            else
                Y = reshape( X, sz );
            end
        end
        
        function pd = get_plot_description( obj )
            pd = struct();
            pp = 1;

            pd(pp).title = 'Color mechanisms';
            pd(pp).id = 'col_mech';
            pp = pp+1;
            

            pd(pp).title = 'Sustained and transient response';
            pd(pp).id = 'sust_trans';
            pp = pp+1;

            pd(pp).title = 'Peak sensitivity';
            pd(pp).id = 'peak_s';
            pp = pp+1;

        end
        
        function plot_mechanism( obj, plt_id )
            switch( plt_id )
                case 'col_mech' 
                    figure,
                    M = obj.get_lms2acc();                    
                    cm_lms = inv(M)*eye(3);
                    cm_dkl = lms2dkl_d65(cm_lms);

                    mech_label = { 'achromatic', 'red-green', 'violet-yellow' };
                    dkl_label = { 'L+M', 'L-M', 'L-(M+S)' };
                    for pp=1:2
                        subplot( 1, 2, pp );
                        if pp==1
                            dd = [2 3];
                        else
                            dd = [2 1];
                        end
                        COLORs = lines(3);
                        hh = [];
                        for cc=1:3
                            hh(cc) = quiver( 0, 0, cm_dkl(dd(1), cc), cm_dkl(dd(2), cc), 'Color', COLORs(cc,:), 'DisplayName', mech_label{cc} );
                            hold on
                        end
                        xlabel( dkl_label{dd(1)} );
                        ylabel( dkl_label{dd(2)} )
                        xlim( [-2, 2] )
                        ylim( [-2, 2] )
                        legend( hh, 'Location', 'SouthWest' );
                    end

                case 'sust_trans' % sust-trans-response
                    figure,
                    omega = linspace( 0, 60 );
                    lums = [0.1 30 1000];
																								
                    hold on,
                    
                    for ll = 1:length(lums)
                        [R_sust, R_trans] = obj.castleCSF_ach.get_sust_trans_resp(omega, lums(ll));
                        if ll == 1
                           hh(1) = plot( omega, R_sust, '-k', 'DisplayName', 'Sustained (achromatic)');
                        end
                        hh(ll+1) = plot( omega, R_trans, 'LineStyle', '--',...
                            'DisplayName',... 
                            sprintf('Transient (achromatic) (%g cd/m^2)', lums(ll)));                      
                    end
                    
                    R_sust = obj.castleCSF_rg.get_sust_trans_resp(omega);                    
                    hh(5) = plot( omega, R_sust, '-r', 'DisplayName', 'Sustained (red-green)');

                    R_sust = obj.castleCSF_yv.get_sust_trans_resp(omega);                    
                    hh(6) = plot( omega, R_sust, 'Color', [0.6 0 1], 'DisplayName', 'Sustained (yellow-violet)');
                    hold off
                    xlabel( 'Temp. freq. [Hz]' );
                    ylabel( 'Response' );
                    legend( hh, 'Location', 'Best' );
                    grid on;

                case 'peak_s' 
                    figure,
                    obj.castleCSF_ach.plot_mechanism(plt_id);

                otherwise
                    error( 'Wrong plt_id' );
                    
            end
        end
        
        function obj = set_pars( obj, pars_vector )
            obj = obj.set_pars@CSF_base(pars_vector);
            obj = obj.update_parameters();
        end

        % Copy parameters to from this object to individual CSF components
        function obj = update_parameters (obj)

            obj.castleCSF_ach.par = CSF_base.update_struct( obj.par.ach, obj.castleCSF_ach.par );
            obj.castleCSF_rg.par = CSF_base.update_struct( obj.par.rg, obj.castleCSF_rg.par );
            obj.castleCSF_yv.par = CSF_base.update_struct( obj.par.yv, obj.castleCSF_yv.par );
            
        end
        
        function print( obj )
            % Print the model parameters in a format ready to be pasted into
            % get_default_par()
                        
            M_lms2acc = obj.get_lms2acc();
            fh = 1;
            
            fprintf( fh, evalc( 'M_lms2acc' ) );
            
            obj.print@CSF_base(fh)

            % Printed formatted parameters for the component classes
            fprintf(fh, 'Parameters for Ach component:\n');
            obj.castleCSF_ach.print(fh);
            fprintf(fh, '\n');

            fprintf(fh, 'Parameters for RG component:\n');
            obj.castleCSF_rg.print(fh);
            fprintf(fh, '\n');
            
            fprintf(fh, 'Parameters for YV component:\n');
            obj.castleCSF_yv.print(fh);
            fprintf(fh, '\n');
        end
    end
    
    methods ( Static )
       
        function p = get_default_par()
            % The list of trainable parameters and their default values. 

            p = CSF_base.get_dataset_par();

            p.rg.sigma_sust = 16.4325;
			p.rg.beta_sust = 1.15591;
			p.rg.ch_sust.S_max = [ 681.434 38.0038 0.480386 ];
			p.rg.ch_sust.f_max = 0.0178364;
			p.rg.ch_sust.bw = 2.42104;
			p.rg.A_0 = 2816.44;
			p.rg.f_0 = 0.0711058;
			p.rg.ecc_drop = 0.0591402;
			p.rg.ecc_drop_nasal = 2.89615e-05;
			p.rg.ecc_drop_f = 2.04986e-69;
			p.rg.ecc_drop_f_nasal = 0.18108;
			p.yv.sigma_sust = 7.15012;
			p.yv.beta_sust = 0.969123;
			p.yv.ch_sust.S_max = [ 166.683 62.8974 0.41193 ];
			p.yv.ch_sust.f_max = 0.00425753;
			p.yv.ch_sust.bw = 2.68197;
			p.yv.A_0 = 2.82789e+07;
			p.yv.f_0 = 0.000635093;
			p.yv.ecc_drop = 0.00356865;
			p.yv.ecc_drop_nasal = 5.85804e-141;
			p.yv.ecc_drop_f = 0.00806631;
			p.yv.ecc_drop_f_nasal = 0.0110662;
			p.ach.ach_sust.S_max = [ 56.4947 7.54726 0.144532 5.58341e-07 9.66862e+09 ];
			p.ach.ach_sust.f_max = [ 1.78119 91.5718 0.256682 ];
			p.ach.ach_sust.bw = 0.000213047;
			p.ach.ach_sust.a = 0.100207;
			p.ach.ach_sust.A_0 = 157.103;
			p.ach.ach_sust.f_0 = 0.702338;
			p.ach.ach_trans.S_max = [ 0.193434 2748.09 ];
			p.ach.ach_trans.f_max = 0.000316696;
			p.ach.ach_trans.bw = 2.6761;
			p.ach.ach_trans.a = 0.000241177;
			p.ach.ach_trans.A_0 = 3.81611;
			p.ach.ach_trans.f_0 = 3.01389;
			p.ach.sigma_trans = 0.0844836;
			p.ach.sigma_sust = 10.5795;
			p.ach.omega_trans_sl = 2.41482;
			p.ach.omega_trans_c = 4.7036;
			p.ach.ecc_drop = 0.0239853;
			p.ach.ecc_drop_nasal = 0.0400662;
			p.ach.ecc_drop_f = 0.0189038;
			p.ach.ecc_drop_f_nasal = 0.00813619;

        end
    end
end
