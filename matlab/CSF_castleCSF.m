classdef CSF_castleCSF < CSF_base
    % Colour, Area, Spatial frequency, Temporal frequency, Luminance,
    % Eccentricity dependent Contrast Sensitivity Function (CSF)
    %
    % Refer to CSF_base.m for the documentation of the main interface
    % functions. 

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
        
        caslteCSF_ach = [];
        castleCSF_rg = [];
        castleCSF_yv = [];
    end
    
    methods
       
        function obj = CSF_castleCSF(  )
            
            obj.caslteCSF_ach = CSF_stelaCSF_lum_peak();
            obj.castleCSF_rg = CSF_chrom('rg');
            obj.castleCSF_yv = CSF_chrom('yv');
            obj.par = obj.get_default_par();
            
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
%                 colmat = [ 0.00123883 0.229778 0.932581 1.07013 6.41585e-07 0.0037047 ];
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
            
%             obj = obj.update_parameters();
            
            C_A_n = C_A.*obj.caslteCSF_ach.sensitivity(csf_pars);
            C_R_n = C_R.*obj.castleCSF_rg.sensitivity(csf_pars);
            C_Y_n = C_Y.*obj.castleCSF_yv.sensitivity(csf_pars);
            
            C = (C_A_n.^obj.chrom_ch_beta + C_R_n.^obj.chrom_ch_beta + C_Y_n.^obj.chrom_ch_beta).^(1/obj.chrom_ch_beta);            

            k_thr = (C.^(-1));
        end
        
        function S = sensitivity_edge(obj, csf_pars)
            
            t_freqs = csf_pars.t_frequency;
            tf_sel = t_freqs == 0;   % Check whether disk is static or temporally modulated
            
%             csf_pars_orig = csf_pars;
            
            if isempty(tf_sel)
                S = [];
                return;
            end
            
                csf_pars.s_frequency = logspace( log10(0.125), log10(16), 100 )';
%                 csf_pars.t_frequency = t_freqs(~tf_sel);
                csf_pars = obj.test_complete_params(csf_pars, { 'luminance', 'ge_sigma' } );
                radius = csf_pars.ge_sigma; % Store ge_sigma for multiple receptor circumference model
                radius = radius(:);
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
                S1 = S_gabor.* (radius'.^(1/beta));
                S1 = max(S1);

                S = permute(S1, circshift(1:numel(size(S1)), -1)); 

        end
        
        
        function [C_A, C_R, C_Y] = csf_chrom_directions( obj, LMS_mean, LMS_delta)
           
            M_lms2acc = obj.get_lms2acc();
            
%            lum = sum(LMS_mean,ndims(LMS_mean));
            
%             if (size(freq, 1) ~= size(LMS_mean, 1)) && (numel(freq)~=1) && (numel(LMS_mean)~=3)
            if (numel(size(LMS_mean)) > 2)  && (numel(LMS_mean)~=3)
                dim3_size = size(LMS_mean);
                dim1_size = dim3_size(1:end-1);
            else
                dim1_size = [max( size(LMS_mean, 1), size(LMS_delta, 1)), 1];
            end
            
            if 0
                % Cone contrast
                CC_LMS = LMS_delta ./ LMS_mean;            
                CC_ACC = reshape(CC_LMS, numel(CC_LMS)/3, 3) * M_lms2acc';

                C_A = reshape(abs(CC_ACC(:,1)), dim1_size);
                C_R = reshape(abs(CC_ACC(:,2)), dim1_size);
                C_Y = reshape(abs(CC_ACC(:,3)), dim1_size);
            
            else
                % Post-receptoral
                ACC_mean = abs(reshape(LMS_mean, numel(LMS_mean)/3, 3) * M_lms2acc');
                ACC_delta = abs(reshape(LMS_delta, numel(LMS_delta)/3, 3) * M_lms2acc');

                alpha = 0;
                C_A = reshape(abs(ACC_delta(:,1)./...
                    ACC_mean(:,1)),dim1_size);
                C_R = reshape(abs(ACC_delta(:,2)./...
                    (alpha*ACC_mean(:,2) + (1-alpha)*ACC_mean(:,1)) ),dim1_size);
                C_Y = reshape(abs(ACC_delta(:,3)./...
                    (alpha*ACC_mean(:,3) + (1-alpha)*ACC_mean(:,1))),dim1_size);
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
                    clf;
                    html_change_figure_print_size( gcf, 20, 10 );

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
                            %text( cm_dkl(2, cc), cm_dkl(3, cc), mech_label{cc}, 'Color', COLORs(cc,:) )
                            hold on
                        end
                        xlabel( dkl_label{dd(1)} );
                        ylabel( dkl_label{dd(2)} )
                        xlim( [-2, 2] )
                        ylim( [-2, 2] )
                        legend( hh, 'Location', 'SouthWest' );
                    end

                case 'sust_trans' % sust-trans-response
                    clf;
                    html_change_figure_print_size( gcf, 10, 10 );
                    omega = linspace( 0, 60 );
                    [R_sust, R_trans] = obj.caslteCSF_ach.get_sust_trans_resp(omega);
                    hh(1) = plot( omega, R_sust, '-k', 'DisplayName', 'Sustained (achromatic)');
                    hold on
                    hh(2) = plot( omega, R_trans, '--k', 'DisplayName', 'Transient (achromatic)');

                    R_sust = obj.castleCSF_rg.get_sust_trans_resp(omega);                    
                    hh(3) = plot( omega, R_sust, '-r', 'DisplayName', 'Sustained (red-green)');

                    R_sust = obj.castleCSF_yv.get_sust_trans_resp(omega);                    
                    hh(4) = plot( omega, R_sust, 'Color', [0.6 0 1], 'DisplayName', 'Transient (yellow-violet)');
                    hold off
                    xlabel( 'Temp. freq. [Hz]' );
                    ylabel( 'Response' );
                    legend( hh, 'Location', 'Best' );
                    grid on;

                case 'peak_s' 
                    obj.caslteCSF_ach.plot_mechanism(plt_id);

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

            obj.caslteCSF_ach.par = CSF_base.update_struct( obj.par.ach, obj.caslteCSF_ach.par );
            obj.castleCSF_rg.par = CSF_base.update_struct( obj.par.rg, obj.castleCSF_rg.par );
            obj.castleCSF_yv.par = CSF_base.update_struct( obj.par.yv, obj.castleCSF_yv.par );
            
        end
        
        function print( obj, fh )
            % Print the model parameters in a format ready to be pasted into
            % get_default_par()
                        
            M_lms2acc = obj.get_lms2acc();
            
            fprintf( fh, evalc( 'M_lms2acc' ) );
            
            obj.print@CSF_base(fh)

            % Printed formatted parameters for the component classes
            fprintf(fh, 'Parameters for Ach component:\n');
            obj.caslteCSF_ach.print(fh);
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

            p.rg.sigma_sust = 16.541;
            p.rg.beta_sust = 1.15549;
            p.rg.ch_sust.S_max = [ 730.147 38.2307 0.469349 ];
            p.rg.ch_sust.f_max = 0.0169475;
            p.rg.ch_sust.bw = 2.46441;
            p.rg.A_0 = 2879.13;
            p.rg.f_0 = 0.0704339;
            p.rg.ecc_drop = 0.0591431;
            p.rg.ecc_drop_nasal = 2.89648e-05;
            p.rg.ecc_drop_f = 2.04986e-69;
            p.rg.ecc_drop_f_nasal = 0.180118;
            p.yv.sigma_sust = 7.9187;
            p.yv.beta_sust = 0.999363;
            p.yv.ch_sust.S_max = [ 99.6752 62.5513 0.407922 ];
            p.yv.ch_sust.f_max = 0.00363324;
            p.yv.ch_sust.bw = 2.72495;
            p.yv.A_0 = 2.85765e+07;
            p.yv.f_0 = 0.000648528;
            p.yv.ecc_drop = 0.00357397;
            p.yv.ecc_drop_nasal = 5.85804e-141;
            p.yv.ecc_drop_f = 0.0080878;
            p.yv.ecc_drop_f_nasal = 0.0147658;
            p.ach.ach_sust.S_max = [ 56.0639 6.64673 0.144809 5.30245e-07 8.71574e+09 ];
            p.ach.ach_sust.f_max = [ 1.78263 66.992 0.269402 ];
            p.ach.ach_sust.bw = 0.000212711;
            p.ach.ach_sust.a = 0.282947;
            p.ach.ach_sust.A_0 = 157.103;
            p.ach.ach_sust.f_0 = 0.702338;
            p.ach.ach_trans.S_max = [ 0.193435 2793.56 ];
            p.ach.ach_trans.f_max = 0.000326927;
            p.ach.ach_trans.bw = 2.67727;
            p.ach.ach_trans.a = 0.000241177;
            p.ach.ach_trans.A_0 = 3.58559;
            p.ach.ach_trans.f_0 = 2.94741;
            p.ach.sigma_trans = 0.085424;
            p.ach.sigma_sust = 10.467;
            p.ach.omega_trans_sl = 2.36951;
            p.ach.omega_trans_c = 4.60295;
            p.ach.ecc_drop = 0.0259781;
            p.ach.ecc_drop_nasal = 0.0452708;
            p.ach.ecc_drop_f = 0.0217926;
            p.ach.ecc_drop_f_nasal = 0.0068348;
            p.colmat = [ 1.95127 0 0 86.9163 ];

        end
    end
end
