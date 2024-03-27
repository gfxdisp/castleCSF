classdef CSF_Wuerger2020 < CSF_base
    
    % A model from: 
    % Wuerger, S., Ashraf, M., Kim, M., Martinovic, J., Pérez-Ortiz, M., & Mantiuk, R. K. (2020). 
    % Spatio-chromatic contrast sensitivity under mesopic and photopic light levels. 
    % Journal of Vision, 20(4), 23-23.
    % https://doi.org/10.1167/jov.20.4.23

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
        
    end

    methods
        
        function obj = CSF_Wuerger2020( fitted_par_vector_file )
            
            obj.par = obj.get_default_par();
            
            if exist( 'fitted_par_vector_file', 'var' )
                lv = load( fitted_par_vector_file );
                obj.par = obj.param2struct( obj.par, lv.fitted_par_vector );
            end
        end
    
        function name = short_name( obj )
            % A short name that could be used as a part of a file name
            name = 'wuerger2020-csf';
        end

        function name = full_name( obj )
            name = 'Wuerger 2020 JOV CSF';
        end

        function S = sensitivity( obj, csf_pars, col_dir )
                      
            csf_pars = obj.test_complete_params(csf_pars, { 'luminance', 'ge_sigma' }, true );
        
                        
            sigma = csf_pars.ge_sigma;
            rho = csf_pars.s_frequency;
            lum = csf_pars.luminance;

            A = pi*(sigma).^2; % Stimulus area
            
            for cc = 1:3
                switch cc               
                    case 1
                        sel = col_dir == 1;
                        if ~isrow(sel)
                            S(sel, :) = obj.csf_acc(1, rho(sel, :), A(sel, :), lum(sel, :), obj.par.ach);
                        else
                            S(:,sel) = obj.csf_acc(1, rho, A, lum(:,sel), obj.par.ach);
                        end
                    case 2
                        sel = col_dir == 2;
                        if ~isrow(sel)
                            S(sel, :) = obj.csf_acc(2, rho(sel, :), A(sel, :), lum(sel, :), obj.par.rg);
                        else
                            S(:,sel) = obj.csf_acc(2, rho, A, lum(:,sel), obj.par.rg);
                        end
                        
                    case 3
                        sel = col_dir == 3;
                        if ~isrow(sel)
                            S(sel, :) = obj.csf_acc(3, rho(sel, :), A(sel, :), lum(sel, :), obj.par.yv);
                        else
                            S(:,sel) = obj.csf_acc(3, rho, A, lum(:,sel), obj.par.yv);
                        end
                        
                end
            end
            


        end
        

        function S = sensitivity_edge(obj, csf_pars, col_dir)
            
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
                
                % Optimized area and beta params from Gabor vs disc csf paper
                csf_pars.area = 2.42437; % Replace area parameter with fixed optimized area
                beta = 3.01142;
                
                S_gabor = sensitivity(obj, csf_pars, col_dir');
                S1 = S_gabor.* (radius'.^(1/beta));
                S1 = max(S1);

                S = permute(S1, circshift(1:numel(size(S1)), -1)); 

        end


        % Opponent-channel CSF model
        function S = csf_acc( obj, col, freq, area, lum, acc_pars )
            % Internal. Do not call from outside the object.
            % A nested CSF as a function of luminance


            S_peak = obj.get_log_parab(col, freq, lum, acc_pars);
            S_peak_lum = obj.get_log_parab(col, freq, 20, acc_pars);

%             S = S_peak;

            % The stimulus size model from the paper:
            %
            % Rovamo, J., Luntinen, O., & N�s�nen, R. (1993).
            % Modelling the dependence of contrast sensitivity on grating area and spatial frequency.
            % Vision Research, 33(18), 2773�2788.
            %
            % Equation on the page 2784, one after (25)

            if isfield( acc_pars, 'f_0' )
                f0 = acc_pars.f_0;
            else
                f0 = 0.65;
            end

            if isfield( acc_pars, 'Ac_prime' )
                Ac_prime = acc_pars.Ac_prime;
            else
                Ac_prime = 270;
            end

            if isfield( acc_pars, 'gamma' )
                gamma = acc_pars.gamma;
            else
                gamma = 1;
            end

            k = Ac_prime + area.*f0;
            S = S_peak .* sqrt( area.^gamma.*freq.^2 ./ (k + area.^gamma.*freq.^2) );

        end

        function S_peak = get_log_parab(obj, col, freq, lum, acc_pars)
            S_max = obj.get_lum_dep( acc_pars.S_max, lum, 'loggauss' );
            if col == 1
                f_max = obj.get_lum_dep( acc_pars.f_max, lum, 'gauss' );
            else
                f_max = obj.get_lum_dep( acc_pars.f_max, lum);
            end

            bw = obj.get_lum_dep(acc_pars.bw, lum);
            gamma = obj.get_lum_dep(acc_pars.gamma, lum);
            

            % Truncated log-parabola for chromatic directions
            S_peak = S_max ./ 10.^( (log10(freq) - log10(f_max)).^2./(0.5*2.^bw) );
            
            % low-pass for chromatic channels
            
            if col ~= 1
                ss = (freq<f_max);
            else
                ss = logical(zeros(size(freq)));
            end

            max_mat = repmat(max(S_peak.*(~ss)), size(freq));
            if numel(max_mat) == 1
                S_peak(ss) = max_mat;
            else
                S_peak(ss) = max_mat(ss);
            end

% 
%             ss = (freq<f_max);
%             max_mat = repmat(max(S_LP), size(freq));
%             if numel(max_mat) == 1
%                 S_LP(ss) = max_mat;
%             else
%                 S_LP(ss) = max_mat(ss);
%             end
% 
%             S_peak = S_max .* S_LP;


        end

    end

    methods( Static )
        
        function p = get_default_par()
            
            p = CSF_base.get_dataset_par();
            
			p.ach.S_max = [ 2.83912 1.28908 43.3254 ];
			p.ach.f_max = [ 2.53012 9.8288 44.2978 ];
			p.ach.bw = 1.57304;
			p.ach.Ac_prime = 115.858;
			p.ach.f_0 = 0.44229;
			p.ach.gamma = 0.879593;
			p.rg.S_max = [ 2.8927 3.17295 43.7516 ];
			p.rg.f_max = [ 0.0911927 3.04969e-24 ];
			p.rg.bw = 1.2271;
			p.rg.Ac_prime = 14.9649;
			p.rg.f_0 = 0.781478;
			p.rg.gamma = 1.45561;
			p.yv.S_max = [ 2.29543 2.99149 24.2307 ];
			p.yv.f_max = 0.0788081;
			p.yv.bw = 3.17419;
			p.yv.Ac_prime = 2.27589;
			p.yv.f_0 = 0.330337;
			p.yv.gamma = 1.38367;
    
        end
        
        function v = get_lum_dep( pars, L, name )
            % A family of functions modeling luminance dependency
            if nargin < 3
                name = '';
            end
            log_lum = log10(L);
            
            switch length(pars)
                case 1
                    % Constant
                    v = ones(size(L)) * pars(1);
                case 2
                    % Linear in log
                    v = 10.^(pars(1)*log_lum + pars(2));
                case 3
                    switch name
                        case 'gauss'
                            v = pars(1) * exp( -(log_lum-pars(2)).^2/pars(3) );
                        case 'loggauss'
                            v = 10.^(pars(1) * exp( -(log_lum-pars(2)).^2/pars(3) ));
                        otherwise
                            % A single hyperbolic function
                            v = pars(1)*(1+pars(2)./L).^(-pars(3));
                    end
                case 5
                    % Two hyperbolic functions
                    v = pars(1)*(1+pars(2)./L).^(-pars(3)) .* (1-(1+pars(4)./L).^(-pars(5)));
                otherwise
                    error( 'not implemented' );
            end
            
        end

    end


end