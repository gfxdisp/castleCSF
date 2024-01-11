classdef CSF_PostReceptContrast < CSF_base
    % A model from: 
    % Mantiuk, R. K., Kim, M., Ashraf, M., Xu, Q., Luo, M. R., Martinovic, J., & Wuerger, S. (2020, January). 
	% Practical Color Contrast Sensitivity Functions for Luminance Levels up to 10000 cd/m 2. 
	% In Color and Imaging Conference (Vol. 28, No. 1, pp. 1-6). Society for Imaging Science & Technology.
    
    properties( Constant )
        % which entries in the meachism matrix should be fixed to 1
        Mones = [ 1 0 0;
            1 0 0;
            0 0 1 ];
        chrom_ch_beta = 2; 
    end
    
    methods
        
        function obj = CSF_PostReceptContrast( fitted_par_vector_file )
            
            obj.par = obj.get_default_par();
            
            if exist( 'fitted_par_vector_file', 'var' )
                lv = load( fitted_par_vector_file );
                obj.par = obj.param2struct( obj.par, lv.fitted_par_vector );
            end
            
        end
        
        function name = short_name( obj )
            % A short name that could be used as a part of a file name
            name = 'post-recept-contrast';
        end
        
        function name = full_name( obj )
            name = 'Postreceptoral contrast';
        end
        
        function M_lms2acc = get_lms2acc( obj )
            % Get the colour mechanism matrix
            
            M_lms2acc = ones(3,3);
            % Set only the cells in the array that should be ~= 1
            M_lms2acc(~CSF_PostReceptContrast.Mones(:)) = obj.par.colmat;
            % Get the right sign
            M_lms2acc =  M_lms2acc .* [ 1 1 1; 1 -1 1; -1 -1 1];
        end
        
        function [P, C] = pdet( obj, csf_pars, LMS_mean, LMS_delta )
            % Predict the probility of detecting a Gabour patch of certain chromatic
            % direction and amplitide
            %
            % [P, C] = camliv_colour_difference( freq, LMS_mean, LMS_delta, area, params )
            %
            % freq - spatial frequency in cpd
            % LMS_mean - LMS of the background colour (CIE2006 CMF)
            % LMS_delta - colour direction vector in the LMS space (LMS_peak-LMS_mean)
            % area - area in deg^2
            %
            % The method returns:
            % P - The probability of detection
            % C - Normalized detection contrast (1 when P=0.5)            
            %            M_lms2acc = [ 1.0000    1.0000         0
            %                1.0000   -2.3112         0
            %                -1.0000   -1.0000   50.9875 ];
            
            M_lms2acc = obj.get_lms2acc();
            
%             ACC_mean = abs(reshape(LMS_mean, numel(LMS_mean)/3, 3) * M_lms2acc');
%             ACC_delta = abs(reshape(LMS_delta, numel(LMS_delta)/3, 3) * M_lms2acc');
            
            if (numel(size(LMS_mean)) > 2)  && (numel(LMS_mean)~=3)
                dim3_size = size(LMS_mean);
                dim1_size = dim3_size(1:end-1);
            else
                dim1_size = [max( size(LMS_mean, 1), size(LMS_delta, 1)), 1];
            end

            % Post-receptoral
            ACC_mean = abs(reshape(LMS_mean, numel(LMS_mean)/3, 3) * M_lms2acc');
            ACC_delta = abs(reshape(LMS_delta, numel(LMS_delta)/3, 3) * M_lms2acc');

            alpha = 0;
            C_A = reshape(abs(ACC_delta(:,1)./ACC_mean(:,1)),dim1_size);
            C_R = reshape(abs(ACC_delta(:,2)./...
                (alpha*ACC_mean(:,2) + (1-alpha)*ACC_mean(:,1)) ),dim1_size);
            C_Y = reshape(abs(ACC_delta(:,3)./...
                (alpha*ACC_mean(:,3) + (1-alpha)*ACC_mean(:,1))),dim1_size);


            csf_pars = obj.test_complete_params(csf_pars, { 'luminance', 'ge_sigma' }, true );
        
            lms_bkg = csf_pars.lms_bkg;
            lms_delta = csf_pars.lms_delta;

            sigma = csf_pars.ge_sigma;
            rho = csf_pars.s_frequency;
            lum = csf_pars.luminance;

            A = pi*(sigma).^2; % Stimulus area

            
            C_A_n = C_A.*obj.csf_acc(1, rho, A, lum, obj.par.ach);
            C_R_n = C_R.*obj.csf_acc(2, rho, A, lum, obj.par.rg);
            C_Y_n = C_Y.*obj.csf_acc(3, rho, A, lum, obj.par.yv);
            
            C = (C_A_n.^obj.beta + C_R_n.^obj.beta + C_Y_n.^obj.beta).^(1/obj.beta);
            
            P = 1 - exp( log(0.5)*C );
            
        end
        
        function S = sensitivity( obj, csf_pars )
            
            csf_pars = obj.test_complete_params(csf_pars, { 'lms_bkg', 'ge_sigma' }, true );

            k = obj.det_threshold( csf_pars );

            LMS_delta_thr = k .* csf_pars.lms_delta;
            last_dim = numel(size(LMS_delta_thr));  % the LMS channels are stored in the last dimension of the tensor
            S = 1./ (sqrt(sum((LMS_delta_thr ./ csf_pars.lms_bkg).^2, last_dim))/sqrt(3));

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
                
                % Optimized area and beta params from Gabor vs disc csf paper
                csf_pars.area = 2.42437; % Replace area parameter with fixed optimized area
                beta = 3.01142;
                
                S_gabor = sensitivity(obj, csf_pars);
                S1 = S_gabor.* (radius'.^(1/beta));
                S1 = max(S1);

                S = permute(S1, circshift(1:numel(size(S1)), -1)); 

        end

        % Returns the detection threshold as k * lms_delta
        % The function retuns k
        function k_thr = det_threshold( obj, csf_pars )
            
            csf_pars = obj.test_complete_params(csf_pars, { 'luminance', 'ge_sigma' }, true );
        
            lms_bkg = csf_pars.lms_bkg;
            lms_delta = csf_pars.lms_delta;
                        
            [C_A, C_R, C_Y] = csf_chrom_directions( obj, lms_bkg, lms_delta);
            
%             obj = obj.update_parameters();
            
            sigma = csf_pars.ge_sigma;
            rho = csf_pars.s_frequency;
            lum = csf_pars.luminance;

            A = pi*(sigma).^2; % Stimulus area

            C_A_n = C_A.*obj.csf_acc(1, rho, A, lum, obj.par.ach);
            C_R_n = C_R.*obj.csf_acc(2, rho, A, lum, obj.par.rg);
            C_Y_n = C_Y.*obj.csf_acc(3, rho, A, lum, obj.par.yv);
            
            C = (C_A_n.^obj.chrom_ch_beta + C_R_n.^obj.chrom_ch_beta + C_Y_n.^obj.chrom_ch_beta).^(1/obj.chrom_ch_beta);            

            k_thr = (C.^(-1));
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
            
            % Post-receptoral
            ACC_mean = abs(reshape(LMS_mean, numel(LMS_mean)/3, 3) * M_lms2acc');
            ACC_delta = abs(reshape(LMS_delta, numel(LMS_delta)/3, 3) * M_lms2acc');

            alpha = 0;
            C_A = reshape(abs(ACC_delta(:,1)./ACC_mean(:,1)),dim1_size);
            C_R = reshape(abs(ACC_delta(:,2)./...
                (alpha*ACC_mean(:,2) + (1-alpha)*ACC_mean(:,1)) ),dim1_size);
            C_Y = reshape(abs(ACC_delta(:,3)./...
                (alpha*ACC_mean(:,3) + (1-alpha)*ACC_mean(:,1))),dim1_size);
            

        end
        
        
        
        function S = csf_acc( obj, col, freq, area, lum, acc_pars )
            % Internal. Do not call from outside the object.
            % A nested CSF as a function of luminance

            S_max = obj.get_lum_dep( acc_pars.S_max, lum );
            f_max = obj.get_lum_dep( acc_pars.f_max, lum);
           gamma = obj.get_lum_dep( acc_pars.gamma, lum );
            bw = obj.get_lum_dep(acc_pars.bw, lum);

            % Truncated log-parabola for chromatic directions
            S_LP = 10.^( -abs(log10(freq) - log10(f_max)).^2./(2.^bw) );
            if col ~= 1
                ss = (freq<f_max);
            else
                ss = logical(zeros(size(freq)));
            end

            max_mat = repmat(max(S_LP), size(freq));
            if numel(max_mat) == 1
                S_LP(ss) = max_mat;
            else
                S_LP(ss) = max_mat(ss);
            end

            S_peak = S_max .* S_LP;


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
            if isfield( acc_pars, 'A_0' )
                A0 = acc_pars.A_0;
            else
                A0 = 270;
            end
% 
%             Ac = A0./(1+(freq/f0).^2);
% 
%             S = S_peak .* sqrt( Ac ./ (1+Ac./area)).*(freq);
            
            if 0
                Ac = A0./(1+(freq/f0).^2);
                S = S_peak .* sqrt( Ac ./ (1+Ac./area)).*(freq);
            else
                k = acc_pars.Ac_prime + area.*f0;
                S = S_peak .* sqrt( area.^gamma.*freq.^2 ./ (k + area.^gamma.*freq.^2) );
            end
        end

        
        
        
    end
    
    
    methods( Static )
        
        function p = get_default_par()
            
            p = CSF_base.get_dataset_par();

			p.ach.S_max = [ 521.522 4.28047 0.205589 14370 0.727125 ];
			p.ach.f_max = [ 1.63716 208.669 0.217275 ];
			p.ach.bw = 0.100024;
			p.ach.gamma = 0.939377;
			p.ach.Ac_prime = 58.0308;
			
			p.rg.S_max = [ 1656.93 21.116 0.516516 ];
			p.rg.f_max = 0.0437323;
			p.rg.bw = 2.25658;
			p.rg.gamma = 1.61491;
			p.rg.Ac_prime = 1.98645;
			
			p.yv.S_max = [ 61559.5 34.6744 0.456428 ];
			p.yv.f_max = 0.000904801;
			p.yv.bw = 3.02109;
			p.yv.gamma = 1.40366;
			p.yv.Ac_prime = 2.30455;
			
			p.colmat = [ 0.0345069 3.05716 2.25236 0.0625121 3.89577e-06 0.853869 ]; 
			
        end
        
    end
end

