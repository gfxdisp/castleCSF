classdef CSF_Pyramid_comb < CSF_base

    % A model from: 
    % Watson A.B. (2021). Stanford Talks: The Chromatic Pyramid of Visibility. 
	% https://talks.stanford.edu/watson-apple-chromatic-pyramid-of-visibililty/
    
    properties( Constant )
       % which entries in the meachism matrix should be fixed to 1
        Mones = [ 1 1 0;
            1 0 0;
            1 1 0 ];
        
        % LMS to DKL matrix parameters used in Wuerger, 2020 (JoV)
        colmat = [2.3112 0 0 50.9875];
        chrom_ch_beta = 2;
        
    end

    properties
        use_gpu = true;

        pyramid_ach = [];
        pyramid_rg = [];
        pyramid_yv = [];
    end

    methods
        
        function obj = CSF_Pyramid_comb( fitted_par_vector_file )
            
            obj.pyramid_ach = CSF_Pyramid();
            obj.pyramid_rg = CSF_Pyramid_chrom('rg');
            obj.pyramid_yv = CSF_Pyramid_chrom('yv');

            obj.par = obj.get_default_par();   

            if exist( 'fitted_par_vector_file', 'var' )
                lv = load( fitted_par_vector_file );
                obj.par = obj.param2struct( obj.par, lv.fitted_par_vector );
            end
        end
        
        function name = short_name( obj )
            % A short name that could be used as a part of a file name
            name = 'csf-pyramid-comb';
        end
        
        function name = full_name( obj )
            name = 'Combined Chromatic Pyramid of Visibility CSF';
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

            csf_pars = obj.test_complete_params(csf_pars, { 'lms_bkg', 'area' }, true );
            
            k = obj.det_threshold( csf_pars );

            LMS_delta_thr = k .* csf_pars.lms_delta;
            last_dim = numel(size(LMS_delta_thr));  % the LMS channels are stored in the last dimension of the tensor
            S = 1./ (sqrt(sum((LMS_delta_thr ./ csf_pars.lms_bkg).^2, last_dim))/sqrt(3));

        end

        function S = sensitivity_edge( obj, csf_pars )
            csf_pars.s_frequency = 0;
            S1 = obj.sensitivity(csf_pars);
            S = permute(S1, circshift(1:numel(size(S1)), -1)); 
        end

         % Returns the detection threshold as k * lms_delta
        % The function retuns k
        function k_thr = det_threshold( obj, csf_pars )
            
            csf_pars = obj.test_complete_params(csf_pars, { 'lms_bkg', 'ge_sigma' }, true );
        
            lms_bkg = csf_pars.lms_bkg;
            lms_delta = csf_pars.lms_delta;
                        
            [C_A, C_R, C_Y] = csf_chrom_directions( obj, lms_bkg, lms_delta);
            
            obj = obj.update_parameters();
            
            C_A_n = C_A.*obj.pyramid_ach.sensitivity(csf_pars);
            C_R_n = C_R.*obj.pyramid_rg.sensitivity(csf_pars);
            C_Y_n = C_Y.*obj.pyramid_yv.sensitivity(csf_pars);
            
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

            if 1
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
                C_A = reshape(abs(ACC_delta(:,1)./ACC_mean(:,1)),dim1_size);
                C_R = reshape(abs(ACC_delta(:,2)./...
                    (alpha*ACC_mean(:,2) + (1-alpha)*ACC_mean(:,1)) ),dim1_size);
                C_Y = reshape(abs(ACC_delta(:,3)./...
                    (alpha*ACC_mean(:,3) + (1-alpha)*ACC_mean(:,1))),dim1_size);
            end

        end

        function obj = update_parameters (obj)
            fn_ach = fieldnames(obj.pyramid_ach.par);
            fn_rg = fieldnames(obj.pyramid_rg.par);
            fn_yv = fieldnames(obj.pyramid_yv.par);
            
            for ff = 1:length(fn_ach)
               if isfield(obj.par, ([fn_ach{ff},'_ach'])) 
                  obj.pyramid_ach.par.(fn_ach{ff}) = obj.par.([fn_ach{ff},'_ach']);
               end
            end

            for ff = 1:length(fn_rg)
               if isfield(obj.par, ([fn_rg{ff},'_rg'])) 
                  obj.pyramid_rg.par.(fn_rg{ff}) = obj.par.([fn_rg{ff},'_rg']);
               end
            end
            
            for ff = 1:length(fn_yv)
               if isfield(obj.par, ([fn_yv{ff},'_yv'])) 
                  obj.pyramid_yv.par.(fn_yv{ff}) = obj.par.([fn_yv{ff},'_yv']);
               end
            end
            
        end

        function obj = print( obj, fh )
            % Print the model parameters in a format ready to be pasted into
            % get_default_par()
            
            % Printed formatted parameters for subclasses
            
            obj = obj.update_parameters();
            
            fn = fieldnames( obj.par );
            for ff=1:length(fn)
                if ismember( fn{ff}, { 'ds' } )
                    continue;
                end
                fprintf( fh, '\t\t\tp.%s = ', fn{ff} );
                obj.print_vector( fh, obj.par.(fn{ff}) );
                fprintf( fh, ';\n' );
            end

            fprintf(fh, 'Parameters for Ach subclass\n');
            obj.pyramid_ach.print(1);

            fprintf(fh, 'Parameters for RG subclass\n');
            obj.pyramid_rg.print(1);
            
            fprintf(fh, 'Parameters for YV subclass\n');
            obj.pyramid_yv.print(1);
            

        end

    end
    
    methods( Static )
        
        function p = get_default_par( colour )
            
            p = CSF_base.get_dataset_par();
            

            % Inital parameters (from the paper)
%             p.P = 1;
%             p.c_w_neg = 0.06;
%             p.c_f_neg = 0.05;
%             p.c_i = 0.5;
%             p.c_0 = 2;
            
            % Fitted to the data

            p.P_ach = 0.681333;
			p.c_w_neg_ach = 0.0396845;
			p.c_f_neg_ach = 0.0538477;
			p.c_i_ach = 0.403671;
			p.c_0_ach = 0.552272;
			p.P_rg = 0.895066;
			p.c_w_neg_rg = 0.0436539;
			p.c_f_neg_rg = 0.0544095;
			p.c_i_rg = 0.374631;
			p.c_0_rg = 0.402751;
			p.P_yv = 0.0198153;
			p.c_w_neg_yv = 0.0472675;
			p.c_f_neg_yv = 0.122009;
			p.c_i_yv = 0.32123;
			p.c_0_yv = 0.0620767;
		
        end
        
    end
end

