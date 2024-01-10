classdef CSF_Pyramid_chrom < CSF_base

    % A model from: 
    % Watson and Ahumada, "The pyramid of visibility", 2016
    % https://www.researchgate.net/publication/305492593_The_pyramid_of_visibility
    
    methods
        
        function obj = CSF_Pyramid_chrom( colour )
            
            obj.par = obj.get_default_par( colour );                       
        end
        
        function name = short_name( obj )
            % A short name that could be used as a part of a file name
            name = 'csf-pyramid-chrom';
        end
        
        function name = full_name( obj )
            name = 'Chromatic Pyramid of Visibility CSF';
        end
                        
        function S = sensitivity( obj, csf_pars )

            csf_pars = obj.test_complete_params(csf_pars, { 'luminance', 'area' }, true );
            
            s_freq = csf_pars.s_frequency;
            t_freq = csf_pars.t_frequency;
            l_adapt = csf_pars.luminance;
            area = csf_pars.area;

            %[s_freq, t_freq, orientation, LMS_mean, LMS_delta, area, eccentricity] = obj.test_and_expand_pars(s_freq, t_freq, orientation, LMS_mean, LMS_delta, area, eccentricity);
            
            % Best-fit parameters from Fig. 4 of Pyramid of Visibility paper
%             c_w = -0.06;
%             c_f = -0.05;
%             c_i = 0.5;
%             c_0 = 2;
            
            % to be fixed
            c_w = -obj.par.c_w_neg;
            c_f = -obj.par.c_f_neg;
            c_i = obj.par.c_i;
            c_0 = obj.par.c_0;
            
            % default value from paper
            age = 30; 
            % pyramid operates on illuminance in Trolands
            d = pupil_d_unified( l_adapt, area, age );
            Illuminance = cdms2trolands(l_adapt, d);
            
            % using formula (4) from paper
            S = c_0 + c_w.*abs(t_freq) + c_f.*abs(s_freq) + c_i.*log10(Illuminance);

            % converting from log-sensitivity to linear
            S = 10.^S;
            
            % scaling by parameter
            P = obj.par.P;
            S = S.*P;
                        
        end

        function S = sensitivity_edge( obj, csf_pars )
            csf_pars.s_frequency = 0;
            S1 = obj.sensitivity(csf_pars);
            S = permute(S1, circshift(1:numel(size(S1)), -1)); 
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
%             p.P = 1.00272;
%             p.c_w_neg = 0.0202649;
%             p.c_f_neg = 0.0429971;
%             p.c_i = 0.258398;
%             p.c_0 = 1.01655;
		
			switch colour
				case 'rg'
					p.P = 0.948289;
	                p.c_w_neg = 0.0481838;
	                p.c_f_neg = 0.0529471;
	                p.c_i = 0.299123;
	                p.c_0 = 0.570596;
				case 'yv'
					p.P = 0.0156404;
	                p.c_w_neg = 0.039198;
	                p.c_f_neg = 0.0695055;
	                p.c_i = 0.301805;
	                p.c_0 = 0.06059;
				otherwise
					error('Invalid colour direction supplied');
			end
		

        end
        
    end
end

