classdef CSF_A_Daly < CSF_base
    
    % A model from: 
    % Daly, S.J. “Visible Differences Predictor: An Algorithm for the Assessment of Image Fidelity.” 
    % In Digital Images and Human Vision, edited by Andrew B. Watson, 1666:179–206. MIT Press, 1993. 
    % https://doi.org/10.1117/12.135952.
    %
    % The model contains a few corrections (typos in the paper). 
    
    %     properties( Constant )
    %     end
    
    methods
        
        function obj = CSF_A_Daly( )
            
            obj.par = obj.get_default_par();
        end
        
        function name = short_name( obj )
            % A short name that could be used as a part of a file name
            name = 'daly-csf';
        end

        function name = full_name( obj )
            name = 'VDP CSF';
        end
                        
        function S = sensitivity( obj, csf_pars )

            csf_pars = obj.test_complete_params(csf_pars, { 'luminance', 'area' }, true );

            freq = csf_pars.s_frequency;
            d = 0.5; % Viewing distance in meters
            l_adapt = csf_pars.luminance;
            theta = csf_pars.orientation; % Orientation in degs
            area = csf_pars.area;
                        
            P=obj.par.P; % Peak CSF sensitivity, from Daly's paper
            epsilon=obj.par.epsilon;
            a_l_m = get_field_def( obj.par, 'a_l_m', 0.7 );
            A_l = 0.801*(1+a_l_m./l_adapt).^(-0.2);
            b_l_m = get_field_def( obj.par, 'b_l_m', 100 );
            B_l = 0.3*(1+b_l_m./l_adapt).^0.15;
            S = ( (3.23*(freq.^2.*area).^(-0.3)).^5 + 1).^(-1/5) .* ...
                A_l.*epsilon.*freq.*exp(-(B_l.*epsilon.*freq)).* ...
                sqrt( 1+0.06*exp(B_l.*epsilon.*freq) );
            
            r_a = 0.856 * d^0.14;
            ecc = csf_pars.eccentricity;   % Eccentrity in visual degrees
            k = get_field_def( obj.par, 'k', 0.24 );
            r_e = 1 ./ (1+k*ecc);
            ob = obj.par.ob;
            r_theta = (1-ob)/2 * cos(4*theta) + (1+ob)/2;
            
            ro_prime = freq ./ (r_a .* r_e .* r_theta);
            
            S_prime = ( (3.23*(ro_prime.^2.*area).^(-0.3)).^5 + 1).^(-1/5) .* ...
                A_l.*epsilon.*ro_prime.*exp(-(B_l.*epsilon.*ro_prime)).* ...
                sqrt( 1+0.06*exp(B_l.*epsilon.*ro_prime) );
            
            S = P*min( S, S_prime );
            S( isnan(S) | isinf(S) ) = 0;
                                    
        end   

        function S = sensitivity_edge(obj, csf_pars)

            csf_pars.s_frequency = logspace( log10(0.125), log10(16), 100 )';
            csf_pars             = obj.test_complete_params(csf_pars, { 'luminance', 'ge_sigma' } );
            radius               = csf_pars.ge_sigma; % Store ge_sigma for multiple receptor circumference model
            radius               = radius(:);
            csf_pars             = rmfield(csf_pars, 'ge_sigma');
            csf_pars.area        = 3.09781; % Replace area parameter with fixed optimized area
            beta                 = 4;
            S_gabor              = sensitivity(obj, csf_pars);
            S                    = S_gabor.* (radius'.^(1/beta));
            S                    = max(S);
            S                    = permute(S, circshift(1:numel(size(S)), -1)); 

        end
        
    end
    
    methods( Static )
        
        function p = get_default_par()
            
            p = CSF_base.get_dataset_par();
            p.P = 133.9;
			p.ob = 1.44757;
			p.k = 0.206707;
			p.epsilon = 1.39533;
			p.a_l_m = 1.967;
			p.b_l_m = 13.5589;
        end
        
    end
end

