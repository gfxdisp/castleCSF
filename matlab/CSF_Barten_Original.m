classdef CSF_Barten_Original < CSF_base

    % Original implementation of the spatiotemporal Barten CSF derived from
    % his Ph.D. thesis (1999) including the extention to peripheral and temporal vision
    % Barten, P. G. (1999). Contrast sensitivity of the human eye and its effects on image quality. SPIE press.
    
    
    methods
        
        function obj = CSF_Barten_Original( )
            
            obj.par = obj.get_default_par();                       
        end
        
        function name = short_name( obj )
            % A short name that could be used as a part of a file name
            name = 'original-Barten-1999';
        end
        
        function name = full_name( obj )
            name = 'original Barten''s CSF (1999)';
        end
                        
        function S = sensitivity( obj, csf_pars )
            % u: Spatial frequency range in cpd
            % w: Temporal frequency range in HZ
            % e: Eccentricity in degrees
            % L: Average luminance of the observed object in cd/m^2
            % X0: Angular field size of object in degrees
            % k: Signal to noise ratio
            % eta0: Constant for quantom efficiency
            % sigma0: Constant for the eye MTF
            % eg: Eccentricity constant (can be different for various subjects)
            % u00: Spatial frequency above which the lateral inhibition ceases in the fovea
            % Phi00: Spectral density of the neural noise in the fovea
            % T: Eye integration time in sec
            % Xmax0: Maximum angular size of the integration area of the noise in the fovea
            % Nmax: Maximum number of cycles over which the eye can integrate the information
            % tau10: Time constant for the temporal filtering of the cones signal
            % tau20: Time constant for the temporal filtering of the lateral inhibition signal
            % n1: Asymptotic slope on double logarithmic scale for the cones signal
            % n2: Asymptotic slope on double logarithmic scale for the lateral inhibition signal

            csf_pars  = obj.test_complete_params(csf_pars, { 'luminance', 'area' }, true );
            u         = csf_pars.s_frequency;
            w         = csf_pars.t_frequency;
            e         = csf_pars.eccentricity;
            L         = csf_pars.luminance;
            X0        = sqrt(csf_pars.area);
            k         = get_field_def( obj.par, 'k', 3 );
            eta0      = get_field_def( obj.par, 'eta0', 0.03 );
            sigma0    = get_field_def( obj.par, 'sigma0', 0.50 );
            eg        = get_field_def( obj.par, 'eg', 3.3 );
            u00       = get_field_def( obj.par, 'u00', 7 );
            Phi00     = get_field_def( obj.par, 'phi00', 3e-8 );
            T         = get_field_def( obj.par, 'T', 0.1 );
            Xmax0     = get_field_def( obj.par, 'Xmax0', 12 );
            Nmax      = get_field_def( obj.par, 'Nmax', 15 );
            tau10     = get_field_def( obj.par, 'tau10', 0.032 );
            tau20     = get_field_def( obj.par, 'tau20', 0.018 );
            n1        = get_field_def( obj.par, 'n1', 7 );
            n2        = get_field_def( obj.par, 'n2', 4 );

            
            % Density of parasol Retinal ganglion cells in fovea
            Ng0      = 0.05*36000;
            
            % Density of parasol Retinal ganglion cells as a function of
            % eccentricity (estimated from primate retina)
            Ng       = Ng0.*(0.85./(1+(e/0.45).^2)+0.15./(1+(e./eg).^2));
            
            % Constants for calculation of sigma based on pupil diameter
            Cab       = 0.08;
                        
            % Maximum angular size of the integration area of the noise as a function of eccentricity
            Xmax      = Xmax0*(0.85./(1+(e/4).^2)+0.15./(1+(e/12).^2)).^-0.5;
            Ymax      = Xmax;
            
            % Quantum efficiency of the eye as a function of eccentricity
            eta       = eta0*(0.4./(1+(e/7).^2)+0.48./(1+(e/20).^2)+0.12);
            
            % Spectral density of the neural noise as a function of eccentricity
            Phi0      = Phi00*(Ng0./Ng);
            
            % Spatial frequency above which the lateral inhibition ceases as a function of eccentricity
            u0        = u00.*(Ng/Ng0).^0.5.*(0.85./(1+(e/4).^2)+0.13./(1+(e/20).^2)+0.02).^-0.5;
            
            % Photon conversion factor
            p         = 1.240*10^6;
            
            % Tangential field size
            Y0        = X0;

            % Diameter of the stimulus - assumed to be circular
            D         = sqrt(4/pi*X0.^2);
            
            % Pupil diameter
            d         = 5-3.*tanh(0.4*log10(L.*(X0.^2/(40^2))));
            
            % Retinal illuminance in Troland
            E         = (pi*(d.^2)/4).*L.*(1-(d/9.7).^2+(d/12.4).^4);
            
            % Standard deviation of the line-spread function caused by the discrete sturcture of the retina
            sigmaret  = 1./(sqrt(7.2.*sqrt(3)*Ng/3600));
            
            % Standard deviation of the line-spread function caused by the other parts of eye
            sigma00   = sqrt(sigma0^2-0.4^2);
            
            % Standard deviation of the line-spread function resulting from the convolution  of the different elements of convolution process
            sigma     = sqrt(sigmaret.^2+sigma00.^2+(Cab*d).^2);
            
            % Moduation transfer function of the optics of the eye 
            Mopt      = exp(-2*pi^2*sigma.^2.*u.^2.*(1/3600));
            
            % Temporal filter constants: time constants
            tau1      = tau10./(1+0.55.*log(1+(1+D/1.0).^0.6.*E/3.5));
            tau2      = tau20./(1+0.37.*log(1+(1+D/3.2).^0.6.*E/120));
            
            % Modulation transfer function of the temporal response of the eye
            H1        = 1./(1+(2*pi.*w.*tau1).^2).^(n1/2); % cone responses
            H2        = 1./(1+(2*pi.*w.*tau2).^2).^(n2/2); % lateral inhibation
            
            % Modulation transfer function of the lateral inhibition process
            F         = 1-sqrt(1-exp(-(u./u0).^2)); % high-pass filter characterizing lateral inhibation
            Mlat      = H1.*(1-H2.*F);
            
            % Spectral Density of photon noise
            PhiPhoton = 1./(eta.*p.*E);
            
            % Spatial equations
            X         = (X0.^-2+Xmax.^-2+(((0.5.*X0).^2+4*e.^2)./((0.5.*X0).^2)+e.^2).*(u.^2)./(Nmax^2)).^-0.5;
            Y         = (Y0.^-2+Ymax.^-2+(((0.5.*X0).^2)./((0.5.*X0).^2)+e.^2).*(u.^2)./(Nmax^2)).^-0.5;
            
            % Contrast Sensitivity Function
            S         = sqrt(2)*(Mopt./(2*k)).*sqrt((X.*Y.*T)./(PhiPhoton+Phi0./(Mlat.^2)));
        end

        function S = sensitivity_edge(obj, csf_pars)

			if 1
				% Maximum peak frequency model for discs from 
				% Ashraf, M., Mantiuk, R., & Chapiro, A. (2023). Modelling contrast sensitivity of discs. 
				%Electronic Imaging, 35, 1-8.
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
			else
				% Barten's reccomendation for handling disks
				% Fundumental frequency and third harmonic are inserted in F(u)
				first_harmonic       = 1./(sqrt(pi).*2.*csf_pars.ge_sigma);
				third_harmonic       = 3./(sqrt(pi).*2.*csf_pars.ge_sigma);
				small_disc           = csf_pars.ge_sigma < 2.5;
				large_disc           = csf_pars.ge_sigma >= 2.5;
				csf_pars.s_frequency = small_disc.*first_harmonic + large_disc.*third_harmonic;
										 
				S                    = sensitivity(obj, csf_pars);
																	 
											  
				S                    = permute(S, circshift(1:numel(size(S)), -1));
			end
        end
        
        
    end
    
    methods( Static )
        
        function p = get_default_par()
            
            p         = CSF_base.get_dataset_par();
            
            p.k = 6.85118;
			p.eta0 = 0.0354405;
			p.sigma0 = 0.52239;
			p.eg = 97.8251;
			p.u00 = 3.21885;
			p.Phi00 = 3e-08;
			p.T = 0.0661693;
			p.Xmax0 = 7.37859;
			p.Nmax = 15.7171;
			p.tau10 = 0.0441063;
			p.tau20 = 0.03434;
			p.n1 = 4.46403;
			p.n2 = 2.53013;

        end
        
    end
end

