classdef CSF_Barten_HF < CSF_base

    % spatiotemporal Barten CSF modified for high temporal frequencies (HF)
    % Spatiotemporal contrast sensitivity functions: predictions for the
    % critical flicker frequency, Human Vision Electronic Imaging 2024

    methods

        function obj = CSF_Barten_HF( )
            obj.par = obj.get_default_par();
        end

        function name = short_name( obj )
            % A short name that could be used as a part of a file name
            name = 'Barten-HF';
        end

        function name = full_name( obj )
            name = 'Barten CSF HF';
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
            % tau10: Time constant for the temporal filtering of the cones signal
            % tau20: Time constant for the temporal filtering of the lateral inhibition signal
            % n1: Asymptotic slope on double logarithmic scale for the cones signal (Temporal)
            % n2: Asymptotic slope on double logarithmic scale for the lateral inhibition signal (Temporal)

            csf_pars  = obj.test_complete_params(csf_pars, { 'luminance', 'ge_sigma' }, true );

            % Gabor or disk (based on input dimensions)
            dimension = size(csf_pars.s_frequency);
            if dimension(1,2) == 1
                aperture = "gabor";
            elseif dimension(1,1) == 1
                aperture = "disc";
            end

            u         = csf_pars.s_frequency;
            w         = csf_pars.t_frequency;
            e         = csf_pars.eccentricity;
            L         = csf_pars.luminance;
            X0        = sqrt(pi*(csf_pars.ge_sigma).^2);
            k         = get_field_def( obj.par, 'k', 3 );
            eta0      = get_field_def( obj.par, 'eta0', 0.03 );
            sigma0    = get_field_def( obj.par, 'sigma0', 0.50 );
            u0        = get_field_def( obj.par, 'u0', 7 );
            Phi00     = get_field_def( obj.par, 'Phi00', 3e-8 );
            T         = get_field_def( obj.par, 'T', 0.1 );
            Xmax0     = get_field_def( obj.par, 'Xmax0', 12 );
            Nmax      = get_field_def( obj.par, 'Nmax', 15 );
            tau10     = get_field_def( obj.par, 'tau10', 0.032 );
            tau20     = get_field_def( obj.par, 'tau20', 0.018);
            n1        = get_field_def( obj.par, 'n1', 7 );
            n2        = get_field_def( obj.par, 'n2', 4 );
            k_tem     = get_field_def( obj.par, 'k_tem', 0.24 );
            k_nas     = get_field_def( obj.par, 'k_nas', 0.24 );


            % Density of ganglion cells based on Watson's formula
            Ngm        = obj.ganglion_density(csf_pars, "midget");

            % Constants for calculation of sigma based on pupil diameter
            Cab       = 0.08;

            % Maximum angular size of the integration area of the noise as a function of eccentricity
            Xmax      = Xmax0*(0.85./(1+(e/4).^2)+0.15./(1+(e/12).^2)).^-0.5;
            Ymax      = Xmax;

            % Quantum efficiency of the eye as a function of eccentricity
            eta       = eta0*(0.4./(1+(e/7).^2)+0.48./(1+(e/20).^2)+0.12);

            % Spectral density of the neural noise as a function of eccentricity
            alpha     = min(1, abs(csf_pars.vis_field-180)/90 );
            M_slope   = alpha .* k_tem + (1-alpha) .* k_nas;
            Phi0      = Phi00*(1+M_slope.*e);

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
            sigmaret  = 1./(sqrt(7.2.*sqrt(3)*Ngm/3600));

            % Standard deviation of the line-spread function caused by the other parts of eye
            sigma00   = sqrt(sigma0^2-0.0312^2);

            % Standard deviation of the line-spread function resulting from the convolution  of the different elements of convolution process
            sigma     = sqrt(sigmaret.^2+sigma00.^2+(Cab*d).^2);

            % Moduation transfer function of the optics of the eye
            if aperture == "gabor"
                Mopt  = exp(-2*pi^2*sigma.^2.*u.^2.*(1/3600));
            elseif aperture == "disc"
                Mopt  = 1;
            end

            Ngt = obj.ganglion_density(csf_pars,'total');
            
            % Normalized receptive field size based on Watson's formula
            rf_size = (33163.2/2)./Ngt;

            % Temporal filter constants: time constants
            tau1      = tau10./(1+0.55.*log(1+((1+D).^0.6).*E./obj.par.t.*rf_size));
            tau2      = tau20./(1+0.37.*log(1+((1+D/3.2).^5).*E/120));

            % Modulation transfer function of the temporal response of the eye
            H1        = 1./(1+(2*pi.*w.*tau1).^2).^(n1/2); % cone responses
            H2        = 1./(1+(2*pi.*w.*tau2).^2).^(n2/2); % lateral inhibation

            % Modulation transfer function of the lateral inhibition process
            F         = 1-(1-exp(-(u./u0).^2)).^0.5; % high-pass filter characterizing lateral inhibation
            Mlat      = H1.*(1-H2.*F);

            % Spectral Density of photon noise
            PhiPhoton = 1./(eta.*p.*E);

            % Spatial equations
            X         = (X0.^-2+Xmax.^-2+(((0.5.*X0).^2+4*e.^2)./((0.5.*X0).^2)+e.^2).*(u.^2)./(Nmax.^2)).^-0.5;
            Y         = (Y0.^-2+Ymax.^-2+(((0.5.*X0).^2)./((0.5.*X0).^2)+e.^2).*(u.^2)./(Nmax.^2)).^-0.5;

            % Contrast Sensitivity Function
            S         =  sqrt(2) .* (Mopt./(2*k)).*sqrt((X.*Y.*T)./(PhiPhoton+Phi0./(Mlat.^2)));
        end

        function S = sensitivity_edge(obj, csf_pars)
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

        function Ng = ganglion_density(obj, csf_pars, cell_type)

            meridian  = csf_pars.vis_field;
            e         = csf_pars.eccentricity;

            % Find out if the stimulus is a gabor or disc
            dimension = size(meridian);

            % Density of total Retinal ganglion cells in fovea
            Ng0       = 33163.2;

            % Scale factor for decline in midget fraction with eccentricity 
            rm        = 41.03; 

            % Cell fraction - assumed that 50% of a cell type form an independent mosaic
            switch cell_type
                case "midget"
                    f0 = 0.8928/2 .* (1+e./rm).^-1;
                case "parasol"
                    % parasol/midget ratio (rough estimate from Dacey (1992) - is not included in watson's formula)
                    f0 = (1-0.8928)/2 .* (1+e./rm).^-1 .* (0.006.*e+0.015);
                case "total"
                    f0 = 1/2;
            end

            % Convert retinal meridians to visual field
            meridian  = rem(meridian+180,360);
            tem = meridian == 0;   sup = meridian == 90;
            nas = meridian == 180; inf = meridian == 270;

            if dimension(1,2) == 1 % Gabor
                
                % Parameters based on visual meridians for the Watson's formula
                temPar = tem.*[0.9851 1.058 22.14]; supPar = sup.*[0.9935 1.035 16.35];
                nasPar = nas.*[0.9729 1.084 7.633]; infPar = inf.*[0.9960 0.9932 12.13];
                visPar = temPar+supPar+nasPar+infPar;
    
                % Watson's formula for ganglion density
                a  = visPar(:,1); r2 = visPar(:,2); re = visPar(:,3);
                Ng = Ng0.*f0.*(a.*(1+(e./r2)).^-2+(1-a).*exp(-e./re));

            elseif dimension(1,1) == 1 % Disk
    
                % Parameters based on visual meridians for the Watson's formula
                temPar = tem.*[0.9851 1.058 22.14]'; supPar = sup.*[0.9935 1.035 16.35]';
                nasPar = nas.*[0.9729 1.084 7.633]'; infPar = inf.*[0.9960 0.9932 12.13]';
                visPar = temPar+supPar+nasPar+infPar;
    
                % Watson's formula for ganglion density
                a  = visPar(1,:); r2 = visPar(2,:); re = visPar(3,:);
                Ng = Ng0.*f0.*(a.*(1+(e./r2)).^-2+(1-a).*exp(-e./re));
            end
        end

    end

    methods( Static )

        function p = get_default_par()
            
            p         = CSF_base.get_dataset_par();

            p.k       = 7.13817;
            p.sigma0  = 0.384239;
            p.u0      = 2.45004;
            p.tau10   = 0.0366639;
            p.tau20   = 0.0917171;
            p.k_tem   = 0.53431;
            p.k_nas   = 0.20627;
            p.t       = 1.00901;

        end

    end
end

