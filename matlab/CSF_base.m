classdef CSF_base
    %Super class of all Spatio-chromatic CSF models
    
    properties
        par; % Model parameters as a structure
        mean_error;
        mean_std;
    end
    
    methods( Abstract )
        
        % A short name that could be used as a part of a file name
        name = short_name( obj )
                

        % A general interface to compute sensitivity, which takes as input a
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
        % All parameters can be specified using multidemensional arrays as
        % long as of they have compatible (broadcastable) sizes. For example 
        % if lum has the size [10 1] and s_freq has the size [1 5], calling:
        %
        % csf_pars = struct( 'luminance', lum, 's_frequency', s_freq, 'area', 1 );          
        % S = csf_model.sensitivity( csf_pars );        
        % 
        % will return S of the size [10 5] computed for all combinations of
        % lum and s_freq. 
        %
        % 'lms_bkg' and 'lms_delta' must have their last dimension of size 3. 
        %
        % For best performance, pass vectors with the the required number
        % of parameters. Avoid calling the sensitivity() function in a loop.      
        %
        % You must specify 'luminance' or 'lms_bkg' but not both. 
        % You must specify 'area' or 'ge_sigma' but not both. 
        %
        % Example: 
        %
        % csf_model = <model_name>()
        % csf_pars = struct( 's_frequency', 4, 't_frequency', 1, 'orientation', 0, 'lms_bkg', [0.7443 0.3054 0.0157], 'area', 1, 'eccentricity', 0 );          
        % S = csf_model.sensitivity( csf_pars );
        S = sensitivity( obj, pars );
                
    end
    
    methods
        
        function str = full_name( obj )
            str = obj.short_name();
        end
        
        % Internally used for training CSF
        function obj = set_pars( obj, pars_vector )
            % Set the parameters of the model, supplied as a row vector
            % (used for optimizing the parameters of the model)
            
            assert( ~isempty( obj.par ) ); % obj.par must be initialized before calling this function
            
            obj.par = obj.param2struct( obj.par, pars_vector );
        end
        
        % Internally used for training CSF
        function pars_vector = get_pars( obj )
            % Get the parameters of the model as a row vector
            % (used for optimizing the parameters of the model)
            
            pars_vector = obj.struct2param( obj.par );
        end
        
        % Predict the sensitivity for a detection of a Gabour patch of certain chromatic
        % direction and amplitide. 
        %
        % S = sensitivity_stolms( obj, s_freq, t_freq, orientation, LMS_bkg, LMS_delta, area, eccentricity );
        %
        % Important: all parameters must be column vectors. LMS_bkg and
        %            LMS_delta are Nx3 matrices
        %
        % freq - spatial frequency in cpd
        % LMS_bkg - LMS of the background colour (CIE2006 CMF)
        % LMS_delta - colour direction vector in the LMS space (LMS_peak-LMS_mean)
        % area - area in deg^2
        %
        % The method returns:
        % S - Sensitivity (the inverse of cone contrast at the threshold)
        function S = sensitivity_stolms( obj, s_freq, t_freq, orientation, LMS_bkg, LMS_delta, area, eccentricity )
            csf_pars = struct( 's_frequency', s_freq, 't_frequency', t_freq, 'orientation', orientation, 'lms_bkg', LMS_bkg, 'lms_delta', LMS_delta, 'area', area, 'eccentricity', eccentricity );
            S = obj.sensitivity( csf_pars );
        end  

        function S = sensitivity_stolms_jov( obj, s_freq, t_freq, orientation, LMS_bkg, LMS_delta, area, eccentricity, col_dir )
            csf_pars = struct( 's_frequency', s_freq, 't_frequency', t_freq, 'orientation', orientation, 'lms_bkg', LMS_bkg, 'lms_delta', LMS_delta, 'area', area, 'eccentricity', eccentricity );
            if strcmp(class(obj.csf_model), 'CSF_JOV')
                S = obj.csf_model.sensitivity( csf_pars, col_dir );
            else
                S = obj.csf_model.sensitivity( csf_pars );
            end
        end

        function S = sensitivity_stolmsv( obj, s_freq, t_freq, orientation, LMS_bkg, LMS_delta, area, eccentricity, vis_field )
            csf_pars = struct( 's_frequency', s_freq, 't_frequency', t_freq, 'orientation', orientation, 'lms_bkg', LMS_bkg, 'lms_delta', LMS_delta, 'area', area, 'eccentricity', eccentricity, 'vis_field', vis_field );
            S = obj.sensitivity( csf_pars );
        end        
        
        function S = sensitivity_stolms_edge( obj, t_freq, orientation, LMS_bkg, LMS_delta, ge_sigma, eccentricity )
            csf_pars = struct( 't_frequency', t_freq, 'orientation', orientation, 'lms_bkg', LMS_bkg, 'lms_delta', LMS_delta, 'ge_sigma', ge_sigma, 'eccentricity', eccentricity);
            S = obj.sensitivity_edge( csf_pars );
        end 
        
        % Test whether all the parameters are correct size, that the
        % names are correct, set the default values for the missing
        % parameters. 
        % pars - the csf_params structure
        % requires - a cell array with the selected parameters that are
        %            required. Used for the multually exclusive
        %            parameters, such as 'luminance'/'lms_bkg' and 
        %            'area'/'ge_sigma' so that one of them is computed as
        %            needed, but not both.
        % expand - if true, all the parameters are expanded to have the
        %          same size.
        function pars = test_complete_params(obj, pars, requires, expand )

            if ~exist( 'expand', 'var' )
                expand = false;
            end

            valid_names = { 'luminance', 'lms_bkg', 'lms_delta', 's_frequency', 't_frequency', 'orientation', 'area', 'ge_sigma', 'eccentricity', 'vis_field' };            
            fn = fieldnames( pars );
%             N = 1; % The size of the vector
            cur_par = 1;
            for kk=1:length(fn)
                if ~ismember( fn{kk}, valid_names )
                    error( 'Parameter structure contains unrecognized field ''%s''', fn{kk} );
                end

                % Check whether the parameters can be broadcasted                
                try 
                    param = pars.(fn{kk});
                    if ismember( fn{kk}, { 'lms_bkg', 'lms_delta' } )
                        p_sz = size(param);
                        if p_sz(end) ~= 3
                            error( 'The last dimension of ''%s'' must have size 3', fn{kk} );
                        end
                        %param = reshape( param, [p_sz(1:(end-1)) 1 3] );
                    end                    
                    cur_par = cur_par .* param;
                catch
                    error( 'Parameter %s cannot be broadcasted', fn{kk});
                end
                
%                 if numel(pars.(fn{kk})) > 1
%                     Nc = numel(pars.(fn{kk}))/par_len;
%                     if N==1
%                         N = Nc;
%                     else
%                         if Nc~=1 && N ~= Nc
%                             error( 'Inconsistent size of the parameter ''%s''', fn{kk} );
%                         end
%                     end
%                 end
            end

            if ismember( 'luminance', requires )
                if ~isfield( pars, 'luminance')
                    if ~isfield( pars, 'lms_bkg')
                        error( 'You need to pass either luminance or lms_bkg parameter.')
                    end
                    pars.luminance = CSF_base.last_dim(pars.lms_bkg,1) + CSF_base.last_dim(pars.lms_bkg,2);
                end
            end

            if ismember( 'lms_bkg', requires )
                if ~isfield( pars, 'lms_bkg')
                    if ~isfield( pars, 'luminance')
                        error( 'You need to pass either luminance or lms_bkg parameter.')
                    end
                    %error( 'Not implemented' )
                    pars.lms_bkg = [0.6991 0.3009 0.0198] .* pars.luminance;
                end
            end

            if ismember( 'ge_sigma', requires )
                if ~isfield( pars, 'ge_sigma')
                    if ~isfield( pars, 'area')
                        error( 'You need to pass either ge_sigma or area parameter.')
                    end
                    pars.ge_sigma = sqrt(pars.area/pi);
                end
            end

            if ismember( 'area', requires )
                if ~isfield( pars, 'area')
                    if ~isfield( pars, 'ge_sigma')
                        error( 'You need to pass either ge_sigma or area parameter.')
                    end
                    pars.area = pi*pars.ge_sigma.^2;
                end
            end
            
            % Default parameter values
            def_pars = struct( 'eccentricity', 0, 'vis_field', 180, 'orientation', 0, 't_frequency', 0, 'lms_delta', [0.6855 0.2951 0.0194] );
            fn_dp = fieldnames( def_pars );
            for kk=1:length(fn_dp)
                if ~isfield(pars, fn_dp{kk})
                    pars.(fn_dp{kk}) = def_pars.(fn_dp{kk});
                end
            end
            
%             if expand && N>1
%                 % Make all parameters the same height 
%                 fn = fieldnames( pars );
%                 for kk=1:length(fn)
%                     if size(pars.(fn{kk}),1)==1
%                         pars.(fn{kk}) = repmat( pars.(fn{kk}), [N 1]);
%                     end
%                 end                
%             end            

        end
        
                
       function print( obj, fh )
            % Print the model parameters in a format ready to be pasted into
            % get_default_par()
            
            obj.print_struct( fh, 'p.', obj.par );            
        end
         
        function print_struct( obj, fh, struct_name, s )
            % Print the model parameters in a format ready to be pasted into
            % get_default_par()
            
            fn = fieldnames( s );
            for ff=1:length(fn)
                if ismember( fn{ff}, { 'cm', 'ds', 'sust', 'trans' } )
                    continue;
                end
                if isstruct(s.(fn{ff}))
                    obj.print_struct( fh, strcat( struct_name, fn{ff}, '.' ), s.(fn{ff}) );
                else
                    fprintf( fh, '\t%s%s = ', struct_name, fn{ff} );
                    obj.print_vector( fh, s.(fn{ff}) );
                    fprintf( fh, ';\n' );
                end
            end
            
        end
                
    end
    
    methods(Static)


        function to = update_struct( from, to )
            fn = fieldnames(from);
            for ff = 1:length(fn)
               if isfield(to, fn{ff})                   
                   if isstruct(to.(fn{ff}))
                       to.(fn{ff}) = CSF_base.update_struct( from.(fn{ff}), to.(fn{ff}) );
                   else
                       assert( all( size(to.(fn{ff})) == size(from.(fn{ff})) ) );
                       to.(fn{ff}) = from.(fn{ff});
                   end
               end
           end            

        end
        
        function Y = sel_dim(X,d)
            cln(1:ndims(X)) = {1};
            if d>1
                cln(d) = {':'};
            else
                cln(end-d+1) = {':'};
            end
            Y = X(cln{:});
        end

        function Y = last_dim(X,d)
            cln(1:ndims(X)) = { ':' };
            cln{end} = d;
            Y = X(cln{:});
        end

        function [s, pos] = param2struct( s, pars_vector )
            % Convert a vector with the parameters to a structure with the
            % same fields as in the structure "s".
            
            pos = 1;
            for cc=1:length(s)
                ff = fieldnames(s(cc));
                for kk=1:length(ff)
                    
                    if isstruct(s(cc).(ff{kk}))
                        [s(cc).(ff{kk}), pos_ret] = CSF_base.param2struct(s(cc).(ff{kk}), pars_vector(pos:end) );
                        pos = pos+pos_ret-1;
                    else
                        N = length(s(cc).(ff{kk}));
                        s(cc).(ff{kk}) = pars_vector(pos:(pos+N-1));
                        pos = pos + N;
                    end
                end
            end
%             if (pos-1) ~= length(pars_vector)
%                 error( 'The parameter vector contains %d elements while the model has %d optimized parameters. Perhaps the optimized for a different set of datasets?', length(pars_vector), (pos-1) );
%             end
        end
        
        
        function pars_vector = struct2param( s )
            % Convert a structure to a vector with the parameters
            
            pars_vector = [];
            
            for cc=1:length(s)
                ff = fieldnames(s(cc));
                for kk=1:length(ff)
                    
                    if isstruct( s(cc).(ff{kk}) )
                        pars_vector = cat( 2, pars_vector, CSF_base.struct2param(s(cc).(ff{kk})) );
                    else
                        pars_vector = cat( 2, pars_vector, s(cc).(ff{kk}) );
                    end
                    
                end
            end
            
        end
        
        
        function v = get_lum_dep( pars, L )
            % A family of functions modeling luminance dependency
            
            %log_lum = log10(L);
            
            switch length(pars)
                case 1
                    % Constant
                    v = ones(size(L)) * pars(1);
                case 2
                    % Linear in log
                    v = pars(2)*L.^pars(1);
                    %v = 10.^(pars(1)*log_lum + log10(pars(2)));
                case 3
                    % Log parabola
                    %        v = pars(1) * 10.^(exp( -(log_lum-pars(2)).^2/pars(3) ));
                    
                    % A single hyperbolic function
                    v = pars(1)*(1+pars(2)./L).^(-pars(3));
                case 5
                    % Two hyperbolic functions
                    v = pars(1)*(1+pars(2)./L).^(-pars(3)) .* (1-(1+pars(4)./L).^(-pars(5)));
                otherwise
                    error( 'not implemented' );
            end
            
        end

        function v = get_lum_dep_dec( pars, L )
            % A family of functions modeling luminance dependency. 
            % The same as abobe but the functions are decreasing with
            % luminance
            
            log_lum = log10(L);
            
            switch length(pars)
                case 1
                    % Constant
                    v = ones(size(L)) * pars(1);
                case 2
                    % Linear in log
                    v = 10.^(-pars(1)*log_lum + pars(2));
                case 3        
                    % A single hyperbolic function
                    v = pars(1)* (1-(1+pars(2)./L).^(-pars(3)));
                case 5
                    % Two hyperbolic functions
                    error( 'TODO' );
                    v = pars(1)*(1+pars(2)./L).^(-pars(3)) .* (1-(1+pars(4)./L).^(-pars(5)));
                otherwise
                    error( 'not implemented' );
            end
            
        end
        
        function p = get_dataset_par()
            p = struct();            
        end
        
        function print_vector( fh, vec )
            if length(vec)>1
                fprintf( fh, '[ ' );
                fprintf( fh, '%g ', vec );
                fprintf( fh, ']' );
            else
                fprintf( fh, '%g', vec );
            end
        end        
        
    end
end

