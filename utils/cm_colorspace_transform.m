function out = cm_colorspace_transform( in, M )
% Transform an image or color vector (Nx3) to another color space using matrix M
% 
% out = cm_colorspace_transform( in, M )
%
% "in" could be an image [width x height x 3] or [rows x 3] matrix of colour
%     vectors, or [width x height x 3 x frames] colour video
% M is a colour transformation matrix so that
%     dest = M * src (where src is a column vector)
%

if any( size(M) ~= [3 3] )
    error( 'Colour transformation must be a 3x3 matrix' );
end

if length(size(in)) == 2 && size(in,2) == 3
    % Colour vector
    out = in * M';
elseif length(size(in)) == 3 && size(in,3) == 3
    % Color image
    out = reshape( (M * reshape( in, [size(in,1)*size(in,2) 3] )')', ...
        [size(in,1) size(in,2) 3] );
elseif length(size(in)) == 4 && size(in,3) == 3
    % Colour video    
    in_p = permute(in, [1 2 4 3] );
    out_p = reshape( (M * reshape( in_p, [size(in,1)*size(in,2)*size(in,4) 3] )')', ...
        [size(in,1) size(in,2) size(in,4) 3] );    
    out = permute( out_p, [1 2 4 3] );
else
    error( 'Unsupported colour representation' );
end

end