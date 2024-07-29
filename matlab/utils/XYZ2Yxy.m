function [ Yxy ] = XYZ2Yxy( XYZ )
%XYZ2YXY covert CIE XYZ to Yxy. 
% input vector must have a length of 3.
% when transforming multiple values at once, the last dimension should have
% size 3

if ndims(XYZ) == 1
    if size(XYZ, 1) ~=3
        error('input vector must have 3 elements');
    end
    
    normTerm = sum(XYZ);
    Yxy = [XYZ(2), XYZ(1) / normTerm, XYZ(2) / normTerm];
    
elseif ndims(XYZ) == 2
    if size(XYZ, 2) ~=3
        error('input matrix must be of size Nx3');
    end
    
    normTerm = sum(XYZ, 2);
    Yxy = [XYZ(:, 2), XYZ(:, 1) ./ normTerm, XYZ(:, 2) ./ normTerm];
    
elseif ndims(XYZ) == 3
    if size(XYZ, 3) ~=3
        error('input matrix must be of size NxMx3');
    end
    
    normTerm = sum(XYZ, 3);
    Yxy = cat(3, XYZ(:, :, 2), XYZ(:, :, 1) ./ normTerm, XYZ(:, :, 2) ./ normTerm);
    
else
    error('xyz2Yxy does not support inputs with more than 3 dimensions')
end
    


end

