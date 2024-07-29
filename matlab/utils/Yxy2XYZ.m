function XYZ = Yxy2XYZ( Yxy, x, y )

sz = size(Yxy);

if( length(sz) == 2 && sz(2) == 3 ) % In a vector format
    
    XYZ = cat( 2, Yxy(:,1).*Yxy(:,2)./Yxy(:,3), Yxy(:,1), ...
        Yxy(:,1)./Yxy(:,3).*(1-Yxy(:,2)-Yxy(:,3)) );
    
else % in an image or other format        
    
    if( length( Yxy ) == 1 && ~exist( 'y', 'var' ) )
        % use D65 if only one param
        Yxy = [Yxy 0.31271 0.32902];        
        XYZ = Yxy(1) * [ Yxy(2)/Yxy(3) 1 1/Yxy(3)*(1-Yxy(2)-Yxy(3)) ];
        return;
    end
    
    if( exist( 'y', 'var' ) )
        Yxy = cat( 3, Yxy, repmat( x, sz ), repmat( y, sz ) );
    elseif( length( sz ) == 2 )
        % use D65 if only one param
        Yxy = cat( 3, Yxy, repmat( 0.31271, sz ), repmat( 0.32902, sz ) );
    end
    
    XYZ = cat( 3, Yxy(:,:,1).*Yxy(:,:,2)./Yxy(:,:,3), Yxy(:,:,1),  Yxy(:,:,1)./Yxy(:,:,3).*(1-Yxy(:,:,2)-Yxy(:,:,3)) );
    
end

end
