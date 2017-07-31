function ML = comp_quad(Bl, x0, y0, s , phi , R, areal , p_I , w_I , p_T)
x = (Bl*p_T')' + [x0.*ones(length(p_T),1) , y0.*ones(length(p_T),1)];
aux = x(:,1)*cos(2*pi*p_I') + x(:,2)*sin(2*pi*p_I');
weight = ( ( -aux + sqrt( aux.^2 + R^2 - ( x(:,1).^2 + x(:,2).^2 )*ones(1,length(p_I)) ) ).^(-2*s) )*w_I;
ML = (areal*2*pi/s).*reshape( phi*weight , 3 , 3);
end
