function VL = fquad( areal, xl , yl , f )
VL = zeros(3,1);
xmid = [(xl(2)+xl(3))/2, (xl(1)+xl(3))/2, (xl(1)+xl(2))/2];
ymid = [(yl(2)+yl(3))/2, (yl(1)+yl(3))/2, (yl(1)+yl(2))/2];
for i=1:3
    for j=1:3
        if j~=i
            VL(i) = VL(i) + areal/6 * f(xmid(j), ymid(j));
        end
    end
end
end

