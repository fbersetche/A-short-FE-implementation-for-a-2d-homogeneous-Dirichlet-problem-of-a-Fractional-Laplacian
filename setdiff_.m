function e = setdiff_( A , B )
e = A;
b = B - A(1) + 1;
b( b<1 )=[];
e(b) = [];
end


