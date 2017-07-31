function [x,y] = circleincircle(bs,s)
% Create  a circle centered at (0,0) with radius R using four segments.
% and a central circle with radius a
a=1;
R=1.1;
switch nargin
    case 0
        x = 8; % eigth edge segments
        return
    case 1
        A = [0,pi/2,pi,3*pi/2,0,-pi/2,-pi,-3*pi/2; % start parameter values
             pi/2,pi,3*pi/2,2*pi,-pi/2,-pi,-3*pi/2,-2*pi; % end parameter values
             1,1,1,1,1,1,1,1; % region label to left
             0,0,0,0,2,2,2,2]; % region label to right
        x = A(:,bs); % return requested columns
        return
    case 2
        x = zeros(size(s));
        y = zeros(size(s));
        if numel(bs) == 1 % Does bs need scalar expansion?
            bs = bs*ones(size(s)); % Expand bs
        end
        cbs=find(bs<=4); % outer circle
        x(cbs) = R*cos(s(cbs));
        y(cbs) = R*sin(s(cbs));
        cbs=find(bs>=5);
        x(cbs) = a*cos(s(cbs));
        y(cbs) = a*sin(s(cbs));
end
