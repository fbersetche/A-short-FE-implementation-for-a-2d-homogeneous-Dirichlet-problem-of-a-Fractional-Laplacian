function [x,y] = squareincircle(bs,s)
% Create  a circle centered at (0,0) with radius R using four segments.
% and a central square [-a a]x[-a a]
a=1;
R=1.6;
switch nargin
    case 0
        x = 8; % eigth edge segments
        return
    case 1
        A = [0,pi/2,pi,3*pi/2,a,a,-a,-a; % start parameter values
             pi/2,pi,3*pi/2,2*pi,-a,-a,a,a; % end parameter values
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
        cbs=find(bs==5); %segment 1
        x(cbs)=a;
        y(cbs)=s(cbs);
        cbs=find(bs==6); %segment 2
        y(cbs)=-a;
        x(cbs)=s(cbs);
        cbs=find(bs==7); %segment 3
        x(cbs)=-a;
        y(cbs)=s(cbs);
        cbs=find(bs==8); %segment 4
        y(cbs)=a;
        x(cbs)=s(cbs);
end
