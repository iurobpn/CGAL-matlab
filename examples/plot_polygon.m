function plot_polygon(Y,c,a)

if nargin < 2
	c = [0,0,1];
end
if nargin < 3
	a = .15;
end

N = length(Y);

if iscell(Y)
	for i=1:N
		vs = Y{i};
		patch('XData',vs(:,1),'YData',vs(:,2),'FaceColor',c,'FaceAlpha',a)
		hold on
	end
else
	patch('XData',Y(:,1),'YData',Y(:,2),'FaceColor',c,'FaceAlpha',a)
end
	