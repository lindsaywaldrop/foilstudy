function foil_midline = find_midline(X,Y,smoothed)
npts = length(X);
foil_midline_x = zeros(npts/2,1);
foil_midline_y = zeros(npts/2,1);
thickness_top_x = zeros(npts/2,1);
thickness_top_y = zeros(npts/2,1);
thickness_bot_x = zeros(npts/2,1);
thickness_bot_y = zeros(npts/2,1);

for i = 1:(npts/2)
    foil_midline_x(i) = (X(i)- X((npts + 1) - i))/2 + X((npts + 1) - i);
    foil_midline_y(i) = (Y(i)- Y((npts+1)-i))/2 + Y((npts + 1) - i);
end

for i = 1:(npts/2)
    thickness_top_x(i) = X(i) - foil_midline_x(i);
    thickness_top_y(i) = Y(i) - foil_midline_y(i);
    thickness_bot_x(i) = foil_midline_x(i) - X((npts + 1) - i);
    thickness_bot_y(i) = foil_midline_y(i) - Y((npts + 1) - i);
end

if smoothed==1

end

foil_midline = table(foil_midline_x, foil_midline_y, thickness_top_x, ...
    thickness_top_y, thickness_bot_x, thickness_bot_y, ...
    'VariableNames', {'x', 'y', 'thickness_top_x', 'thickness_top_y',...
    'thickness_bot_x', 'thickness_bot_y'});

end



