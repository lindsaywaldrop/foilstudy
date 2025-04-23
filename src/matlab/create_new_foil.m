function new_foil = create_new_foil(new_midline)

npts = size(new_midline,1)*2;
half_thickness = sqrt((new_midline{:, 'thickness_top_x'}).^2 + ...
    (new_midline{:, 'thickness_top_y'}).^2);
top_angles = calc_midline_angles(new_midline) - 0.5*pi;
bot_angles = calc_midline_angles(new_midline) + 0.5*pi;
new_foil_x = zeros(npts,1);
new_foil_y = zeros(npts,1);
new_foil_x(1:(npts/2)) = new_midline{:,'x'} - half_thickness.*cos(top_angles);
new_foil_y(1:(npts/2)) = new_midline{:,'y'} - half_thickness.*sin(top_angles);
new_foil_x((npts/2 + 1):npts) = flip(new_midline{:, 'x'} - half_thickness.*cos(bot_angles));
new_foil_y((npts/2 + 1):npts) = flip(new_midline{:, 'y'} - half_thickness.*sin(bot_angles));
new_foil_x(1) = 1;
new_foil_y(1) = 0;
new_foil_x(npts) = new_foil_x(1) - 5e-4;
new_foil_y(npts) = 0;
new_foil = table(new_foil_x, new_foil_y, 'VariableName', {'x', 'y'});
end

function angles = calc_midline_angles(new_midline)

angles = zeros(size(new_midline{:,'x'}));
for i = 2:length(angles)
    p1 = [new_midline{i, 'x'}, new_midline{i, 'y'}];
    p2 = [new_midline{i-1, 'x'}, new_midline{i-1, 'y'}];
    p3 = [new_midline{i, 'x'}, new_midline{i-1, 'y'}];
    angles(i) = calc_angle(p1, p2, p3);
end
angles(1) = angles(2);

end

function alpha = calc_angle(p1, p2, p3)

a = [p3(1) - p2(1), p3(2) - p2(2)];
b = [p1(1) - p2(1), p1(2) - p2(2)];

alpha = acos((a(1).*b(1) + a(2).*b(2))./(mag(a).* mag(b)));

end

function m = mag(x)

m = sqrt((x(1))^2 + (x(2))^2);

end