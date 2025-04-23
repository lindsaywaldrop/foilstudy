function new_midline = adjust_midline(foil_midline, camber_new)
    new_midline = foil_midline;
    new_midline(:, 'y') = ((foil_midline(:, 'y') - min(foil_midline(:, 'y'))) .* camber_new)./...
        (max(foil_midline(:,'y')) - min(foil_midline(:,'y')));
    if min(new_midline{:, 'y'}) < 0
        new_midline(:,'y') = new_midline(:, 'y') + abs(min(new_midline(:, 'y')));
    end
    
end