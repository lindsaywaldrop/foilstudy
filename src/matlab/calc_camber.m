function camber = calc_camber(x, y)
    camber = (max(y) - min(y))./ (max(x) - min(x));
end