function err = nmse(xr, xorg)

% calculates NMSE as sum((xr - xorg)^2)/sum(xorg^2)

err = sum((xr(:) - xorg(:)).^2)/sum(xorg(:).^2);