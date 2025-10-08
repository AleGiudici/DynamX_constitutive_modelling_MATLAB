function [output_radius] = radiusFromIincompressibility(input_radius,axial_length,vessel_volume,choice)
% This function uses the incompressibility assumption to calculate the
% inner radius from outer radius (or viceversa). Input parameters are the
% input radius (either outer and inner), the vessel deformed axial length,
% the vessel volume and "choice" which can take values:
%   1 --> determines the inner radius given an input outer radius, or
%   2 --> determines the outer radius given an input innner radius.

    switch(choice)
        case 1 % inner radius from outer radius
            output_radius = sqrt(input_radius.^2-(vessel_volume./(pi*axial_length))); % inner radii in loaded configurations (mm)
        case 2 % outer radius from inner radius
            output_radius = sqrt(input_radius.^2+(vessel_volume./(pi*axial_length))); % inner radii in loaded configurations (mm)
    end
end