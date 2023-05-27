function [map] = crop_circ3d(imageSizeX,imageSizeY,imageSizeZ,centx,centy,centz,slopexy,slopexz,slopeyz)
   [xx,yy,zz]=ndgrid((1:imageSizeX), (1:imageSizeY),(1:imageSizeZ));
   map=round(xx*slopeyz+yy*slopexz+zz*slopexy,6,'significant')==round(centx*slopeyz+centy*slopexz+centz*slopexy,6,'significant');
end