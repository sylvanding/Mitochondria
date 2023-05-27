function [map] = crop_sphere(imageSizeX,imageSizeY,imageSizeZ,centx,centy,centz,radius,radius2)
   [yy,xx,zz]=ndgrid((1:imageSizeX)-centx, (1:imageSizeY)-centy,(1:imageSizeZ)-centz);
   zfact=radius2/radius;
   map=(yy.^2+((xx.^2)/1+(zz.^2)/(zfact^2))<=radius^2);
end