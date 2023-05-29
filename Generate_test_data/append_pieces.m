filename4='Untreated_2'; %Name of 3d tiff
cell_img=[];
for n=1:8
   
   filename=[filename4,"_" + n]; % Specify the output file name
   filename=append(filename(1),filename(2));
   tempsave=SPLMload([convertStringsToChars(filename),'.tif'],'tiff');
   cell_img=cat(3,cell_img,tempsave);
end
%%
filename3=[filename4,'.tif']; % Specify the output file name
[~,~,d]=size(cell_img);
for idx = 1:d
    %[A,map] = rgb2ind(temp2(:,:,:,idx),256);
    if idx == 1
        imwrite(uint16(cell_img(:,:,idx)),filename3);
    else
        imwrite(uint16(cell_img(:,:,idx)),filename3,'WriteMode','append');
    end
end