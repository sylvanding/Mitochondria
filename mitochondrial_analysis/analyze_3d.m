%The goal is to extract morphological parameters from 3D mitochondria images. The sample data consists of "untreated_2.tif", which is a tiff stack of a fluorescent microscopy image captured using 3D SMLM with a pixel size of 25x25x25 nm and "untreated_2_segmented.tif" which is a 2D projection of the former image that has been manually segmented using the marker tool in Fiji imageJ.

%The code is ready to run as is. To analyze your own samples, you need to change the variables "filename" and "filename2" to the names of your 3d and 2d segmented tiff images, respectively. You may also need to change px to match your data's pixel size in nm. You may need to change CL_Dilationx to account for image distortions induced in the x-direction by the cylindrical lens. This can be measured by taking an image of a resolution chart or nanohole array and measuring any unilateral difference in size with and without the cylindrical lens in the detection path. If the image becomes smaller after insertion of the cylindrical lens in the x-direction, CL_Dilationx should be <1. It is calculated by S_CL/S_nCL, where S_CL is the distance between two features in the x-direction with the cylindrical lens, and S_nCL is the distance between the same features without the cylindrical lens.
%You may also need to change CL_Dilationz to accound for image distortions induced in the z-direction by the refractive index mismatch. The method for calculating this factor is described in the supplementary materials of Huang B, Wang W, Bates M, Zhuang X. Three-dimensional super-resolution imaging by stochastic optical reconstruction microscopy. Science. 2008 Feb 8;319(5864):810-3. doi: 10.1126/science.1153529. Epub 2008 Jan 3. PMID: 18174397; PMCID: PMC2633023.
%Finally, you may need to change the variables "thresh_fact" and "thresh_fact_slice". These variables are the factor by which the threshhold of the 2D image is multiplied and the factor by which threshhold of the 3D image is multiplied, respectively. In the course of running the analysis code, two images will be displayed. The first is called "Sum of masked frames", which shows the sum of the binzarized 3D image slices multiplied by the binarized 2D projection image. If this image excludes mitochondria, "thresh_fact" may be too high. If it includes background "thresh_fact" may be too low.
%The second is called "Masked frame __" where __ is the central frame of the 3D data. It shows a slice of 3D data after binarization. This image should match with the corresponding frame of your tiff file. If there is too much included, the variable "thresh_fact_slice" should be lowered. If there is too little, the variable "thresh_fact_slice" should be raised.
%The code has a mat file output and three tiff stack outputs.
%The mat contains the following variables vol: A list of volumes of each mitochondrion in cubic microns. vol_ful: The volume of the full mitochondrial network in cubic microns. lengths: A list of lengths of each mitochondrion in microns. widths: A list of widths of each mitochondrion in microns. blank_skel: A tiff stack showing a binarized mitochondrial skeleton minsz: The user input value for minimum mitochondrioon area in pixels thresh: The user input value for threshhold.
%The tiff outputs are: The mitochondrial skeleton A binarized image of the mitochondria An 8-bit image of the mitochondria with background subtraction.
%All three tiffs should be reviewed to ensure that the mitochondria were processed correctly.

clear;
filename='Untreated_2'; %Name of 3d tiff
filename2='Untreated_2_segmeted'; %Name of 2d segmented tiff
px=.025; %pixel size in um
thresh_fact=1.5; %Factor by which the threshhold of the 2D image is multiplied (lower if "Sum of masked frames" excludes mitochondria. Raise if it includes background)
thresh_fact_slice=1/8; %Factor by which threshhold of the 3D image is multiplied (lower if "Masked frame __" excludes mitochondria. Raise if it includes background)
CL_Dilationx=.8743; %the dilation factor in the x-direction 
CL_Dilationz=1.3889; %The dilation factor in the z-direction (Huang et. al)
cell_img=SPLMload([filename,'.tif'],'tiff');
[y,x,z]=size(cell_img);
[Y,X,Z] = (meshgrid(1:1:x,1:1:y,1:1:z));
[Yq,Xq,Zq] = meshgrid(.8743:.8743:x,1:y,1.3889:1.3889:z);
cell_img=interp3(Y,X,Z,double(cell_img),Yq,Xq,Zq);
cell_img=uint8(cell_img);
cell_img2=SPLMload([filename2,'.tif'],'tiff');
[y2,x2]=size(cell_img2);
[Y2,X2] = meshgrid(1:1:x2,1:1:y2);
[Yq2,Xq2] = meshgrid(.8743:.8743:x2,1:y2);
cell_img2=interp2(Y2,X2,double(cell_img2),Yq2,Xq2);
[cell_img,cell_img2]=im_shift((cell_img),cell_img2);
%%
cell_img2=uint8(cell_img2);

thresh=((mean(mean(cell_img2(cell_img2>0)))))/255;
thresh=thresh*thresh_fact;
minsz=200;
cell_img2bin = imbinarize(cell_img2(:,:), thresh);
cell_img2bin = bwareafilt(cell_img2bin(:,:),[minsz 500000000]);
cell_img2bin = ~bwareaopen(~cell_img2bin, 2000);

thresh2=thresh*thresh_fact_slice;

%cell_img2=imgaussfilt(cell_img,2);
%mask_generate_temp = (cell_img2>254);
for n=1:size(cell_img,3)
  
    ce_temp=cell_img(:,:,n);
    mean(mean(ce_temp(ce_temp>0)))
   % mask_generate_temp(:,:,n) = imbinarize(imgaussfilt(cell_img(:,:,n),1), 0.00344803563250706*double(mean(mean(ce_temp(ce_temp>0)))));
    mask_generate_temp(:,:,n) = imbinarize(imgaussfilt(cell_img(:,:,n),1), thresh2);
    mask_generate_full(:,:,n)=mask_generate_temp(:,:,n);
    mask_generate_temp(:,:,n)=mask_generate_temp(:,:,n).*cell_img2bin;
    %mask_generate_temp(:,:,n)= imbinarize(cell_img(:,:,n), 'adaptive','Sensitivity', 0.25);
    %mask_generate_temp = imbinarize(temp_mask, 0.05);
    % se = strel('disk',3);
    % mask_generate = imerode(mask_generate_temp, se);
    % figure
    mask_generate(:,:,n)=mask_generate_temp(:,:,n);
    if n==9
        ffdjgnf=4
    end
    mask_generate(:,:,n) = bwareafilt(mask_generate(:,:,n),[minsz 500000000]);
    mask_generate_full(:,:,n)=bwareafilt(mask_generate_full(:,:,n),[minsz 500000000]);
    %mask_generate = imfill(mask_generate, 'holes');
     mask_generate2(:,:,n) = ~bwareaopen(~mask_generate(:,:,n), 2000);
    mask_generate_full(:,:,n) = ~bwareaopen(~mask_generate_full(:,:,n), 2000);

end
frameshow=round(size(mask_generate,3)/2);
figure; imshow(mask_generate2(:,:,frameshow))
title("masked frame " + round(frameshow*CL_Dilationz))
cell_img_mask1=uint8(mask_generate).*cell_img;
cell_img_mask2=uint8(mask_generate2).*cell_img;


figure; imagesc(sum(mask_generate2,3))
title('sum of masked frames')
%%

%imagesc(mask_generate>0)
con6=(mask_generate2>0);
con7=(sum(mask_generate2,3)>0);
%boundaries = I>0 & L == 0;




B = bwboundaries(con7);
% B3=[];
% for n=1:size(con6,3)
%     B3=[B3; bwboundaries(con6(:,:,n))];
% 
% end
figure; imshow(sum(cell_img_mask1,3),[0,900]);
colormap gray
hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1)
   area(k)=length(B{k,1});
end

% figure; imshow(sum(cell_img_mask1,3),[0,900]);
% colormap gray
% hold on
% for k = 1:length(B3)
%    boundary = B3{k};
%    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1)
%    area(k)=length(B3{k,1});
% end
for n=1:size(cell_img,3)
   cell_img_3d2(:,:,n)=cell_img(:,:,n).* uint8(mask_generate2(:,:,n));
end
%%
objects=bwconncomp(mask_generate2);
for n=1:length(objects.PixelIdxList)
   volume(n)=px*px*px*length(objects.PixelIdxList{1,n}); %volume in microns 
    vol_full=px*px*px*length(find(mask_generate_full(mask_generate_full>0)));
end

I2=con7;
objects2=bwconncomp(I2);
for n=1:length(objects2.PixelIdxList)        
list=objects2.PixelIdxList{1,n};
[row,col]=ind2sub(size(I2),list);
dist=rand(1,1,3);
dist=dist/norm(squeeze(dist));
for n2=1:length(row)
    blank(row(n2),col(n2),1:3)=single(I2(row(n2),col(n2))).*dist;
end
end


[x,y,z]=size(cell_img);
blank=zeros(x,y,z);
blank_skel=blank;
for n=1:objects.NumObjects
    blank=zeros(x,y,z);
    n
    blank(objects.PixelIdxList{1,n})=1;
    indices=find(blank);
    [x2,y2,z2]=ind2sub([x,y,z],indices);
   
    if length(unique(z2))<6
       continue 
    end
    for n2=1:size(blank,3)
       blank(:,:,n2)=imfill(blank(:,:,n2));
    end
    [skel_fin,crackWidthImage,length_skel]=shortpathpls3d(blank(min(x2):max(x2),min(y2):max(y2),min(z2):max(z2)));
    blank_skel(min(x2):max(x2),min(y2):max(y2),min(z2):max(z2))=blank_skel(min(x2):max(x2),min(y2):max(y2),min(z2):max(z2))+skel_fin;
    widths(n)=mean(crackWidthImage(~isnan(crackWidthImage)))*px;
    lengths(n)=length_skel*px;
    vol(n)=sum(sum(sum(blank))).*(px^3);

 
end
lengths=(nonzeros(lengths));
vol=(nonzeros(vol));
widths=(nonzeros(widths));
blank_skel=logical(blank_skel);
%%
filename35=[filename,'_morph'];
save(filename35,'vol','vol_full','lengths','widths','blank_skel','minsz','thresh')
%%
%prompt='What do you want to call the output gif?:\n';
%filename2 = input(prompt,'s');




filename3=[filename,'_binarized','.tif']; % Specify the output file name
[~,~,d]=size(cell_img);
for idx = 1:d
    %[A,map] = rgb2ind(temp2(:,:,:,idx),256);
    if idx == 1
        imwrite(uint8(con6(:,:,idx)),filename3);
    else
        imwrite(uint8(con6(:,:,idx)),filename3,'WriteMode','append');
    end
end


filename4=[filename,'_skeletonized','.tif']; % Specify the output file name
[~,~,d]=size(cell_img);
for idx = 1:d
    %[A,map] = rgb2ind(temp2(:,:,:,idx),256);
    if idx == 1
        imwrite(uint8(blank_skel(:,:,idx)),filename4);
    else
        imwrite(uint8(blank_skel(:,:,idx)),filename4,'WriteMode','append');
    end
end
filename4=[filename,'_no_background','.tif']; % Specify the output file name
[~,~,d]=size(cell_img);
for idx = 1:d
    %[A,map] = rgb2ind(temp2(:,:,:,idx),256);
    if idx == 1
        imwrite(uint8(cell_img_mask2(:,:,idx)),filename4);
    else
        imwrite(uint8(cell_img_mask2(:,:,idx)),filename4,'WriteMode','append');
    end
end
