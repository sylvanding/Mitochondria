# Mitochondria
This repository contains the analysis described in the paper: Quantifying nanoscopic alterations associated with mitochondrial dysfunction using three-dimensional single-molecule localization microscopy. 

The repository consists of 3 folders

Automated segmentation:

Instructions for automated segmentation:

1. Download the file 'temporal-Color_Code1' and place it in the FIJI/imagej folder in the file path 'fiji-win64\Fiji.app\plugins\Scripts\Image\Hyperstacks'

2. Open your 3D super-resolution image in Thunderstorm in FIJI

3. In THUNDERSTORM: Results, click "Plot histogram" and select the 'z' parameter

4. Select the region with the high signal and click "Apply ROI to filter"

5. In THUNDERSTORM: Results, click 'Apply'

6. In THUNDERSTORM: Results, click visualization

7. Make sure the '3D' option is checked. Click 'Auto size by results' and click 'OK'

8. Save the image stack as a tiff

9. In FIJI select 'Image'>'Hyperstack'>'Temporal-Color_Code1

10. Save the image as [your file name]_segmented.tif

11. Run the '[your file name]_segmented' through mitochondria_segmentation.ipynb and save the output numpy file to your working matlab folder

12. Download all .m files and save them to your working matlab folder

13. Download the tif and the zipped file. Extract the zipped file

14. Open 'make_boundary_auto.m' If you are using your own data, change line 1 to your 3d stack filename

15. run the code.

16. Open '[your file name]_segmented.tif in FIJI. If the code missed any boundaries between mitochondria, select the pencil tool and a line thickness of '2' and go over the image manually

17. Proceed to the mitochondrial analysis folder

Mitochondrial analysis:
Contains main script, "analyze_3d.m" and all dependent functions, as well as sample data. The goal is to extract morphological parameters from 3D mitochondria images. The sample data consists of
  "untreated_9.tif", which is a tiff stack of a fluorescent microscopy image captured using 3D SMLM with a pixel size of 25x25x25 nm
  and "untreated_9_segmented.tif" which is a 2D projection of the former image that has been automatically segmented using our neural network and then further segmented using the marker tool in in Fiji imageJ.
  
The code is ready to run as is. To analyze your own samples, you need to change the variables "filename" and "filename2" to the names of your 3d and 2d segmented tiff images, respectively. You may also need to change px to match your data's pixel size in nm. 
You may need to change CL_Dilationx to account for image distortions induced in the x-direction by the cylindrical lens. This can be measured by taking an image of a resolution chart or nanohole array and measuring any unilateral difference in size with and without the cylindrical lens in the detection path. If the image becomes smaller after insertion of the cylindrical lens in the x-direction, CL_Dilationx should be <1. It is calculated by S_CL/S_nCL, where S_CL is the distance between two features in the x-direction with the cylindrical lens, and S_nCL is the distance between the same features without the cylindrical lens.

You may also need to change CL_Dilationz to accound for image distortions induced in the z-direction by the refractive index mismatch. The method for calculating this factor is described in the supplementary materials of Huang B, Wang W, Bates M, Zhuang X. Three-dimensional super-resolution imaging by stochastic optical reconstruction microscopy. Science. 2008 Feb 8;319(5864):810-3. doi: 10.1126/science.1153529. Epub 2008 Jan 3. PMID: 18174397; PMCID: PMC2633023.

Finally, you may need to change the variables "thresh_fact" and "thresh_fact_slice". These variables are the factor by which the threshhold of the 2D image is multiplied and the factor by which threshhold of the 3D image is multiplied, respectively. In the course of running the analysis code, two images will be displayed. The first is called "Sum of masked frames", which shows the sum of the binzarized 3D image slices multiplied by the binarized 2D projection image. If this image excludes mitochondria, "thresh_fact" may be too high. If it includes background "thresh_fact" may be too low.

The second is called "Masked frame __" where __ is the central frame of the 3D data. It shows a slice of 3D data after binarization. This image should match with the corresponding frame of your tiff file. If there is background included, the variable "thresh_fact_slice" should be raised. If there is too little signal, the variable "thresh_fact_slice" should be lowered.

The code has a mat file output and three tiff stack outputs.

The mat contains the following variables:

vol: A list of volumes of each mitochondrion in cubic microns.

vol_ful: The volume of the full mitochondrial network in cubic microns.

lengths: A list of lengths of each mitochondrion in microns.

widths: A list of widths of each mitochondrion in microns.

blank_skel: A tiff stack showing a binarized mitochondrial skeleton

minsz: The user input value for minimum mitochondrion area in pixels

thresh: The user input value for threshhold.


The tiff outputs are:

The mitochondrial skeleton

A binarized image of the mitochondria

An 8-bit image of the mitochondria with background subtraction.

All three tiffs should be reviewed to ensure that the mitochondria were processed correctly.

Stitching:
Contains main script, "combnew4.m" and all dependent functions, as well as sample data. The goal is to create a single image from multiple overlapping images in cases where the mitochondrial network is larger than the imaging system field of view. As inputs, the script takes a list of filenames of csv files from analysis of SMLM datasets on thunderSTORM. The user must define the variable "changing" which is a corresponding list of the relative axial objective positions in each csv file. If all data was acquired at the same focus, "changing" should just be a list of zeros. 
As an output, the code saves a csv file made from combining all csv files.

The code also plots a graph of points from two stitched images. In the unlikely event that the stitching does not look correct, the user should put a break point at the line
"figure; imagesc(hey1)"

at the line
The user should then type the following commands in the command prompt:
"figure; imagesc(hey1)"
They should then find the point of maximum correlation, which should look like a symmetric local maximum.
They should then identify roughly the x- and y-coordinates of the maximum.

For example, if the x-coordinate of the maximum is between 200 and 280 and the y coordinate is between 1880 and 2020, they should enter the following in the command line:

TR1(:,3) = (TR1(:,3))-32*yshift;

TR1(:,4) = (TR1(:,4))-32*xshift;

[row,column] = find(hey1==max(max(hey1(1880:2020,200:280))));

yshift = row-ceil(sz1);

xshift = column-ceil(sz2);

TR1(:,3) = (TR1(:,3))+32*yshift;

TR1(:,4) = (TR1(:,4))+32*xshift;


figure; plot(TR0(1:round(length(x0)/60000):length(x0),3),TR0(1:(round(length(x0)/60000)):(length(x0)),4),'.','markersize',.01)

hold on; plot(TR1(1:60000,3),TR1(1:60000,4),'.','markersize',.01)

Simulation:
Contains main script, "call_mito.m" and all dependent functions. Simulates mitochondrial images with various parameters for the user to change.


