# Mitochondria
This document contains the analysis described in the paper: Quantifying nanoscopic alterations associated with mitochondrial dysfunctions using three-dimensional single-molecule localization microscopy. All code is written in MatLab

The code consists of 4 folders

Mitochondrial analysis:
Contains main script, "analyze_3d.m" and all dependent functions, as well as sample data. The sample data consists of
  "untreated_2.tif", which is a tiff stack of a fluorescent microscopy image captured using 3D SMLM with a pixel size of 25x25x25 nm
  and "untreated_2_segmented.tif" which is a 2D projection of the former image that has been manually segmented using the marker tool in Fiji imageJ.
  
The code is ready to run as is. To analyze your own samples, you need to change the variables "filename" and "filename2" to the names of your 3d and 2d segmented tiff images, respectively. You may also need to change px to match your data's pixel size in nm. 
You may need to change CL_Dilationx to account for image distortions induced in the x-direction by the cylindrical lens. This can be measured by taking an image of a resolution chart or nanohole array and measuring any unilateral difference in size with and without the cylindrical lens in the detection path. If the image becomes smaller after insertion of the cylindrical lens in the x-direction, CL_Dilationx should be <1. It is calculated by S_CL/S_nCL, where S_CL is the distance between two features in the x-direction with the cylindrical lens, and S_nCL is the distance between the same features without the cylindrical lens.

You may also need to change CL_Dilationz to accound for image distortions induced in the z-direction by the refractive index mismatch. The method for calculating this factor is described in the supplementary materials of Huang B, Wang W, Bates M, Zhuang X. Three-dimensional super-resolution imaging by stochastic optical reconstruction microscopy. Science. 2008 Feb 8;319(5864):810-3. doi: 10.1126/science.1153529. Epub 2008 Jan 3. PMID: 18174397; PMCID: PMC2633023.

Finally, you may need to change the variables "thresh_fact" and "thresh_fact_slice". These variables are the factor by which the threshhold of the 2D image is multiplied and the factor by which threshhold of the 3D image is multiplied, respectively. In the course of running the analysis code, two images will be displayed. The first is called "Sum of masked frames", which shows the sum of the binzarized 3D image slices multiplied by the binarized 2D projection image. If this image excludes mitochondria, "thresh_fact" may be too high. If it includes background "thresh_fact" may be too low.

The second is called "Masked frame __" where __ is the central frame of the 3D data. It shows a slice of 3D data after binarization. This image should match with the corresponding frame of your tiff file. If there is too much included, the variable "thresh_fact_slice" should be lowered. If there is too little, the variable "thresh_fact_slice" should be raised.



