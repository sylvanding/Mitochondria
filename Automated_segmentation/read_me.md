To run this code:
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
13. Download the tif and the zipped file. Extract them
14. Open 'make_boundary_auto.m' If you are using your own data, change line 1 to your 3d stack filename
15. run the code.
16. Open '[your file name]_segmented.tif in FIJI. If the code missed any boundaries between mitochondria, select the pencil tool and a line thickness of '2' and go over the image manually
17. Proceed to the mitochondrial analysis folder

