%The goal is to create a single image from multiple overlapping images in cases where the mitochondrial network is larger than the imaging system field of view. As inputs, the script takes a list of filenames of csv files from analysis of SMLM datasets on thunderSTORM. The user must define the variable "changing" which is a corresponding list of the relative axial objective positions in each csv file. If all data was acquired at the same focus, "changing" should just be a list of zeros. As an output, the code saves a csv file made from combining all csv files.
%The code also plots a graph of points from two stitched images. If the stitching does not look correct, the user should put a break point at the line "figure; imagesc(hey1)"
%at the line The user should then type the following commands in the command prompt: "figure; imagesc(hey1)" They should then find the point of maximum correlation, which should look like a symmetric local maximum. They should then identify roughly the x- and y-coordinates of the maximum.
%For example, if the x-coordinate of the maximum is between 200 and 280 and the y coordinate is between 1880 and 2020, they should enter the following in the command line (without the % symbols):
%TR1(:,3) = (TR1(:,3))-32*yshift;
%TR1(:,4) = (TR1(:,4))-32*xshift;
%[row,column] = find(hey1==max(max(hey1(1880:2020,200:280))));
%yshift = row-ceil(sz1);
%xshift = column-ceil(sz2);
%TR1(:,3) = (TR1(:,3))+32*yshift;
%TR1(:,4) = (TR1(:,4))+32*xshift;
%figure; plot(TR0(1:round(length(x0)/60000):length(x0),3),TR0(1:(round(length(x0)/60000)):(length(x0)),4),'.','markersize',.01)
%hold on; plot(TR1(1:60000,3),TR1(1:60000,4),'.','markersize',.01)

fnam=[{'untreated_9_min0_1'};{'untreated_9_min500_1'};{'untreated_9_min0_2'};{'untreated_9_min500_2'};{'untreated_9_min500_3'}]
changing=[0,-500,0,-500,-500];
savefile = 'untreated_9_fin.csv';
for n2=1:(length(fnam)-1)
    fact=n2;
    add=changing(n2);
    if n2==1
    fnam0=fnam{1};
    fnam1=fnam{2};
    
    fname0 = sprintf('%s.csv',fnam0);
    fname1 = sprintf('%s.csv',fnam1);
    % [num_data, text_data] = xlsread('Cell_1a_min0_2010.csv');
    % xind0=find((strcmp(text_data(:,1),'dXPos')));
    % yind0=find((strcmp(text_data(:,1),'dYPos')));
    % xpos0=num_data(xind0-1,1);
    % ypos0=num_data(yind0-1,1);
    % [num_data, text_data] = xlsread('Cell_1a_p800ish_2012.csv');
    % xind1=find((strcmp(text_data(:,1),'dXPos')));
    % yind1=find((strcmp(text_data(:,1),'dYPos')));
    % xpos1=num_data(xind1-1,1);
    % ypos1=num_data(yind1-1,1);
    
    TR0 = csvread(fname0,1,0);
    TR1 = csvread(fname1,1,0);
    else
        fnam1=fnam{n2+1};
        fname1 = sprintf('%s.csv',fnam1);
        TR0 = TRfin;
        TR1 = csvread(fname1,1,0);
    end
    TR0(:,5)=TR0(:,5)-0*add;
    TR1(:,5)=TR1(:,5)-1*add;
    
    x0 = TR0(:,3);
    y0 = TR0(:,4);
    x1 = TR1(:,3);
    y1 = TR1(:,4);
    minnie = round(min([x0;y0;x1;y1]));
    x0 = x0+64-minnie;
    y0 = y0+64-minnie;
    x1 = x1+64-minnie;
    y1 = y1+64-minnie;
    img1 = zeros(round(max(x1)/32),round(max(y1)/32));
    img0 = zeros(round(max(x0)/32),round(max(y0)/32));
    step1 = round(length(x1)/190000);
    step0 = round(length(x0)/190000);
    for n = 1:step1:length(x1)
        img1(round(x1(n)/32),round(y1(n)/32))=img1(round(x1(n)/32),round(y1(n)/32))+1;
    end
    for n = 1:step0:length(x0)
        img0(round(x0(n)/32),round(y0(n)/32))=img0(round(x0(n)/32),round(y0(n)/32))+1;
    end
    [a,b] = size(img0);
    [a2,b2] = size(img1);
    t1 = ones(a,b); t2=ones(a2,b2);
    % t12 = xcorr2(t1,t2);
    img0 = imgaussfilt(img0,2);
    img1 = imgaussfilt(img1,2);
    img0(find(img0>10))=2.5;
    img1(find(img1>10))=2.5;
    
    hey1 = xcorr2(img0,img1);
    %  hey = hey1./t12;
    %  hey2 = gradient(gradient(hey));
    %  hey2 = (hey2.*(hey.^2));
    [sz1,sz2] = size(img1);
    % [A1,B1] = find(hey==max(max(hey(:))));
    % [A2,B2] = find(hey2==min(min(hey2(:))));
    
    figure;
    % findpeaks(hey1(:),'Threshold',30,'Npeaks',10,'MinPeakHeight',0.2*max(hey1(:)))
 
    % [pks,locs] = findpeaks(hey1(:),'Threshold',30,'Npeaks',10,'MinPeakHeight',0.2*max(hey1(:)))
 
    [row,column] = find(hey1==max(max(hey1(:,:))));
    
    yshift = row-ceil(sz1);
    xshift = column-ceil(sz2);
    %yshift2 = (1000/32)*(xpos1-xpos0)
    %xshift2 = -(1000/32)*(ypos1-ypos0)
    TR1(:,3) = (TR1(:,3))+32*yshift;
    TR1(:,4) = (TR1(:,4))+32*xshift;
    
    %TR1(:,12)=TR1(:,12)+120;
    frames = 10000;
    TR1(:,2 ) = TR1(:,2)+fact*frames;
    
    figure; plot(TR0(1:round(length(x0)/60000):length(x0),3),TR0(1:(round(length(x0)/60000)):(length(x0)),4),'.','markersize',.01)
    hold on; plot(TR1(1:60000,3),TR1(1:60000,4),'.','markersize',.01)
    figure; imagesc(hey1)
    %TR1(:,12)=TR1(:,12)+120;
    
    
    TRfin = [TR0;TR1];
    TRfin(:,3) = (TRfin(:,3))-min(TRfin(:,3));
    TRfin(:,4) = (TRfin(:,4))-min(TRfin(:,4));
    TRfin(:,1) = 1:length(TRfin(:,1));
    if n2==1
     im0c=ones(a,b).*mean(mean(img0));
    else 
        im0c=double(summa2)/10000;
    end
 im1c=ones(a2,b2).*mean(mean(img1));
  [sz3,sz4] = size(img0);
 yshift3 = row-ceil(sz3);
xshift3 = column-ceil(sz4); 

 [summa,summa2]= imcombo(im0c,im1c,[yshift3,xshift3]);
%place='F:\Ben Data\20190208_COS7\Cell1001\Mito_Seg\From_nader_chambers'
figure; imagesc(double(summa2)); colormap gray 

end
imwrite((summa2(:,:)),[ savefile '.tif']);

 %%
 A6 = {'id','frame','x [nm]','y [nm]','z [nm]','sigma1 [nm]','sigma2 [nm]','intensity [photon]','offset [photon]','bgstd','chi2','uncertainty [nm]'};

 writecell(A6,savefile)
% 
 dlmwrite(savefile,TRfin,'delimiter',',','-append','precision','%.7f');
 %%
 %im0c=imread('buffalo_1a_2_fii.csv.tif');
%%
[summa3,summa4]= imcombo(img0,img1,[yshift3,xshift3]);
%place='F:\Ben Data\20190208_COS7\Cell1001\Mito_Seg\From_nader_chambers'
figure; imagesc(double(summa4)); colormap gray
%%
function [summa,summa2]= imcombo(A,B,RB)
[xa,ya]=size(A);
[xb,yb]=size(B);
if xa>xb
    B=padarray(B,[xa-xb,0],'pre');
elseif xb>xa
    A=padarray(A,[xb-xa,0],'pre');
end
if ya>yb
    B=padarray(B,[0,ya-yb],'pre');
elseif yb>ya
    A=padarray(A,[0,yb-ya],'pre');
end

if (RB(1)<0)&&(RB(2)<0)
    B=padarray(B,[abs(RB(1)),abs(RB(2))],'post');
    A=padarray(A,[abs(RB(1)),abs(RB(2))],'pre');
elseif (RB(1)<0)
    B=padarray(B,[abs(RB(1)),0],'post');
    A=padarray(A,[abs(RB(1)),0],'pre');
    B=padarray(B,[0,RB(2)],'pre');
    A=padarray(A,[0,RB(2)],'post');
elseif (RB(2)<0)
    B=padarray(B,[abs(RB(1)),0],'pre');
    A=padarray(A,[abs(RB(1)),0],'post');
    B=padarray(B,[0,abs(RB(2))],'post');
    A=padarray(A,[0,abs(RB(2))],'pre');
else 
    B=padarray(B,[RB(1),RB(2)],'pre');
    A=padarray(A,[RB(1),RB(2)],'post');
end
summa(:,:,1)=uint16(A*10000);
summa(:,:,3)=uint16(B*10000);
summa2=(summa(:,:,1)+summa(:,:,3));
%summa2=summa2*255/(max(max(summa2)));

end
