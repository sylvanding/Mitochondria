for k=1:1
simnames={'simulated_normal_1'; 'simulated_normal_2';'simulated_normal_3'; 'simulated_normal_4';'simulated_normal_5'; 'simulated_normal_6';'simulated_normal_7'; 'simulated_normal_8';'simulated_normal_9'; 'simulated_normal_10';'simulated_normal_11'; 'simulated_normal_12';'simulated_normal_13'; 'simulated_normal_14';'simulated_normal_15'; 'simulated_normal_16';'simulated_normal_17'; 'simulated_normal_18';'simulated_normal_19'; 'simulated_normal_20'; 'simulated_normal_21';; 'simulated_normal_22'};
gtnames={'gt_normal_fin1'; 'gt_normal_fin2';'gt_normal_fin3'; 'gt_normal_fin4';'gt_normal_fin5'; 'gt_normal_fin6';'gt_normal_fin7'; 'gt_normal_fin8';'gt_normal_fin9'; 'gt_normal_fin10';'gt_normal_fin11'; 'gt_normal_fin12';'gt_normal_fin13'; 'gt_normal_fin14';'gt_normal_fin15'; 'gt_normal_fin16';'gt_normal_fin17'; 'gt_normal_fin18';'gt_normal_fin19'; 'gt_normal_fin20'; 'gt_normal_fin21'; 'gt_normal_fin22'};

px=160; %nm
%dist=10; %nm
NA=1.4; %Numerica aperture
%nblinks=25;
%ndot=9;
%nit=35;
dens=1; %Labelling density
FOV=12; %Field of view in microns
nframes=500; %number of frames
anti_length=20; %Length of antibody in nm; double it if you are using primary and secondary
mito_diam=0.5; %Circular diameter of mitochondria (probably should not change)
ep_dens=15.49; %1 alpha/beta mitochondria epitope per 14.49 nm slong length of microtubule
k2=0.0135; %Chance of moving from on state to triplet state ms^-1
k3=0.333e-04; %Chance of moving from triplet state to on state ms^-1
k4= 0.0011;%Chance of moving from on state to bleached state ms^-1
frame_rate=20; %ms
zeroth_int=500; %Desired photon count of zeroth order
first_int=1500; %Desired photon count of first order
Ibgp = 15; %Bacground noise of zeroth order (photons*100)
Ibg = 3*Ibgp; %Background noise of first order (photons*100)
RN=100; %readout noise (photons*100)
mito_length=1.0; %um
heterogeneity=1; %Heterogeneity of mitochondria sizes; scales from 0 to 1
numChains=25; %number of mitochondria
dissociation=0;
%dist=15;
    %[a,b,tempx,tempy,zpos,spec_hetr]=sim_grid(px,dist);
    [mito,mito_edge,mito_viewable,tempx,tempy,zpos,spec_hetr]=sim_mito_3D(dens,px,numChains,anti_length,mito_diam,ep_dens,FOV,mito_length,heterogeneity);
    gee=(mito_viewable-mito_edge)>0;

%%
% for n=51:(size(mito_viewable,3))
%     if n==50
%         ffff=3
%     end
%     inters=(mito_viewable(:,:,n)>1);
%     %line=bwskel(inters);
%     mito_fin(:,:,n)=(mito_viewable(:,:,n)>0);
%     mito_fin=double(mito_fin);
%     guff(:,:,n)=edge(mito_fin(:,:,n),'canny');
%     objects=bwconncomp(inters);
%     for n4=1:length(objects.PixelIdxList)
%         inters2=zeros(size(mito_fin(:,:,n),1),size(mito_fin(:,:,n),2));
%         
%         [x1,y1]=ind2sub([size(mito_fin(:,:,n),1),size(mito_fin(:,:,n),2)],objects.PixelIdxList{1,n4})
%         for n5=1:length(x1)
%            inters2(x1(n5),y1(n5))=1; 
%         end
%         
%       
%         x=[]; y=[];
%         if length(x1)>1
%             if length(unique(x1))>length(unique(y1))
%                 x=unique(x1);
%                 for n2=min(x1):max(x1)
%                     y=[y;mean(find(inters2(n2,:)))];
% 
%                 end
%             else
%                 y=unique(y1);
%                 for n2=min(y1):max(y1)
%                     x=[x;mean(find(inters2(:,n2)))];
% 
%                 end
% 
%             end
% 
%         end
% 
%         if length(x)>0
%             for n3=1:length(x)
%                 mito_fin(round(x(n3)),round(y(n3)),n)=2;
%             end
%         end
%     end
% end
%%
% filename2='conf_psf2';
% cell_img=SPLMload([filename2,'.tif'],'tiff');
% [x,y,z]=size(cell_img);
% %cell_img=cell_img(90:170,90:170,150:350);
% %cell_img2=interp3(double(cell_img),-1);
% mito_edge2=imdilate(mito_edge,ones(2,2,2));
% mito_edge2=interp3(mito_edge2,-1);
% for n2=1:size(mito_edge2,3)
% mito_edge2(:,:,n2)=double(bwskel(logical(mito_edge2(:,:,n2)))).*mito_edge2(:,:,n2);
% end
% 
% in_mit_fact=1;
% edge_mit=1;
% for n=1:size(mito_edge2,3)
% inner=imfill(mito_edge2(:,:,n));
% inner=inner.*(rand(size(mito_edge2,1),size(mito_edge2,2))>0.875);
% mito_edge2(:,:,n)=mito_edge2(:,:,n).*(rand(size(mito_edge2,1),size(mito_edge2,2))>0.0);
% inner(inner>0)=in_mit_fact;
% mito_edget=mito_edge2(:,:,n)*edge_mit;
% backgr=rand(size(mito_edge2,1),size(mito_edge2,2))>0.995;
% fin_gt(:,:,n)=mito_edget+inner+backgr;
% end
% conf_img=convn(fin_gt,cell_img,'same');
% % Nsp = random('Poisson',conf_img);
% % Nr = random('norm', 0,RN,n11,n11,Nf);
% % conv3(mito_edge,)
%%




in_mit_fact=1;
edge_mit=1;
for n=1:size(mito_edge,3)
inner=imfill(mito_edge(:,:,n));
inner=inner.*(rand(size(mito_edge,1),size(mito_edge,2))>0.975);
mito_edge3(:,:,n)=mito_edge(:,:,n).*(rand(size(mito_edge,1),size(mito_edge,2))>0.55);
inner(inner>0)=in_mit_fact;
mito_edget=mito_edge3(:,:,n)*edge_mit;
backgr=rand(size(mito_edge3,1),size(mito_edge3,2))>0.998;
fin_gt2(:,:,n)=mito_edget+inner+backgr;
end
%%
vals=unique(fin_gt2);
gt=[];
GT_list=[];
chance=[];
for n=1:length(vals)
    chance(n)=length(find(fin_gt2==vals(n)))*double(vals(n));
    GT_list=[GT_list;repmat(find(fin_gt2==vals(n)),[vals(n),1])];
end
    
    [x,y,z] = ind2sub([size(fin_gt2,1),size(fin_gt2,2),size(fin_gt2,3)],GT_list);

x=16*(x+rand-.5); y=16*(y+rand-.5); z=16*(z+rand-.5);
x=x+dissociation*randn(length(x),1);
y=y+dissociation*randn(length(x),1);
z=z+dissociation*randn(length(x),1);
z0=mean(z);

chance2=sum(chance);
nblinks=300000;
order=randperm(length(y));
mu=3;
blinks=round(exprnd(mu,nblinks,1));
indices=[];
for n=1:(nblinks/mu)
    indices=[indices;repmat(n,blinks(n),1)];
end
xn=x(order);
yn=y(order);
zn=z(order);
% indices=round(pool*(rand(nlabels,1)));
% indices=nonzeros(indices);
TR0gt=(1:length(indices))';
TR0gt(:,2)=ones(length(indices),1);
TR0=TR0gt;
TR0gt(:,3:5)=[xn(indices),yn(indices),zn(indices)];
lat_unc=(15.96+10.37*randn(length(TR0gt(:,3)),1));%+abs((TR0gt(:,5)-z0))/60;
ax_unc=2.5641*lat_unc;
lat_unc=lat_unc;
TR0(:,3)=TR0gt(:,3)+lat_unc.*randn(length(TR0gt(:,3)),1);TR0(:,4)=TR0gt(:,4)+lat_unc.*randn(length(TR0gt(:,3)),1);TR0(:,5)=TR0gt(:,5)+ax_unc.*randn(length(TR0gt(:,3)),1);
% TR0(:,3:5)=[x0(indices),y0(indices),z0(indices)];

% 
filename = sprintf('%s', gtnames{k});
filename2=[filename,'.tif']; % Specify the output file name
gee=uint16(gee);
[~,~,d]=size(gee);
for idx = 1:d
    %[A,map] = rgb2ind(temp2(:,:,:,idx),256);
    if idx == 1
        imwrite(gee(:,:,idx),filename2);
    else
        imwrite(gee(:,:,idx),filename2,'WriteMode','append');
    end
end

A6={'id','frame','x [nm]','y [nm]','z [nm]'}%,'sigmay [nm]','intensity [photon]','offset [photon]','uncertainty [nm]', 'centroid [nm]','bgstd', 'z [nm]', 'z uncertainty','ydif','centroid2','fwhm0','fwhm1'};
%savefile='simulated_short_1.csv';
savefile = sprintf('%s', simnames{k});
writecell(A6,[savefile,'.csv'])
%writecell(A6,'eu_2.csv')

dlmwrite([savefile,'.csv'],TR0,'delimiter',',','-append');
clear;
end
%%
% filename='conf_succ2';
% filename2=[filename,'.tif']; % Specify the output file name
% factor=max(max(max(conf_img)))/65535;
% conf_imgf=uint16(conf_img/factor);
% [~,~,d]=size(conf_imgf);
% for idx = 1:d
%     %[A,map] = rgb2ind(temp2(:,:,:,idx),256);
%     if idx == 1
%         imwrite(conf_imgf(:,:,idx),filename2);
%     else
%         imwrite(conf_imgf(:,:,idx),filename2,'WriteMode','append');
%     end
% end
% 
% filename='gt_fin';
% filename2=[filename,'.tif']; % Specify the output file name
% fin_gt2=uint16(fin_gt2);
% [~,~,d]=size(fin_gt2);
% for idx = 1:d
%     %[A,map] = rgb2ind(temp2(:,:,:,idx),256);
%     if idx == 1
%         imwrite(fin_gt2(:,:,idx),filename2);
%     else
%         imwrite(fin_gt2(:,:,idx),filename2,'WriteMode','append');
%     end
% end