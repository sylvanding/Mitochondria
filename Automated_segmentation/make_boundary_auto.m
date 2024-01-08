filename='mito_uncut_conf_17';
b = readNPY([filename,'.npy']);
orig=imread([filename,'.tif']);
b2=permute(b,[2,3,1]);
b3=sum(b2,3);
b4=zeros(size(b3,1),size(b3,2),3);
b44=b4;
boundary=zeros(size(b3,1),size(b3,2));
boundary2=boundary;
b5=b4;
se=strel('disk',2);
px=.025;
minsz=5/px;
for n2=1:size(b2,3)
    cmap=rand(3,1);
    b2(:,:,n2) = bwareafilt(b2(:,:,n2),[minsz 500000000]);
    b2(:,:,n2)=imfill(b2(:,:,n2),'holes');
    object=bwconncomp(b2(:,:,n2));
    samples=object.PixelIdxList;
    clear maxtest;
    if length(object.PixelIdxList)==0
        continue
    end
    for n0=1:length(object.PixelIdxList)
        maxtest(n0)=length(samples{n0});
    end
    indo=find(maxtest==max(maxtest));
    sub=samples{indo};
    [row,col]=ind2sub(size(b3),sub);

    %[row,col]=find(b2(:,:,n2));
    
    for n3=1:3
        for n4=1:length(row)
            b4(row(n4),col(n4),n3)=b4(row(n4),col(n4),n3)+cmap(n3);
        
        end
    end
        for n4=1:length(row)
            norm=sum(b4(row(n4),col(n4),:));
            b4(row(n4),col(n4),:)=b4(row(n4),col(n4),:)/norm;
            b5(row(n4),col(n4),n2)=1;
        end

    B1=imdilate(b5(:,:,n2),se);
    B2=B1-b5(:,:,n2);
    boundary=boundary+B2;
    object5=bwconncomp(b5(:,:,n2));
    obsz5(n2)=length(object5.PixelIdxList{1});
end
    [obsz5srt,ind]=sort(obsz5);
    b6=zeros(size(b3));
    b5=logical(b5);
    for n3=1:size(b5,3)
        cmap=rand(3,1);
        b6=b6+double(b5(:,:,ind(n3)));
        if max(max(b6))>1
            b6=b6-double(b5(:,:,ind(n3)));
            b5t=b5(:,:,ind(n3))-b6;
            b5t(b5t<0)=0;
            b5(:,:,ind(n3))=b5t;
            b5(:,:,ind(n3)) = bwareafilt(b5(:,:,ind(n3)),[minsz 500000000]);
            b6=b6+double(b5(:,:,ind(n3)));
        end
        B7=imdilate(b5(:,:,ind(n3)),se);
        B8=B7-b5(:,:,ind(n3));
        boundary2=boundary2+B8;
        object5=bwconncomp(b5(:,:,ind(n3)));
        samples=object5.PixelIdxList;
        if length(samples)>0
            sub=samples{1};
            [row,col]=ind2sub(size(b3),sub);
        for n33=1:3
            for n4=1:length(row)
                b44(row(n4),col(n4),n33)=b44(row(n4),col(n4),n33)+cmap(n33);

            end
        end
        end
    end
    
    se=strel('disk',1);
    boundary3=imdilate(boundary2,se);
    boundary_inv=1-boundary3;

    final=orig.*uint8(boundary_inv);
    %%
    for k=1:size(b5,3)
        B{k} = bwboundaries(b5(:,:,k));

    end
    
figure; imagesc((orig),[0,900]);
colormap gray
hold on
blank=zeros(size(B1));
for k = 1:length(B)
   if ~isempty(B{1,k})
       boundary = B{1,k}{1,1};
       plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
       co=B{1,k}{1,1};
       for k2=1:size(co,1)
          blank(co(k2,1),co(k2,2))=1;
       end
       area(k)=length(B{1,k}{1,1});
   end
end
    %%
    imwrite(uint8(boundary2),[filename,'_boundary.tif']);
    final2=sum(final,3);
    norm=max(max(final2))./255;
    imwrite(uint8((final2./norm)),[filename,'_segmented.tif'])
    imwrite(b44,[filename,'_multi_col.tif'])
    %%
