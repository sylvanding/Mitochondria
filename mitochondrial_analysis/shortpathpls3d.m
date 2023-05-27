function [skel_fin,crackWidthImage,Length_fin]=shortpathpls3d(cell_img)
%cell_img=SPLMload(['test_skel','.tif'],'tiff');
cell_img=logical(cell_img);
[skel_test,link,node2]=bwskel2(cell_img,'MinBranchLength', 1);


blank1=zeros(size(cell_img,1)+2,size(cell_img,2)+2,size(cell_img,3)+2);
blank2=blank1;
% for n=1:length(link)
%    blank(link(n).point)=1;
% end
% for n=1:length(node2)
%    blank2(node2(n).idx)=1;
% end
% indices=find((blank2-blank)==1);
if ~isempty(link)
    s=[];
    t=[];
    weight=[];
    for n=1:length(link)
        s(n)=link(n).n1;
        t(n)=link(n).n2;
        weight(n)=length(link(n).point);
    end
    
    %endpoints =  bwmorph(skel_test,'endpoints');
    %nodes =  bwmorph(skel_test,'branchpoints');
    % [nx,ny]=find(nodes);
    % nodes=[nx,ny];
    % [ex,ey]=find(endpoints);
    % endpoints=[ex,ey]
    %figure; imagesc(uint8(skel_test)+uint8(cell_img)+uint8(endpoints)+uint8(nodes)); colormap gray
    mx_node=max([s,t]);
    s=[s,t];
    t=[t,s(1:(length(s)/2))];
    weight=[weight,weight];
    G = digraph(s,t,weight);
    ind3=0
    for ind1=1:mx_node
        for ind2=1:mx_node
            ind3=ind3+1;
            [no,pathlength(ind3)]=shortestpath(G,ind1,ind2);
            pathlength(isinf(pathlength))=0;
            if pathlength(ind3)==max(pathlength)
                noperm=no;
                pathperm=pathlength(ind3);
            end
        end
    end
    noperm;
    pathperm;
    blank=zeros(size(cell_img,1)+2,size(cell_img,2)+2,size(cell_img,3)+2);
    index=[1:length(link),1:length(link)];
    for n=2:length(noperm)
        indices=noperm(n-1:n);
        npath=find((indices(1)==s).*(indices(2)==t));
        npath=npath(1);
        
        blank(link(index(npath)).point)=1;
    end
    skel_fin=blank(2:(end-1),2:(end-1),2:(end-1));
else
    skel_fin=skel_test;
end
[x,y,z]=size(skel_fin);
indices=find(skel_fin);
surf=cell_img - imerode(cell_img, true(3));
disttrans=bwdist(surf);
[x2,y2,z2]=ind2sub([x,y,z],indices);
if length(indices)<=2
    Length_fin=1;
    indices2=find(surf);
    [x3,y3,z3]=ind2sub([x,y,z],indices2);
    x2=mean(x2); y2=mean(y2); z2=mean(z2);
    radii=mean(sqrt((x3-x2).^2+(y3-y2).^2+(z3-z2).^2));
    crackWidthImage = 2 * radii;
    
else
    skel_length(1)=sqrt(abs(x2(2)-x2(1))+abs(y2(2)-y2(1))+abs(z2(2)-z2(1)));
    for n=2:(length(indices)-1)
        
        tl=sqrt((x2(n+1)-x2(n-1))^2+(y2(n+1)-y2(n-1))^2+(z2(n+1)-z2(n-1))^2);
        slopexy(n-1)=(z2(n+1)-z2(n-1))/tl;
        slopeyz(n-1)=(x2(n+1)-x2(n-1))/tl;
        slopexz(n-1)=(y2(n+1)-y2(n-1))/tl;
        skel_length(n)=sqrt(abs(x2(n+1)-x2(n))+abs(y2(n+1)-y2(n))+abs(z2(n+1)-z2(n)));
        plane=crop_circ3d(size(cell_img,1),size(cell_img,2),size(cell_img,3),x2(n),y2(n),z2(n),slopexy(n-1),slopexz(n-1),slopeyz(n-1));
        ring=plane.*surf;
        indices2=find(ring);
        
        [x3,y3,z3]=ind2sub([x,y,z],indices2);
        radii(n-1)=mean(sqrt((x3-x2(n)).^2+(y3-y2(n)).^2+(z3-z2(n)).^2));
    end
    % for n=1:z
    %    if isempty(find(skel_fin(:,:,n)))
    %    else
    %        boundaries=bwboundaries(cell_img(:,:,n));
    %        boundaries=boundaries{1,1};
    %       [compx,compy]=find(skel_fin(:,:,n));
    %
    %        for n3=1:length(compx)
    %            distance=0;
    %             for n2=1:length(boundaries)
    %                 distance=distance+sqrt((boundaries(n2,1)-compx(n3))^2+(boundaries(n2,2)-compy(n3))^2);
    %             end
    %             distance=2*(distance/length(boundaries));
    %             dist_save(n)=distance;
    %        end
    %
    %    end
    % end
    % crackWidthImage=dist_save;
    % boundaries=(bwboundaries(cell_img));
    % boundaries=boundaries{1,1};
    % blank=zeros(size(cell_img,1),size(cell_img,2));
    % for n=1:length(boundaries)
    %     blank(boundaries(n,1),boundaries(n,2))=1;
    % end
    
    
    crackWidthImage = 2 * radii;
    Length_fin=sum(skel_length);
    %figure; imagesc(uint8(cell_img)+uint8(skel_fin)+uint8(nodes))
end
end
