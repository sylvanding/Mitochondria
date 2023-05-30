function [summa,summa2]=im_shift(cell_img,cell_img2)

%%
[xa,ya,ml]=size(cell_img);
[xb,yb,ml]=size(cell_img2);
if xa>xb
    cell_img2=padarray(cell_img2,[xa-xb,0],'pre');
  
elseif xb>xa
    cell_img=padarray(cell_img,[xb-xa,0,0],'pre');

end
if ya>yb
    cell_img2=padarray(cell_img2,[0,ya-yb],'pre');
 
elseif yb>ya
    cell_img=padarray(cell_img,[0,yb-ya,0],'pre');

end
%%
cell_imgs=sum(cell_img,3);
[a,b]=size(cell_imgs);
[a2,b2]=size(cell_img2);
fact=1;
c1=imresize(cell_imgs,[a/fact,b/fact]);
c2=imresize(cell_img2,[a2/fact,b2/fact]);
%  c1=([zeros(100,300),ones(100,100);zeros(200,200),zeros(200,200)])
%  c2=([zeros(100,70),zeros(100,70);ones(100,100),zeros(100,40)])
t1=ones(a/fact,b/fact); t2=ones(a2/fact,b2/fact);
a4=round(a/2); b4=round(b/2);
a5=round(a2/2); b5=round(b2/2);
t12= xcorr2(t1((a4-400):(a4+400),(b4-400):(b4+400)),t2((a4-400):(a4+400),(b4-400):(b4+400)));
hey=xcorr2(uint8(c1((a4-400):(a4+400),(b4-400):(b4+400))),uint8(c2((a4-400):(a4+400),(b4-400):(b4+400))));
hey=hey./t12;
hey2=gradient(gradient(hey));
sz1=801;
sz2=801;
[A,B]=find(hey2==min(min(hey2(sz1-300:sz1+300,sz2-300:sz2+300))));

yshift=A-ceil(sz1);
xshift=B-ceil(sz2); 

[summa,summa2]= imcombo(cell_img,cell_img2,[yshift*fact,xshift*fact]);
%place='F:\Ben Data\20190208_COS7\Cell1001\Mito_Seg\From_nader_chambers'
red(:,:,1)=sum(summa,3);
red(:,:,2)=sum(summa2,3);
red(:,:,3)=zeros(size(summa2,(1)),size(summa2,(2)));
end

function [summa,summa2]= imcombo(A,B,RB)
[xa,ya,ml]=size(A);
[xb,yb,ml]=size(B);
if xa>xb
    B=padarray(B,[xa-xb,0],'pre');
  
elseif xb>xa
    A=padarray(A,[xb-xa,0,0],'pre');

end
if ya>yb
    B=padarray(B,[0,ya-yb],'pre');
 
elseif yb>ya
    A=padarray(A,[0,yb-ya,0],'pre');

end

if (RB(1)<0)&&(RB(2)<0)
    B=padarray(B,[abs(RB(1)),abs(RB(2))],'post');
    A=padarray(A,[abs(RB(1)),abs(RB(2)),0],'pre');
    
elseif (RB(1)<0)
    B=padarray(B,[abs(RB(1)),0],'post');
    A=padarray(A,[abs(RB(1)),0,0],'pre');
    B=padarray(B,[0,RB(2)],'pre');
    A=padarray(A,[0,RB(2),0],'post');
elseif (RB(2)<0)
    B=padarray(B,[abs(RB(1)),0],'pre');
    A=padarray(A,[abs(RB(1)),0,0],'post');
    B=padarray(B,[0,abs(RB(2))],'post');
    A=padarray(A,[0,abs(RB(2)),0],'pre');
else 
    B=padarray(B,[RB(1),RB(2)],'pre');
    A=padarray(A,[RB(1),RB(2),0],'post');
end
summa=A;
summa2=B;
%summa2=summa2*255/(max(max(summa2)));

end
