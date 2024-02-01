function [verim]=specokuma(filename,scan_numarasi)
%written by MMT in Istanbul, ITU on July 13, 2016
%reads spec data

arast1=sprintf('#S %d',scan_numarasi);uz1=length(arast1);
arast2=sprintf('#N');uz2=length(arast2);
arast3=sprintf('#L');uz3=length(arast3);

fin=fopen(filename,'r');
oku1=fgetl(fin);
while(strncmpi(oku1,arast1,uz1)==0) %strncmpi return 1 if true, 0 if not
 oku1=fgetl(fin);
 if(feof(fin))
    error('\n bir hata var (1)');
    break;
 end;
end;

%%%%now arast1 bulundu %%%%%%
bos1=strfind(oku1,' '); % burda scan komutunu buluyoruz. 
scantipi1=oku1(bos1(3)+1:bos1(4)-1); %herzaman 3 ile 4 arasinda komut var.
scankomutu1=oku1(bos1(3)+1:end);% tum scan komutu

oku2=fgetl(fin);
while(strncmpi(oku2,arast2,uz2)==0)
 oku2=fgetl(fin);
  if(feof(fin))
    error('\n bir hata var (2)');
    break;
 end;
end;
bos2=strfind(oku2,' '); %colum sayisini okundugu satirdayiz
columsayisi=str2num(oku2(bos2:end));

%%%now Label satirini okunmasi %%%%%%
oku3=fgetl(fin);
while(strncmpi(oku3,arast3,uz3)==0)
 oku3=fgetl(fin);
  if(feof(fin))
    error('\n bir hata var (3)');
    break;
 end;
end;
bos3=strfind(oku3,' '); %#L satirinda ki bosluklar
basarasi=find(diff(bos3)>1);%Label satirindaki iki bosluk olan araliklar 
basliksayisi=length(basarasi); %bu easesn #N ile ayni olmali...
%%%% #L okundu ve simdi data okuma sirasi %%%%%%
veri=fscanf(fin,'%g',[columsayisi inf]);
verim=veri';

fclose(fin);
return;
%fprintf('\n komut:%s ve colum sayisi:%d',scankomutu1,columsayisi);
%fprintf('\n label:%s',oku3);
%pause;

