%%%%%%%load data
Sample=importdata('SampleNameFile.txt');
for i=1:size(Sample,1)
    fileID = fopen(['./RNA/',Sample{i,1},'.txt']);
    C = textscan(fileID,'%s %f32');
    fclose(fileID);
    Exp(:,i)=double(C{1,2});
    Symbol=C{1,1};
    fileID = fopen(['./ATAC/',Sample{i,1},'.bed']);
    C = textscan(fileID,'%s %f32');
    fclose(fileID);
    Opn(:,i)=double(C{1,2});
    PeakName=C{1,1};
end
Exp=log2(1+quantilenorm(Exp));
Opn=log2(1+quantilenorm(Opn));
%%%%%%%%%%%RE_TG
fileID = fopen('peak_gene.txt');
C = textscan(fileID,'%s %s');
fclose(fileID);
[d1 f1]=ismember(C{1,2},Symbol);
A1=zeros(length(d1),size(Exp,2));
A1(d1==1,:)=Exp(f1(d1==1),:);
[d2 f2]=ismember(C{1,1},PeakName);
Opn1=Opn(f2,:);
Z1=zscore(A1');
Z2=zscore(Opn1');
r=sum(Z1.*Z2)/(size(Z1,1)-1);
r=r';
dlmwrite('corr',r,'\t')
dis=dlmread('dis');
RE_TG=[f1(d1) f2(d1) r(d1) dis(d1)];
Element_name=PeakName;
%save('Exp_Opn.mat','Exp','Opn','Symbol','Element_name','Sample','RE_TG')
%%%%TF-TG
TFName=importdata('TFName.txt');
[d,f]=ismember(TFName,Symbol);
TFExp=Exp(f(d==1),:);
TFName=TFName(d==1);
R2=corr(TFExp',Exp');
List=Symbol;
Exp_median=median(Exp')';
TFExp_median=median(TFExp')';
save('../Prior/TFTG_corr.mat','TFName','List','R2','Exp_median','TFExp_median')
