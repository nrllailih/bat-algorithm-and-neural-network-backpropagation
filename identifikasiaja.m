clc;
jb = input  ('masukkan jumlah kelelawar : ');
iter= input ('masukkan banyak iterasi: ');
A= input ('loudness awal: ');
r= input ('pulse rate awal: ');
g= input ('koefisien peningkatan pulse rate: ');
a= input ('koefisien peningkatan loudness: ');
alfa= input ('masukkan nilai learning rate: ');
maxit= input('masukkan max iterasi identifikasi: ');



if jb<1 
    disp ('jumlah kelelawar harus lebih dari satu')
 
end
fmin = input ('masukkan frekuensi minimal : ');
fmax = input ('masukkan frekuesi maksimal : ');
posisiawalkelelawar=zeros(jb,5);
kecepatanawalkelelawar=zeros(jb,5);
for i=1:jb
    for j=1:5
        x(i,j)=rand();
        kecepatanawalkelelawar(i,j) = rand();
        posisiawalkelelawar(i,j)= 1/(1+exp(-(x(i,j))));
    end
    
end
posisiawalkelelawar;
kecepatanawalkelelawar
    delta= posisiawalkelelawar(:,1);
    omega= posisiawalkelelawar(:,2);
    miu= posisiawalkelelawar(:,3);
    epsilon= posisiawalkelelawar(:,4);
    teta= posisiawalkelelawar(:,5);
    Data= load('datafull.txt');
    S= Data(:,1);
    I= Data(:,2);
    R= Data(:,3);
    M= Data(:,4);
    Sh= zeros(jb,12);
    Ih= zeros(jb,12);
    Rh= zeros(jb,12);
 
for a=1:iter
    disp('iterasi ke: '), disp (a);
    for i=1:jb
        disp('kelelawar'), disp (i);
        posisiawalkelelawar(i,:)
        kecepatanawalkelelawar(i,:)
    end
        
    %proses runge kuttes 
        for j=1:60-1
            if j==1
       
                Sh(i,j)= S(1,1);
                Ih(i,j)= I(1,1);
                Rh(i,j)= R(1,1);
            end
            k1(i,j)= 10000*delta(i,1)-omega(i,1)*Ih(i,j)*Sh(i,j)/(Sh(i,j)+Ih(i,j)+Rh(i,j))+miu(i,1)*Sh(i,j);
            l1(i,j)= omega(i,1)*Ih(i,j)*Sh(i,j)/(Sh(i,j)+Ih(i,j)+Rh(i,j));
            m1(i,j)= epsilon(i,1)*Ih(i,j)-miu(i,1)*Rh(i,j);
            k2(i,j)= 10000*delta(i,1)-omega(i,1)*(Ih(i,j)+l1(i,j)/2)*(Sh(i,j)+k1(i,j)/2)/(Sh(i,j)+k1(i,j)/2+Ih(i,j)+l1(i,j)/2+Rh(i,j)+m1(i,j)/2)+miu(i,1)+(Sh(i,j)+k1(i,j)/2);
            l2(i,j)= omega(i,1)*(Ih(i,j)+l1(i,j)/2)*(Sh(i,j)+k1(i,j)/2)/(Sh(i,j)+k1(i,j)/2+Ih(i,j)/2+Rh(i,j)+m1(i,j)/2)-(epsilon(i,1)+miu(i,1)+teta(i,1))*(Ih(i,j)+l1(i,j)/2);
            m2(i,j)= epsilon(i,1)*(Ih(i,j)+l1(i,j)/2)-miu(i,1)*(Rh(i,j)+m1(i,j)/2);
            k3(i,j)= 10000*delta(i,1)-omega(i,1)*(Ih(i,j)+l2(i,j)/2)*(Sh(i,j)+k2(i,j)/2)/(Sh(i,j)+k2(i,j)/2+Ih(i,j)+l2(i,j)/2+Rh(i,j)+m2(i,j)/2)+miu(i,1)*(Sh(i,j)+k2(i,j)/2);
            l3(i,j)= omega(i,1)*(Ih(i,j)+l2(i,j)/2)*(Sh(i,j)+k2(i,j)/2)/(Sh(i,j)+k2(i,j)/2+Ih(i,j)+l2(i,j)/2+Rh(i,j)+m2(i,j)/2)-(epsilon(i,1)+miu(i,1)+teta(i,1))*(Ih(i,j)+l2(i,j)/2);
            m3(i,j)= epsilon(i,1)*(Ih(i,j)+l2(i,j)/2)-miu(i,1)*(Rh(i,j)+m2(i,j)/2);
            k4(i,j)= 10000*delta(i,1)-omega(i,1)*(Ih(i,j)+l3(i,j))*(Sh(i,j)+k3(i,j))/(Sh(i,j)+k3(i,j)+Ih(i,j)+l3(i,j)+Rh(i,j)+m3(i,j))+miu(i,1)*(Sh(i,j)+k3(i,j));
            l4(i,j)= omega(i,1)*(Ih(i,j)+l3(i,j))*(Sh(i,j)+k3(i,j))/(Sh(i,j)+k3(i,j)+Ih(i,j)+l3(i,j)+Rh(i,j)+m3(i,j))-(epsilon(i,1)+miu(i,1)+teta(i,1))*(Ih(i,j)+l3(i,j));
            m4(i,j)= epsilon(i,1)*(Ih(i,j)+l3(i,j))-miu(i,1)*(Rh(i,j)+m3(i,j));
    
   
            Sh(i,j+1)= Sh(i,j)+(k1(i,j)+2*k2(i,j)+2*k3(i,j)+k4(i,j))/6;
            Ih(i,j+1)= Ih(i,j)+(l1(i,j)+2*l2(i,j)+2*l3(i,j)+l4(i,j))/6;
            Rh(i,j+1)= Rh(i,j)+(m1(i,j)+2*m2(i,j)+2*m3(i,j)+m4(i,j))/6;
        end
    end
 
%Menentukan MMRE dan fitness
jmh= zeros(jb,1);
jmh(1,1)=0;
fitness= zeros(jb,1);
for i=1:jb
   for j=1:60
        mmres(i,j)= abs(S(j,1)+Sh(i,j))/S(j,1);
        mmrei(i,j)= abs(I(j,1)+Ih(i,j))/I(j,1);
        mmrer(i,j)= abs(R(j,1)+Rh(i,j))/R(j,1);
        rata2(i,j)=(mmres(i,j)+mmrei(i,j)+mmrer(i,j))/3;
        jmh(i,1)=jmh(i,1)+rata2(i,j);
    end
    mmre(i,1)=jmh(i,1)/jb;
    fitness(i,1)=1/mmre(i,1);
end
fitness
%Mengurutan fitness
%fmax= max(fitness(i,1))
fitbest= 0;
 
for i=1:jb
    if  fitbest<fitness(i,1);
        fitbest=fitness(i,1);
          xbest=posisiawalkelelawar(i,:);
          batbest=i;
    end
end
batbest
xbest
betai=zeros(jb,1);  
frekuensi=zeros(jb,1);
for i=1:jb
   betai(i,1)= rand ();
   frekuensi(i,1)= fmin+(fmax-fmin)*betai(i,1);
end
betai
frekuensi
m=zeros(jb,5);
for i=1:jb
    m(i,:)= posisiawalkelelawar(batbest,:);
end
 
for i=1:jb
    for j=1:5
        kb(i,j)=(kecepatanawalkelelawar(i,j)+(posisiawalkelelawar(i,j)-m(i,j)));   
    end
end
posisibaru=zeros(jb,5);
for i=1:jb
    for j=1:5
        v(i,j)=frekuensi(i,1)*kb(i,j);
        kecepatanbaru(i,j)=1/(1+exp(-(v(i,j))));
        xx(i,j)= kecepatanbaru(i,j) + posisiawalkelelawar(i,j);
        posisibaru(i,j)=1/(1+exp(-(xx(i,j))));
   end
end
posisibaru
delta1= posisibaru(:,1);
omega1= posisibaru(:,2);
miu1= posisibaru(:,3);
epsilon1= posisibaru(:,4);
teta1= posisibaru(:,5);
 for i=1:jb
        disp('kelelawar'), disp (i);
        posisibaru(i,:)
        kecepatanbaru(i,:)
        
    %proses runge kuttes 
        for j=1:60-1
            if j==1
       
                Sh(i,j)= S(1,1);
                Ih(i,j)= I(1,1);
                Rh(i,j)= R(1,1);
            end
            k1(i,j)= 10000*delta1(i,1)-omega1(i,1)*Ih(i,j)*Sh(i,j)/(Sh(i,j)+Ih(i,j)+Rh(i,j))+miu1(i,1)*Sh(i,j);
            l1(i,j)= omega1(i,1)*Ih(i,j)*Sh(i,j)/(Sh(i,j)+Ih(i,j)+Rh(i,j));
            m1(i,j)= epsilon1(i,1)*Ih(i,j)-miu1(i,1)*Rh(i,j);
            k2(i,j)= 10000*delta1(i,1)-omega1(i,1)*(Ih(i,j)+l1(i,j)/2)*(Sh(i,j)+k1(i,j)/2)/(Sh(i,j)+k1(i,j)/2+Ih(i,j)+l1(i,j)/2+Rh(i,j)+m1(i,j)/2)+miu1(i,1)+(Sh(i,j)+k1(i,j)/2);
            l2(i,j)= omega1(i,1)*(Ih(i,j)+l1(i,j)/2)*(Sh(i,j)+k1(i,j)/2)/(Sh(i,j)+k1(i,j)/2+Ih(i,j)/2+Rh(i,j)+m1(i,j)/2)-(epsilon1(i,1)+miu1(i,1)+teta1(i,1))*(Ih(i,j)+l1(i,j)/2);
            m2(i,j)= epsilon1(i,1)*(Ih(i,j)+l1(i,j)/2)-miu1(i,1)*(Rh(i,j)+m1(i,j)/2);
            k3(i,j)= 10000*delta1(i,1)-omega1(i,1)*(Ih(i,j)+l2(i,j)/2)*(Sh(i,j)+k2(i,j)/2)/(Sh(i,j)+k2(i,j)/2+Ih(i,j)+l2(i,j)/2+Rh(i,j)+m2(i,j)/2)+miu1(i,1)*(Sh(i,j)+k2(i,j)/2);
            l3(i,j)= omega1(i,1)*(Ih(i,j)+l2(i,j)/2)*(Sh(i,j)+k2(i,j)/2)/(Sh(i,j)+k2(i,j)/2+Ih(i,j)+l2(i,j)/2+Rh(i,j)+m2(i,j)/2)-(epsilon1(i,1)+miu1(i,1)+teta1(i,1))*(Ih(i,j)+l2(i,j)/2);
            m3(i,j)= epsilon1(i,1)*(Ih(i,j)+l2(i,j)/2)-miu1(i,1)*(Rh(i,j)+m2(i,j)/2);
            k4(i,j)= 10000*delta1(i,1)-omega1(i,1)*(Ih(i,j)+l3(i,j))*(Sh(i,j)+k3(i,j))/(Sh(i,j)+k3(i,j)+Ih(i,j)+l3(i,j)+Rh(i,j)+m3(i,j))+miu1(i,1)*(Sh(i,j)+k3(i,j));
            l4(i,j)= omega1(i,1)*(Ih(i,j)+l3(i,j))*(Sh(i,j)+k3(i,j))/(Sh(i,j)+k3(i,j)+Ih(i,j)+l3(i,j)+Rh(i,j)+m3(i,j))-(epsilon1(i,1)+miu1(i,1)+teta1(i,1))*(Ih(i,j)+l3(i,j));
            m4(i,j)= epsilon1(i,1)*(Ih(i,j)+l3(i,j))-miu1(i,1)*(Rh(i,j)+m3(i,j));
    
   
            Sh(i,j+1)= Sh(i,j)+(k1(i,j)+2*k2(i,j)+2*k3(i,j)+k4(i,j))/6;
            Ih(i,j+1)= Ih(i,j)+(l1(i,j)+2*l2(i,j)+2*l3(i,j)+l4(i,j))/6;
            Rh(i,j+1)= Rh(i,j)+(m1(i,j)+2*m2(i,j)+2*m3(i,j)+m4(i,j))/6;
        end
 end
    %Menentukan MMRE dan fitness
jmh= zeros(jb,1);
jmh(1,1)=0;
fitnessbaru=zeros(jb,1);
for i=1:jb
   for j=1:60
        mmres(i,j)= abs(S(j,1)+Sh(i,j))/S(j,1);
        mmrei(i,j)= abs(I(j,1)+Ih(i,j))/I(j,1);
        mmrer(i,j)= abs(R(j,1)+Rh(i,j))/R(j,1);
        rata2(i,j)=(mmres(i,j)+mmrei(i,j)+mmrer(i,j))/3;
        jmh(i,1)=jmh(i,1)+rata2(i,j);
    end
    mmre(i,1)=jmh(i,1)/jb;
    fitnessbaru(i,1)=1/mmre(i,1);
end
mmre
fitnessbaru
pulserate=zeros(jb,1);
nilaiacak=zeros(jb,1);
loudnessawal=zeros(jb,1);
for i=1:jb
    nilaiacak(i,1)= rand ();
    pulserate(i,1)= r;
    loudnessawal(i,1)=A;
end
nilaiacak
pulserate
loudnessawal
jbest=0;
jnonbest=0;
j=0;
%%for i=1:jb
   % for j=1:5
    %if (pulserate(i,1)<nilaiacak(i,1))&&(fitness(i,1)<fitnessbaru(i,1))
        %localsearch
     %   temp(i,j)= posisibaru(i,j);
    %else
     %   temp(i,j)=posisiawalkelelawar(i,j);
    %end
    %end
%end
%temp
        
for i=1:jb 
    if (pulserate(i,1)<nilaiacak(i,1))
        
        jbest= jbest+1;
    else
        jnonbest=jnonbest+1;
    end
end
jbest
jnonbest
best=zeros(jbest,1);
nonbest=zeros(jnonbest,1);
p=0;
for i=1:jb 
    if (pulserate(i,1)<nilaiacak(i,1))
        j=j+1;
        best(j,1)=i;
    else
        p=p+1;
        nonbest(p,1)=i;
        
    end
end
best
nonbest
%bat tidak melakukan localsearch
xxbestnls= zeros(jnonbest,5);
fitnessnls= zeros(jnonbest,1);
fitnessave1= zeros (jnonbest,1);
xawalnls= zeros(jnonbest,5);
for i=1:jnonbest
    xxbestnls(i,:)=posisibaru(nonbest(i,1),:);
    fitnessnls(i,:)= fitnessbaru(nonbest(i,1),1);
    fitnessave1(i,:)= fitness(nonbest(i,1),1);
    xawalnls(i,:)=posisiawalkelelawar(nonbest(i,1),:);
end
xxbestnls
fitnessnls
%membandingkan fitnes bat lama dan bat baru (localsearch)
xxbest=zeros(jbest,5);
fitnessave=zeros(jbest,1);
for i=1:jbest
    if fitness(best(i,1),1)<fitnessbaru(best(i,1),1)
        xxbest(i,:)=posisibaru(best(i,1),:);
        finessave(i,1)= fitnessbaru(best(i,1),1);
    else
        xxbest(i,:)=posisiawalkelelawar(best(i,1),:);
        fitnessave(i,1)=fitness(best(i,1),1);
    end
end
xxbest
fitnessave
 
eps=zeros(jbest,5);
xls=zeros(jbest,5);
for i=1:jbest
    for j=1:5
    eps(i,j) = rand();
    xxx(i,j)= (xxbest(i,j)+eps(i,j))*A;
    xls(i,j)=1/(1+exp(-(xxx(i,j))));
    end
end
xls
%proses rungekutte localsearch
delta2= xls(:,1);
omega2= xls(:,2);
miu2= xls(:,3);
epsilon2= xls(:,4);
teta2= xls(:,5);
for i=1:jbest
        disp('kelelawar localsearch ke'), disp (i);
        xls(i,:)        
    %proses runge kuttes 
        for j=1:60-1
            if j==1
       
                Sh(i,j)= S(1,1);
                Ih(i,j)= I(1,1);
                Rh(i,j)= R(1,1);
            end
            k1(i,j)= 10000*delta2(i,1)-omega2(i,1)*Ih(i,j)*Sh(i,j)/(Sh(i,j)+Ih(i,j)+Rh(i,j))+miu2(i,1)*Sh(i,j);
            l1(i,j)= omega2(i,1)*Ih(i,j)*Sh(i,j)/(Sh(i,j)+Ih(i,j)+Rh(i,j));
            m1(i,j)= epsilon2(i,1)*Ih(i,j)-miu2(i,1)*Rh(i,j);
            k2(i,j)= 10000*delta2(i,1)-omega2(i,1)*(Ih(i,j)+l1(i,j)/2)*(Sh(i,j)+k1(i,j)/2)/(Sh(i,j)+k1(i,j)/2+Ih(i,j)+l1(i,j)/2+Rh(i,j)+m1(i,j)/2)+miu2(i,1)+(Sh(i,j)+k1(i,j)/2);
            l2(i,j)= omega2(i,1)*(Ih(i,j)+l1(i,j)/2)*(Sh(i,j)+k1(i,j)/2)/(Sh(i,j)+k1(i,j)/2+Ih(i,j)/2+Rh(i,j)+m1(i,j)/2)-(epsilon2(i,1)+miu2(i,1)+teta2(i,1))*(Ih(i,j)+l1(i,j)/2);
            m2(i,j)= epsilon2(i,1)*(Ih(i,j)+l1(i,j)/2)-miu2(i,1)*(Rh(i,j)+m1(i,j)/2);
            k3(i,j)= 10000*delta2(i,1)-omega2(i,1)*(Ih(i,j)+l2(i,j)/2)*(Sh(i,j)+k2(i,j)/2)/(Sh(i,j)+k2(i,j)/2+Ih(i,j)+l2(i,j)/2+Rh(i,j)+m2(i,j)/2)+miu2(i,1)*(Sh(i,j)+k2(i,j)/2);
            l3(i,j)= omega2(i,1)*(Ih(i,j)+l2(i,j)/2)*(Sh(i,j)+k2(i,j)/2)/(Sh(i,j)+k2(i,j)/2+Ih(i,j)+l2(i,j)/2+Rh(i,j)+m2(i,j)/2)-(epsilon2(i,1)+miu2(i,1)+teta2(i,1))*(Ih(i,j)+l2(i,j)/2);
            m3(i,j)= epsilon2(i,1)*(Ih(i,j)+l2(i,j)/2)-miu2(i,1)*(Rh(i,j)+m2(i,j)/2);
            k4(i,j)= 10000*delta2(i,1)-omega2(i,1)*(Ih(i,j)+l3(i,j))*(Sh(i,j)+k3(i,j))/(Sh(i,j)+k3(i,j)+Ih(i,j)+l3(i,j)+Rh(i,j)+m3(i,j))+miu2(i,1)*(Sh(i,j)+k3(i,j));
            l4(i,j)= omega2(i,1)*(Ih(i,j)+l3(i,j))*(Sh(i,j)+k3(i,j))/(Sh(i,j)+k3(i,j)+Ih(i,j)+l3(i,j)+Rh(i,j)+m3(i,j))-(epsilon2(i,1)+miu2(i,1)+teta2(i,1))*(Ih(i,j)+l3(i,j));
            m4(i,j)= epsilon2(i,1)*(Ih(i,j)+l3(i,j))-miu2(i,1)*(Rh(i,j)+m3(i,j));
    
   
            Sh(i,j+1)= Sh(i,j)+(k1(i,j)+2*k2(i,j)+2*k3(i,j)+k4(i,j))/6;
            Ih(i,j+1)= Ih(i,j)+(l1(i,j)+2*l2(i,j)+2*l3(i,j)+l4(i,j))/6;
            Rh(i,j+1)= Rh(i,j)+(m1(i,j)+2*m2(i,j)+2*m3(i,j)+m4(i,j))/6;
        end
end
 %Menentukan MMRE dan fitness localsearch
jmh= zeros(jbest,1);
jmh(1,1)=0;
fitnessbaruls=zeros(jbest,1);
for i=1:jbest
   for j=1:60
        mmres(i,j)= abs(S(j,1)+Sh(i,j))/S(j,1);
        mmrei(i,j)= abs(I(j,1)+Ih(i,j))/I(j,1);
        mmrer(i,j)= abs(R(j,1)+Rh(i,j))/R(j,1);
        rata2(i,j)=(mmres(i,j)+mmrei(i,j)+mmrer(i,j))/3;
        jmh(i,1)=jmh(i,1)+rata2(i,j);
    end
    mmre(i,1)=jmh(i,1)/jbest;
    fitnessbaruls(i,1)=1/mmre(i,1);
end
    fitnessbaruls
 
%melakukan perubahan loudness dan pulserate localsearch
Amat=zeros(jbest,1);
ramat=zeros(jbest,1);
Au=zeros(jbest,1);
ru=zeros(jbest,1); 
xakhirls= zeros(jbest,5);
ffixls=zeros(jbest,1);
 for i=1:jbest
  Amat(i,1)=A;
  ramat(i,1)=r;
  disp('kelelawar localsearch ke'), disp (best(i));
  xls(i,:)
  for j=1:5
             if (fitnessbaruls(i,1)<fitnessave(i,1))&& (pulserate(i,1)<loudnessawal(i,1))
             xakhirls(i,j)=xxbest(i,j);
             Au(i,1)=Amat(i,1);
             ru(i,1)=ramat(i,1);
             ffixls(i,1)=fitnessave(i,1);
         else
              xakhirls(i,j)=xls(i,j);
              Au(i,1)=a*Amat(i,1);
              ru(i,1)=ramat(i,1)*(1-exp(-g*a));
              ffixls(i,1)=fitnessbaruls(i,1);
           end
        
  end
 end
     xakhirls
     ffixls
     Au
     ru
 %melakukan perubahan loudness dan pulserate biasa
 Amat1=zeros(jnonbest,1);
 ramat1=zeros(jnonbest,1);
 Au1=zeros(jnonbest,1);
 ru1=zeros(jnonbest,1); 
 xakhirnls= zeros(jnonbest,5);
 ffix1=zeros(jnonbest,1);
 for i=1:jnonbest
  Amat1(i,1)=A;
  ramat1(i,1)=r;
  disp('kelelawar nonlocalsearch ke'), disp (nonbest(i));
  xxbestnls(i,:)
  for j=1:5
             if (fitnessnls(i,1)<fitnessave1(i,1))&& (pulserate(i,1)<loudnessawal(i,1))
             xakhirnls(i,j)=xawalnls(i,j);
             Au1(i,1)=Amat1(i,1);
             ru1(i,1)=ramat1(i,1);
             ffix1(i,1)=fitnessave1(i,1);
         else
              xakhinrls(i,j)=xxbestnls(i,j);
              Au1(i,1)=a*Amat1(i,1);
              ru1(i,1)=ramat1(i,1)*(1-exp(-g*a));
              ffix1(i,1)=fitnessnls(i,1);
           end
        
  end
 end
  xakhirnls
  ffix1
  Au1
  ru1
%menentukan parameter terbaik nls
fitbest1= 0;
for i=1:jnonbest
    if fitbest1<ffix1(i,1)
        fitbest1= ffix1(i,1);
        xbest1= xakhirnls(i,:);
    end
end
%menentukan parameter terbaik ls
fitbest2= 0;
for i=1:jbest
    if fitbest2<ffixls(i,1)
        fitbest= ffixls(i,1);
        xbest2= xakhirls(i,:);
    end
end
%menentukan parameter terbaik nls dan ls
fakhir=0;
xxakhirr= zeros(1,5);
if fitbest<fitbest1
    fakhir= fitbest1;
    xakhirr= xbest1;
else
    fakhir=fitbest;
    xakhirr=xbest2;
 
end
fakhir
xakhirr
 
BB= 10000*xakhirr(1,1);
omegafix=xakhirr(1,2);
miufix= xakhirr(1,3);
epsilonfix= xakhirr(1,4);
tetafix= xakhirr(1,5);

%normalisasi
smax= max(S);
imax= max(I);
rmax= max(R);
smin= min(S);
imin= min(I);
rmin= min(R);
mmax= max(M);
mmin= min(M);
xs= zeros(60,1);
xi= zeros(60,1);
xr= zeros(60,1);
xm= zeros(60,1);
for i=1:60
    xs(i,1)= (2*(S(i,1)-smin)/(smax-smin))-1;
    xi(i,1)= (2*(I(i,1)-imin)/(imax-imin))-1;
    xr(i,1)= (2*(R(i,1)-rmin)/(rmax-rmin))-1;
    xm(i,1)= (2*(M(i,1)-mmin)/(mmax-mmin))-1;
end
xs;
xi;
xr;
xm;
%data pelatihan
SP= [ xs xm xs];
IP= [ xi xm xi];
RP= [ xr xm xr];
%konversi parameter menjadi bobot dan bias
MIH= zeros(2,2);
mih= zeros(2,2);
bmih= [BB/10000 BB/10000];
mihfix= zeros(3,2);
MIHI= zeros(3,2);
mhoi= zeros(3,1);
MIHR= zeros(3,2);
mhor= zeros(3,1);
mho= zeros(3,1);
mhou= zeros(3,1);
mihfixu= zeros(3,2);
dltws= zeros(3,2);
dltvs= zeros(3,1);
dltwi= zeros(3,2);
dltvi= zeros(3,1);
dltwr= zeros(3,2);
dltvr= zeros(3,1);
%untuk S
for j=1:2
   for k=1:2
       MIH(j,k)=2*rand()-1; %merandom bilangan [-1,1]              
   end
end
disp ('Matrik bobot dan bias input ke hidden untuk S')
mihfix= [MIH; bmih]
disp ('Matrik bobot dan bias hidden ke output')
for l=1:3
    for m=1:1
        mho(l,m)= 2*rand()-1;
    end
end
mho
%untuk I
for j=1:3
    for k=1:2
        MIHI(j,k)=2*rand()-1;
    end
end
disp('Matriks bobot dan bias input ke hidden untuk I')
MIHI
for l=1:3
    for m=1:1
        mhoi(l,m)=2*rand()-1;
    end
end
disp('Matriks bobot dan bias hidden ke output untuk I')
mhoi

%untuk R
for j=1:3
    for k=1:2
        MIHR(j,k)=2*rand()-1;
    end
end
disp('Matriks bobot dan bias input ke hidden untuk R')
MIHR
for l=1:3
    for m=1:1
        mhor(l,m)=2*rand()-1;
    end
end
disp('Matriks bobot dan bias hidden ke output untuk R')
mhor

it=1;
%while ertr>err|it<maxit 
for p=1:maxit
    
%disp('Iterasi ke '), disp(p)
%it=it+1 

%disp('_____________________________________________________________');


for i=1:60
    %disp ('Untuk pola ke'), disp (i)
%feedforward dan backpropagation error S
    zin1s(i)=zeros;
    zin1s(i)= mihfix(3,1)+(mihfix(1,1)*SP(i,1)+mihfix(2,1)*SP(i,2));   
    zin2s(i)=zeros;
    zin2s(i)= mihfix(3,2)+(mihfix(1,2)*SP(i,1)+mihfix(2,2)*SP(i,2));
    z1s(i)=zeros;
    z1s(i)= (1-exp(-2*zin1s(i)))/(1+exp(-2*zin1s(i)));
    z2s(i)=zeros;
    z2s(i)= (1-exp(-2*zin2s(i)))/(1+exp(-2*zin2s(i)));
    yin1s(i)=zeros;
    yin1s(i)= mho(3,1)+((z1s(i)*mho(1,1)+z2s(i)*mho(2,1)));
    y1s(i)=zeros;
    y1s(i)= (1-exp(-2*yin1s(i)))/(1+exp(-2*yin1s(i)));
    ds(i)=zeros;
    ds(i)= (SP(i,3)-y1s(i))*(1-y1s(i))*(1+y1s(i));
    dltws = [alfa*ds(i)*z1s(i); alfa*ds(i)*z2s(i); alfa*ds(i)];
    sigmains1(i)=zeros;
    sigmains1(i)= ds(i)*dltws(1,1);
    sigmains2(i)=zeros;
    sigmains2(i)= ds(i)*dltws(2,1);
    sigmas1(i)=zeros;
    sigmas1(i)= sigmains1(i)*(1-z1s(i))*(1+z1s(i));
    sigmas2(i)=zeros;
    sigmas2(i)= sigmains2(i)*(1-z2s(i))*(1+z2s(i));
    dltvs= [ alfa*sigmas1(i)*SP(i,1) alfa*sigmas1(i)*SP(i,2); alfa*sigmas2(i)*SP(i,1) alfa*sigmas2(i)*SP(i,2); alfa*sigmas1(i) alfa*sigmas2(i)];
    mihfix=mihfix+dltvs;
    mho= mho+dltws;
    a(i)= (SP(i,3)-y1s(i)).^2;

%feedforward dan backpropagation error I    
    zin1i(i)=zeros;
    zin1i(i)= MIHI(3,1)+(MIHI(1,1)*IP(i,1)+MIHI(2,1)*IP(i,2));   
    zin2i(i)=zeros;
    zin2i(i)= MIHI(3,2)+(MIHI(1,2)*IP(i,1)+MIHI(2,2)*IP(i,2));
    z1i(i)=zeros;
    z1i(i)= (1-exp(-2*zin1i(i)))/(1+exp(-2*zin1i(i)));
    z2i(i)=zeros;
    z2i(i)= (1-exp(-2*zin2i(i)))/(1+exp(-2*zin2i(i)));
    yin1i(i)=zeros;
    yin1i(i)= mhoi(3,1)+((z1i(i)*mhoi(1,1)+z2i(i)*mhoi(2,1)));
    y1i(i)=zeros;
    y1i(i)= (1-exp(-2*yin1i(i)))/(1+exp(-2*yin1i(i)));
    di(i)=zeros;
    di(i)= (IP(i,3)-y1i(i))*(1-y1i(i))*(1+y1i(i));
    dltwi = [alfa*di(i)*z1i(i); alfa*di(i)*z2i(i); alfa*di(i)];
    sigmaini1(i)=zeros;
    sigmaini1(i)= di(i)*dltwi(1,1);
    sigmaini2(i)=zeros;
    sigmaini2(i)= di(i)*dltwi(2,1);
    sigmai1(i)=zeros;
    sigmai1(i)= sigmaini1(i)*(1-z1i(i))*(1+z1i(i));
    sigmai2(i)=zeros;
    sigmai2(i)= sigmaini2(i)*(1-z2i(i))*(1+z2i(i));
    dltvi= [ alfa*sigmai1(i)*IP(i,1) alfa*sigmai1(i)*IP(i,2); alfa*sigmai2(i)*IP(i,1) alfa*sigmai2(i)*IP(i,2); alfa*sigmai1(i) alfa*sigmai2(i)];
    MIHI=MIHI+dltvi;
    mhoi= mhoi+dltwi;
    
    b(i)= (IP(i,3)-y1i(i)).^2;
    
%feedforward dan backpropagation error R
    zin1r(i)=zeros;
    zin1r(i)= MIHR(3,1)+(MIHR(1,1)*RP(i,1)+MIHR(2,1)*RP(i,2));   
    zin2r(i)=zeros;
    zin2r(i)= MIHR(3,2)+(MIHR(1,2)*RP(i,1)+MIHR(2,2)*RP(i,2));
    z1r(i)=zeros;
    z1r(i)= (1-exp(-2*zin1r(i)))/(1+exp(-2*zin1r(i)));
    z2r(i)=zeros;
    z2r(i)= (1-exp(-2*zin2r(i)))/(1+exp(-2*zin2r(i)));
    yin1r(i)=zeros;
    yin1r(i)= mhor(3,1)+((z1r(i)*mhor(1,1)+z2r(i)*mhor(2,1)));
    y1r(i)=zeros;
    y1r(i)= (1-exp(-2*yin1r(i)))/(1+exp(-2*yin1r(i)));
    dr(i)=zeros;
    dr(i)= (RP(i,3)-y1r(i))*(1-y1r(i))*(1+y1r(i));
    dltwr = [alfa*dr(i)*z1r(i); alfa*dr(i)*z2r(i); alfa*dr(i)];
    sigmainr1(i)=zeros;
    sigmainr1(i)= dr(i)*dltwr(1,1);
    sigmainr2(i)=zeros;
    sigmainr2(i)= dr(i)*dltwr(2,1);
    sigmar1(i)=zeros;
    sigmar1(i)= sigmainr1(i)*(1-z1r(i))*(1+z1r(i));
    sigmar2(i)=zeros;
    sigmar2(i)= sigmainr2(i)*(1-z2r(i))*(1+z2r(i));
    dltvr= [ alfa*sigmai1(i)*RP(i,1) alfa*sigmar1(i)*RP(i,2); alfa*sigmar2(i)*RP(i,1) alfa*sigmar2(i)*RP(i,2); alfa*sigmar1(i) alfa*sigmar2(i)];
    MIHR=MIHR+dltvr;
    mhor= mhor+dltwr;
    c(i)= (RP(i,3)-y1r(i)).^2;
   
end
disp('MSE hasil identifikasi : ')
    ertr=(sum(a)+sum(b)+sum(c))/180
end
ertr

disp ('matriks bobot dan bias optimal dari input ke hidden untuk S')
   mihfix
   mho
   disp ('matriks bobot dan bias optimal dari input ke hidden untuk I')
   MIHI
   mhoi
   disp ('matriks bobot dan bias optimal dari input ke hidden untuk R')
   MIHR
   mhor
%validasi
for i=1:60
%untuk S
  zin1s(i)=zeros;
  zin1s(i)= mihfix(3,1)+(mihfix(1,1)*SP(i,1)+mihfix(2,1)*SP(i,2));   
  zin2s(i)=zeros;
  zin2s(i)= mihfix(3,2)+(mihfix(1,2)*SP(i,1)+mihfix(2,2)*SP(i,2));
  z1s(i)=zeros;
  z1s(i)= (1-exp(-2*zin1s(i)))/(1+exp(-2*zin1s(i)));
  z2s(i)=zeros;
  z2s(i)= (1-exp(-2*zin2s(i)))/(1+exp(-2*zin2s(i)));
  yin1s(i)=zeros;
  yin1s(i)= mho(3,1)+((z1s(i)*mho(1,1)+z2s(i)*mho(2,1)));
  y1s(i)=zeros;
  y1s(i)= (1-exp(-2*yin1s(i)))/(1+exp(-2*yin1s(i)));
  
%denormalisasi S
 % disp('hasil denormalisasi I')
  denos(i)= zeros;
  denos(i)= ((y1s(i)+1)*(smax-smin)/2)+smin;
  av(i)=abs(S(i)-denos(i))/S(i);
%untuk I
  zin1i(i)=zeros;
  zin1i(i)= MIHI(3,1)+(MIHI(1,1)*IP(i,1)+MIHI(2,1)*IP(i,2));   
  zin2i(i)=zeros;
  zin2i(i)= MIHI(3,2)+(MIHI(1,2)*IP(i,1)+MIHI(2,2)*IP(i,2));
  z1i(i)=zeros;
  z1i(i)= (1-exp(-2*zin1i(i)))/(1+exp(-2*zin1i(i)));
  z2i(i)=zeros;
  z2i(i)= (1-exp(-2*zin2i(i)))/(1+exp(-2*zin2i(i)));
  yin1i(i)=zeros;
  yin1i(i)= mhoi(3,1)+((z1i(i)*mhoi(1,1)+z2i(i)*mhoi(2,1)));
  y1i(i)=zeros;
  y1i(i)= (1-exp(-2*yin1i(i)))/(1+exp(-2*yin1i(i)));
  
%denormalisasi I
  %disp('hasil denormalisasi I')
  denoi(i)= zeros;
  denoi(i)= ((y1i(i)+1)*(imax-imin)/2)+imin;
  bv(i)=abs(I(i)-denoi(i))/I(i);
%untuk R
  zin1r(i)=zeros;
  zin1r(i)= MIHR(3,1)+(MIHR(1,1)*RP(i,1)+MIHR(2,1)*RP(i,2));   
  zin2r(i)=zeros;
  zin2r(i)= MIHR(3,2)+(MIHR(1,2)*RP(i,1)+MIHR(2,2)*RP(i,2));
  z1r(i)=zeros;
  z1r(i)= (1-exp(-2*zin1r(i)))/(1+exp(-2*zin1r(i)));
  z2r(i)=zeros;
  z2r(i)= (1-exp(-2*zin2r(i)))/(1+exp(-2*zin2r(i)));
  yin1r(i)=zeros;
  yin1r(i)= mhor(3,1)+((z1r(i)*mhor(1,1)+z2r(i)*mhor(2,1)));
  y1r(i)=zeros;
  y1r(i)= (1-exp(-2*yin1r(i)))/(1+exp(-2*yin1r(i)));
  
  %denormalisasi R
  %disp('hasil denormalisasi R')
  denor(i)= zeros;
  denor(i)= ((y1r(i)+1)*(rmax-rmin)/2)+rmin;
  cv(i)=abs(R(i)-denor(i))/R(i);
end
disp('error hasil validasi: ')
errval= (sum(av)+sum(bv)+sum(cv))/180
XS= [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60];
figure(1)
plot(XS,S,'-r',XS,denos,'-g')
ylabel('Susceptible population')
xlabel('month')
legend('data','result', 'Location','Best')
figure(2)
plot(XS,I,'-r',XS,denoi,'-g')
ylabel('Infected population')
xlabel('month')
figure(3)
plot(XS,R,'-r',XS,denor,'-g')
ylabel('Recovery population')
xlabel('month')
