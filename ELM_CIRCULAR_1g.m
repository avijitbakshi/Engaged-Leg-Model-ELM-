
close all
  
%lt is +y     and     rt is -y
clc;clear
datestr(now)






global ktm mset negmset mrandix velnow4beta velpast4beta mold mnew cold cnew ktmgt4 swcori disply centerborder pauseordisplay_interval generaltempstorageshelf betaloadedengleg disenglegload c_lt_begin c_lt_end c_rt_begin c_rt_end cend flagcom flagactiveleg alphatauml slopebeta gravity alphatauRL vtannedswitch

somecounter=0;count=0;shitfallen=0;
gravity=9.8;
ank2ank=.11;.27;stra2a=num2str(fix(ank2ank*100));

vtannedswitch=1;iterend=10;randinitials=1


slopetau2range=[0:2:20];slopebetarange=[0:10];slopetau1range=[0:100:1000];
slopetau2range=[0:2:20];slopebetarange=[0:15];slopetau1range=[0:100:1500];


byby=5;
ks1row=0;
lk1=length(slopetau1range);
lk2=length(slopetau2range);
lked=length(slopebetarange);


t_fall_matrix=zeros(lk2,lked,lk1);
donematrix=zeros(lk2,lked);





for alphatauml=slopetau1range;
    ks1row=ks1row+1;ks2row=0;
for alphatauRL=slopetau2range
    ks2row=ks2row+1;kdecol=0;
    for slopebeta=slopebetarange
        clear rn rcomo
         disp([alphatauml alphatauRL slopebeta]) 
        kdecol=kdecol+1;
 somecounter=somecounter+1;
        for iter=1:iterend

clear sr sv


%define parameter values reqd
kfig=0;tstop=20;%1.5;
firstplot=0;
c=0;
   %%%Beta
    basem=10;10;
    basec=0.37;mold=basem;cold=basec;
    msett=abs(basem+4*randn(1,100));
    msetucutoff=msett(find(msett<14));
    mset=msetucutoff(find(msetucutoff>1));lm=length(mset); negmset=-mset(find(mset<10));ln=length(negmset);
    mrandix=ceil(lm*rand);
    extraslopeenhancer=1;mset=mset+extraslopeenhancer;
    
    c_std=0.05/20;c_rt_end_mean=0.5;c_lt_end_mean=0.4;c_rt_end=rand*c_std + (c_rt_end_mean-c_std/2);c_lt_end=rand*c_std + (c_lt_end_mean-c_std/2);
    


justaftermidlinecrossing=0;



swcori=0;swm=1;%SWITCHES


ktm=1;
kkx=1/4;                                                                                                                                                                                                                                                                                                                 
kky=1;0.057;1/3;
% Initial position etc values 
% X IS AP; Y IS ML
hcomub=0.9;%ht wrt to center of ankles
xlean0=-0.05*kkx;ylean0=-0.05*kky;





if randinitials==1
                      
       %AVIJIT NOW CIRCULAR rand initials 
       rlower=0.05;rupper=0.07;
                        
                         rini0=rand*(rupper-rlower)+rlower;
                         thetaini0=rand*360;
                         xlean0=rini0*cosd(thetaini0);

                        ylean0=rini0*sind(thetaini0);
                        z0=sqrt(hcomub^2-xlean0^2-ylean0^2);
                        rcomub0=[xlean0 ylean0 z0];
                        rplane0=[xlean0 ylean0];
                        vr=0.08*1/2;
                        sn=sign(rand-0.5);
                        theta = 90;
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
v0= (rplane0./norm(rplane0)*R)*vr*sn;
v0=[v0 0];
                        
                    end


if      ylean0>0 c_lt_end=rand*c_std+ c_lt_end_mean;cend=c_lt_end;
elseif  ylean0<0 c_rt_end=rand*c_std+ c_rt_end_mean;cend=c_rt_end;
end


z0=sqrt(hcomub^2-xlean0^2-ylean0^2);
y=[0 1 0];anklt=(ank2ank/2).*y;ankrt=(ank2ank/2).*(-y);
rcomub0=[xlean0 ylean0 z0];rcomo=rcomub0;
htt(ktm)=norm(rcomo);httlt(ktm)=norm(rcomo-anklt);httrt(ktm)=norm(rcomo-ankrt);
theta=atan(abs(ank2ank/2-abs(ylean0))/hcomub);

vo=v0;for i=1:3 sr(ktm,i)=rcomo(i);sv(ktm,i)=vo(i);end
mub=75 ;wtub=mub*gravity;

% Initial FLAGGING lt rt or center %lt is +y and rt is -y

centerborder =0;0.01/4;0.004/2;
mlproj=dot(rcomo,y);absmlproj=abs(mlproj); if absmlproj>centerborder  if mlproj>0 flag0 =1;else flag0 =-1;end;  else flag0=0;end     
flaglast=flag0;

%Below all is time bin dependent for COMub
ttot=tstop;25;
dt=0.01;h=dt;
Nbintot=ttot/dt;t=(ktm-1)*dt;tm(ktm)=t;
ktm=2;


























while ktm<=Nbintot+1 
t=(ktm-1)*dt;
ktmgt4=1;if ktm<4 ktmgt4=0;end

sync=1;syncdisp=20;forcedisplaylim=2000; 
if abs(mlproj)>centerborder flag_rtlt_lasttime=flag0;end
% RUNGE KUTTA Multivariable integration (r and l for com ; v and k for velocity)
disply=0;pauseordisplay1=20;if sync==1 pauseordisplay1=syncdisp;end;if ktm>forcedisplaylim pauseordisplay1=1; end
if mod(ktm,pauseordisplay1)==0 disply=1;end ;disply=0;
a=accln__CIRCULAR_selectKS(rcomo,vo ,ank2ank,swm,c,mub,flag0,1);
disply=0;
k1 = h *a;
l1 = h *vo;
k2 = h *accln__CIRCULAR_selectKS( rcomo+ l1./2, vo + k1./2,ank2ank,swm,c,mub,flag0,2);
l2 = h *(vo + k1./2 );
k3 = h *accln__CIRCULAR_selectKS( rcomo+ l2./2, vo + k2./2,ank2ank,swm,c,mub,flag0,3);
l3 = h *(vo + k2./2 );
k4 = h *accln__CIRCULAR_selectKS( rcomo+ l3, vo + k3,ank2ank,swm,c,mub,flag0,4);
l4 = h *(vo + k3 );
vn = vo + 1/6 *(k1 + 2*k2 + 2*k3 +  k4); 
rn = rcomo  + 1/6 *(l1 + 2*l2 + 2*l3 +  l4) ;
rncap=rn/norm(rn);
vn=vn-dot(vn,rncap).*rncap;
z=[0 0 1];
vn=vn-dot(vn,z).*z;


% DATA STORAGE

     sr(ktm,:)=rcomo(:);
    sv(ktm,:)=vo(:);

pauseordisplay_interval=100;50;4;100;
if sync==1 pauseordisplay_interval=syncdisp;end;
if ktm>forcedisplaylim pauseordisplay_interval=1;end

if mod(ktm,pauseordisplay_interval)==0 %How often update figure
    pod=pauseordisplay_interval-1;
                                               
   
mod9=mod(somecounter,byby^2);fix9=fix(somecounter/byby^2);if somecounter/byby^2-fix(somecounter/byby^2)~=0 figno=fix9+1;splotno=mod9;else figno=fix9;splotno=byby^2;end

plotrealtime=0;
if plotrealtime
figure(figno);
       
set(gcf,'position',[19         236        1168         727]);

       
 subplot(byby,byby,splotno);hold on
lsr=length(sr);
colour=[0 1 0];

if ktm==pauseordisplay_interval plot(rcomub0(2),rcomub0(1),'k*');end
plot(sr(lsr,2),sr(lsr,1),'r*')
        
 title(['K_{S2}= ' num2str(alphatauRL) '    K_{ED}= ' num2str(slopebeta)     '   K_{S1}= ' num2str(alphatauml)   ]);ylabel('ap (mt)');xlabel('ml (mt)');
 axis([-0.1,0.1, -0.1, 0.1]);%1=x=ap  2=y=ml
set(gca,'XDir','rev');grid on
plot([-0.1 0.1],[0 0]);plot([0 0],[-0.1 0.1])
plot([0+centerborder 0+centerborder], [-0.1 0.1], 'r')
plot([0-centerborder 0-centerborder], [-0.1 0.1], 'r')

end

end







%%%%%%%%%%%%%%%PLOTTER AND ENDER-(plots when done)%%%%%%%%%%%%%%%%%%%%%
outofbound=0.1;
if abs(sr(ktm-1,1))>outofbound | abs(sr(ktm-1,2))>outofbound 
    shitfallen=1;tfall=ktm*dt;
    
    
    ktm=4000;t=tstop; 
   

end 


% FLAGGING revision lt rt or center %lt is +y and rt is -y
if justaftermidlinecrossing==1 %pause;
    justaftermidlinecrossing=0;end  
 mlproj=dot(rn,y);absmlproj=abs(mlproj);centerborderbound =centerborder; if absmlproj>centerborderbound  if mlproj>0 flag1 =1;else flag1 =-1;end;  else flag1=0;end     


if flag1-flag_rtlt_lasttime==2 %%Rt--> Lt so Lt begin 
     c_lt_begin = 1-cend;cold=c_lt_begin; % Now causality ok , cend expected ~ c_rt_end
     c_lt_end=rand*c_std+ c_lt_end_mean;cend=c_lt_end;
elseif flag1-flag_rtlt_lasttime==-2 %%Lt--> Rt so Rt begin 

     c_rt_begin = 1-cend;cold=c_rt_begin;%cend expected ~ c_lt_end
     c_rt_end=rand*c_std+ c_rt_end_mean;cend=c_rt_end;
end




rcomo=rn;vo=vn;flaglast=flag0;flag0=flag1; 

ktm=ktm+1;

end % end of while
if shitfallen==1 count=count+1 ;
    display('Fallen')

    tfallen(count)=tfall;
else    count=count+1   ;  
    

    tfallen(count)=tstop;
    disp('full time')
end
shitfallen=0;

        end %iter

t_fall_matrix(ks2row,kdecol,ks1row)=mean(tfallen);
count=0;
disp('---------------')




donematrix(ks2row,kdecol)=0.1*ks1row;
         fD=1200;hD=figure(fD);
         surf(slopebetarange,slopetau2range,donematrix)%changed x<-->y in surf
         set(gca,'YTick',[0:slopetau2range(end)])
         set(gca,'XTick',[0:slopebetarange(end)])  
         zlim([0 1])
         
         zlabel('DONE')
         if ispc
             set(gcf,'position',[23 512 560 420])
         else
         set(gcf,'position',[20 555 560 420])
         end




fosz=16;
t_fall_matrix20=t_fall_matrix/20;
         figure(2000+ks1row);
        surf(slopebetarange,slopetau2range,squeeze(t_fall_matrix20(:,:,ks1row)))%changed x<-->y in surf
         set(gca,'YTick',slopetau2range)
         set(gca,'XTick', slopebetarange)
          zlim([0 1])
          ylabel('KS2')
          xlabel('KED')
          view([0 90])
           set(gcf,'name',['KS1= ' num2str(slopetau1range(ks1row))])
           set(gca,'fontweight','bold','fontsize',fosz)
caxis([0 1])
set(gcf,'position',[ 1330          88         560         420])











    end 
end


    figure(2000+ks1row);
close(2000+ks1row)




end



fosz=16;
t_fall_matrix20=t_fall_matrix/20;
for iterk1=1:lk1 




         figure(6000+iterk1);
        surf(slopebetarange,slopetau2range,squeeze(t_fall_matrix20(:,:,iterk1)))%changed x<-->y in surf
 
         set(gca,'YTick',slopetau2range)
         set(gca,'XTick', slopebetarange)
        
          zlim([0 1])
%           ylabel('KS2')
%           xlabel('KED')
          view([0 90])
           set(gcf,'name',['KS1= ' num2str(slopetau1range(iterk1))])
            set(gca,'fontweight','bold','fontsize',fosz)
            caxis([0 1])
end
          




for iterk2=1:lk2 
         figure(7000+iterk2);
        surf(slopetau1range,slopebetarange,squeeze(t_fall_matrix20(iterk2,:,:)))%changed x<-->y in surf
  
         set(gca,'YTick',slopebetarange)
         set(gca,'XTick',slopetau1range) 
          zlim([0 1])
%           xlabel('KS1')
%           ylabel('KED')
          view([-270 -90])
          set(gcf,'name',['KS2= ' num2str(slopetau2range(iterk2))])
          set(gca,'fontweight','bold','fontsize',fosz)
end
datestr(now)


endtime=strrep(datestr(now),':','-')


for ifg=1:lk1
    figure(6000+ifg);BLACK_and_WHITE_KS_KED_and_Error;
    saveas(gcf,[num2str( fix(gravity/(9.8*1.2))+1) 'g ' get(gcf,'name') '  a2a=' stra2a 'cm' '  ' endtime],'fig')
end
for ifg=1:lk2
    figure(7000+ifg);BLACK_and_WHITE_KS_KED_and_Error;view([-270 -90])%coz BW code has  view([0 90]) in it so change it back 
    saveas(gcf,[num2str( fix(gravity/(9.8*1.2))+1) 'g ' get(gcf,'name') '  a2a=' stra2a 'cm' '  ' endtime],'fig')
end
beep on; for k=1:20 beep;pause(.1); end
error('end')
        

