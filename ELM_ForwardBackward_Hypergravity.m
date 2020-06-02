
close all
%lt is +y     and     rt is -y
clc;clear all


global ktm mset negmset mrandix velnow4beta velpast4beta mold mnew cold cnew ktmgt4 swcori disply disply2 centerborder pauseordisplay_interval generaltempstorageshelf betaloadedengleg disenglegload c_lt_begin c_lt_end c_rt_begin c_rt_end cend flagcom flagactiveleg alphatauml slopebeta gravity goal APmusclegoalbound MLmusclebound correcting vMLtrigthresh runtillthisktm runtillthisktmATedge randloopdirectionchooseswitchgoingUP randloopdirectionchooseswitchgoingDN rando randnu switch_x_or_z_tork swmlV lots_of_noise amp amplow Ataumid edgeloopcontrol midVcontrol APswgoal fdeltaktm APedgestiffens MLedgestiffens APvthresh4goal brake

starttime=datestr(now)

%CONSTANTS
gravity=9.8*1.8;
mub=75 ;wtub=mub*gravity;
ank2ank=.11;.27;
stra2a=num2str(fix(ank2ank*100));
y=[0 1 0];anklt=(ank2ank/2).*y;ankrt=(ank2ank/2).*(-y);



%PARAM SETTING / DIALS
slopetaurange=0:50:1500;slopebetarange=0:.5:15;


APmusclegoalbound=0.04;%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
MLmusclebound=0.02;%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
APgoalthresh=APmusclegoalbound;%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
APvthresh4goal=0.01/2;



PARAMRANGE1=slopetaurange;
PARAMRANGE2=slopebetarange;



%SWITCHES
swcori=0;
swm=1;
randinitials=1;
lots_of_noise=1;
keepnoise=1;
switch_x_or_z_tork=0;
swmlV=0;
edgeloopcontrol=0;
midVcontrol=0;
APswgoal=1;
APedgestiffens=1;
MLedgestiffens=0;

%SIMULATION SETTINGS
tstop=20;
centerborder =0.01/4;0.004/2;
iterruns=10;
if length(PARAMRANGE1)>1 iterend=iterruns;else iterend=1;end
fdeltaktm=4;




%DISPLAY SETTINGS
disply2=0;
plotrealtime=0;
if length(PARAMRANGE1)==1 plotrealtime=1;end
byby=5;
sync=1;syncdisp=500;forcedisplaylim=2000;
shortperiod=20;
tauxz_vs_swV_testing=0;
tickparam1=100;
tickparam2=1;
















%INITILIZERS
somecounter=0;shitfallen=0;
ksrow=0;
if length(PARAMRANGE1)>1 
    donematrix=zeros(length(PARAMRANGE1),length(PARAMRANGE2));
    P_notfall_matrix_iter=donematrix;
    t_fall_matrix=donematrix;
    var_t_fall_matrix=donematrix;   
end














% Initial position etc values  %TAKE VARIED VALUES  %NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
kkx=1;                                                                                                                                                                                                                                                                                                                 
kky=1/3;
hcomub=0.9;%ht wrt to center of ankles
if randinitials==0
    xlean0=-0.05*kkx;ylean0=-0.05*kky;

    z0=sqrt(hcomub^2-xlean0^2-ylean0^2);
    rcomub0=[xlean0 ylean0 z0];
end
    v0=[0.08*1/2 0.08*0 0]*2;%TAKE VARIED VALUES  %NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@







%VARIABLE VALUES
ttot=tstop;25;
dt=0.01;h=dt;
Nbintot=ttot/dt;



% PARAMETER SPACE SEARCH BEGIN
for  alphatauml=PARAMRANGE1%Ataumid=PARAMRANGE1  %alphatauml=slopetaurange
    ksrow=ksrow+1;kdecol=0;
    for  slopebeta=PARAMRANGE2%tdurktm=PARAMRANGE2   %slopebeta=slopebetarange
                    kdecol=kdecol+1;

     

    countiter=0;clear time2fall_end 
    falliter=0;notfalliter=0;
    [slopebeta alphatauml]

for iterrun=1:iterend
    iterrun               
                    
                    shitfallen=0;onceplot=0;
                    c=0;
                    somecounter=somecounter+1;
                    clear sr sv br
                    
                    
                    
                    
                    
                    if randinitials==1
                        xlean0=(rand-0.5)/10;
                        ylean0=(2*rand-1)*.02;
                        z0=sqrt(hcomub^2-xlean0^2-ylean0^2);
                        rcomub0=[xlean0 ylean0 z0];
                        v0(1)=-sign(xlean0)*abs(v0(1));
                    end

                    rcomo=rcomub0;
                    vo=v0;
                    ktm=1;
                    htt(ktm)=norm(rcomo);httlt(ktm)=norm(rcomo-anklt);httrt(ktm)=norm(rcomo-ankrt);
                    for i=1:3 sr(ktm,i)=rcomo(i);sv(ktm,i)=vo(i);br(ktm,1)=0;end

                    goal=1;
                    correcting=0;
                    ncorrkting_MidREGION_VCORREC=0;
                    runtillthisktmATedge=0;
                    runtillthisktm=0;
                    delktm=0;

                    %Below all is time bin dependent for COMub
                    t=(ktm-1)*dt;tm(ktm)=t;
                    ktm=2;


                    % Initial FLAGGING lt rt or center %lt is +y and rt is -y
                    mlproj=dot(rcomo,y);absmlproj=abs(mlproj); if absmlproj>centerborder  if mlproj>0 flag0 =1;else flag0 =-1;end;  else flag0=0;end     
                    flaglast=flag0;


                    

           































































            while ktm<=Nbintot+1 

            t=(ktm-1)*dt;
            ktmgt4=1;if ktm<4 ktmgt4=0;end

    
            if abs(mlproj)>centerborder flag_rtlt_lasttime=flag0;
            end



            if sync==1 pauseordisplay_interval=syncdisp;end; 
            if ktm==shortperiod pauseordisplay_interval=shortperiod;end
            if length(PARAMRANGE1)==1 pauseordisplay_interval=shortperiod;end

            if ~lots_of_noise
                if keepnoise
                   rando=rand(1,10);randnu=randn(1,10);
                else
                    rando=zeros(1,10);randnu=rando;
                end
            else
            rando=rand(1,16);randnu=randn(1,16);
            end
brake=0;
            % RUNGE KUTTA Multivariable integration (r and l for com ; v and k for velocity)
            disply=0;pauseordisplay1=20;if sync==1 pauseordisplay1=syncdisp;end;if ktm>forcedisplaylim pauseordisplay1=1; end
            if mod(ktm,pauseordisplay1)==0 disply=1;end ;disply=0;
            a=accln_for_no_coriolis_PRE_TRACE(rcomo,vo ,ank2ank,swm,c,mub,flag0,1);
            disply=0;
            k1 = h *a;
            l1 = h *vo;
            k2 = h *accln_for_no_coriolis_PRE_TRACE( rcomo+ l1./2, vo + k1./2,ank2ank,swm,c,mub,flag0,2);
            l2 = h *(vo + k1./2 );
            k3 = h *accln_for_no_coriolis_PRE_TRACE( rcomo+ l2./2, vo + k2./2,ank2ank,swm,c,mub,flag0,3);
            l3 = h *(vo + k2./2 );

          
            dispNplot= mod(ktm,pauseordisplay_interval)==0; 

            if dispNplot disply2=0;else disply2=0; end
            k4 = h *accln_for_no_coriolis_PRE_TRACE(rcomo+ l3, vo + k3,ank2ank,swm,c,mub,flag0,4);
            disply2=0;
            l4 = h *(vo + k3 );
            vn = vo + 1/6 *(k1 + 2*k2 + 2*k3 +  k4); 
            rn = rcomo  + 1/6 *(l1 + 2*l2 + 2*l3 +  l4); 

            





if APswgoal
            % Extra Forward backward kick with Extra ML torque
            if rn(1)>0 & rcomo(1)<0
                goal=+1;
            elseif rn(1)<0 & rcomo(1)>0
                goal=-1;
            else
                    if rn(1)>APgoalthresh
                    goal=-1;
                    elseif rn(1)<-APgoalthresh
                    goal=+1;
                    end
            end
end




if edgeloopcontrol
            % Controlling EDGE Loops
            
            if rn(1)>APmusclegoalbound & rcomo(1)<APmusclegoalbound
                randloopdirectionchooseswitchgoingUP=rand>=0.5;%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                edgedelktm=ifcorrectdothis_inPre0Cor(ktm);%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                runtillthisktmATedge=ktm+edgedelktm;
            elseif rn(1)<-APmusclegoalbound & rcomo(1)>-APmusclegoalbound
                randloopdirectionchooseswitchgoingDN=rand>0.8;%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                edgedelktm=ifcorrectdothis_inPre0Cor(ktm);%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                runtillthisktmATedge=ktm+edgedelktm;
            end
     
            if plotrealtime
                edgerang=[0 0 1];
                if randloopdirectionchooseswitchgoingUP==1&randloopdirectionchooseswitchgoingDN==1 edgerang=[1 1 0];end
                if randloopdirectionchooseswitchgoingUP==1&randloopdirectionchooseswitchgoingDN==0 edgerang=[1 0 0];end
                if randloopdirectionchooseswitchgoingUP==0&randloopdirectionchooseswitchgoingDN==1 edgerang=[0 1 0];end
                if randloopdirectionchooseswitchgoingUP==0&randloopdirectionchooseswitchgoingDN==0 edgerang=[0 0 0];end

                if randloopdirectionchooseswitchgoingUP==1 uparedgerang=[0 1 0];else  uparedgerang=[0 0 1];end

            end

end


if swmlV
            %Below for Velocity Threshold Corrections Control when calculating lateral muscle torques
            rap=rcomo(1);
            rml=rcomo(2);
            vap=vo(1);
            vml=vo(2);
            usetogetsmallvthresh=10;
            usetogetbigvthresh=1;
            if ~correcting


                  

                    vMLtrigthresh=0.02;
                    if abs(vml)>vMLtrigthresh
                        correcting=1;
                        runtillthisktm=ktm+ifcorrectdothis_inPre0Cor(ktm);
                    end

                    ncorrkting_MidREGION_VCORREC=ncorrkting_MidREGION_VCORREC+correcting
            end

             if ktm>runtillthisktm correcting=0;
             end

     end





















            % DATA STORAGE
                           sr(ktm,:)=rcomo(:);
                sv(ktm,:)=vo(:);
            
br(ktm)=brake;

















            % DISPLAY AND PLOTS
            if ktm>forcedisplaylim pauseordisplay_interval=1;end

            if dispNplot  %How often update figure
                pod=pauseordisplay_interval-1;
            
            
            mod9=mod(somecounter,byby^2);fix9=fix(somecounter/byby^2);
            if somecounter/byby^2-fix(somecounter/byby^2)~=0 figno=fix9+1;splotno=mod9;
            else figno=fix9;splotno=byby^2;end


            if plotrealtime
                
                
                
                  
                         
                    
                   
                    
                    
                    figure(figno);

                        if length(PARAMRANGE1)==1 extratit=['   ' 'ktm:' num2str(ktm)];
                        if tauxz_vs_swV_testing 
                           set(gcf,'position',[7 562 1271 421]);
                           subplot(1,2,2)
                            plot(ktm,switch_x_or_z_tork+0.1,'ro');hold on
                            plot(ktm,swmlV+0.2,'go');title('r=txz; g=vel')
                            ylim([-.05 1.5])
                           subplot(1,2,1)
                        else
                            set(gcf,'position',[6 494 530 490]);
                        end;
                    else subplot(byby,byby,splotno);set(gcf,'position',[4          15        1206         907]);extratit=[];end


                     hold on
                    lsr=size(sr,1);
                    colour=[0 1 0];
                    if ktm==pauseordisplay_interval plot(rcomub0(2),rcomub0(1),'k*');end
                    plot(sr(lsr,2),sr(lsr,1),'ro','MarkerSize',2)
                                                  if ~onceplot       
                         title(['K_{S}= ' num2str(alphatauml) '       K_{ED}= ' num2str(slopebeta) '   C:off' extratit]);
                         ylabel('ap (mt)');xlabel('ml (mt)');
                         axis([-0.1,0.1, -0.1, 0.1]*(1+.05));%1=x=ap  2=y=ml
                        set(gca,'XDir','rev');grid on
                        plot([-0.1 0.1],[0 0]);plot([0 0],[-0.1 0.1])
                         if length(PARAMRANGE1)==1
                            plot([0+centerborder 0+centerborder], [-0.1 0.1], 'y')
                            plot([0-centerborder 0-centerborder], [-0.1 0.1], 'y')
                            hv(1,APmusclegoalbound);
                            hv(1,-APmusclegoalbound);
                            hv(2,MLmusclebound);
                            hv(2,-MLmusclebound);
                         end
                                        onceplot=1;
                     end
                     
                  
%                        
            end

           
           
            end






















    
            outofbound=0.1;
            if abs(sr(ktm-1,1))>outofbound | abs(sr(ktm-1,2))>outofbound 
                shitfallen=1;tfall=ktm*dt;


             ktm=4000;t=tstop; 

           
            end %abs(sr(ktm-1,1))>0.1 | abs(sr(ktm-1,2))>0.1

                                                   




















            % FLAGGING revision lt rt or center %lt is +y and rt is -y
             
             mlproj=dot(rn,y);absmlproj=abs(mlproj);centerborderbound =centerborder; if absmlproj>centerborderbound  if mlproj>0 flag1 =1;else flag1 =-1;end;  else flag1=0;end     
            
            if exist('flag_rtlt_lasttime')==0 flag_rtlt_lasttime=0;end
            if flag1-flag_rtlt_lasttime==2 %%Rt--> Lt so Lt begin 
                             rcmlt=(rn-anklt);
                 rcaplt= rcmlt/norm(rcmlt);
                 vn=vn-dot(vn,rcaplt).*rcaplt;

            elseif flag1-flag_rtlt_lasttime==-2 %%Lt--> Rt so Rt begin 
                          rcmrt=(rn-ankrt);
                 rcaprt= rcmrt/norm(rcmrt);
                 vn=vn-dot(vn,rcaprt).*rcaprt;
            end







            


           




            rcomo=rn;vo=vn;flaglast=flag0;flag0=flag1; 
            ktm=ktm+1;
            end % END OF WHILE































            % Things to do after decided Fall or Not Fall 

           
            if shitfallen==1 
                countiter=countiter+1;
                falliter=falliter+1;
                
                time2fall_end(countiter)=tfall; 
            
            else    
                countiter=countiter+1;
                notfalliter=notfalliter+1;
                time2fall_end(countiter)=tstop;
              

            end
                   


            if plotrealtime
                figure(figno)
                  if shitfallen==1        
                      plot(sr(end,2),sr(end,1),'c*') 
                  else 
                      plot(sr(end,2),sr(end,1),'g*')
                  end
            plot(sr(:,2),sr(:,1),'r-')
            end
            
            
            
            
end %END ITERRUN























    if length(PARAMRANGE1)>1
           
         donematrix(ksrow,kdecol)=0.1;
         fD=1200;hD=figure(fD);
         surf(PARAMRANGE2,PARAMRANGE1,donematrix)%changed x<-->y in surf
         set(gca,'YTick',[0:tickparam1:PARAMRANGE1(end)])
         set(gca,'XTick',[0:tickparam2:PARAMRANGE2(end)])  
         zlim([0 1])
         
         zlabel('DONE')
         if ispc
             set(gcf,'position',[23 512 560 420])
         else
         set(gcf,'position',[20 555 560 420])
         end
        
        time2fall_end
        t_fall_matrix(ksrow,kdecol)=mean(time2fall_end);
        var_t_fall_matrix(ksrow,kdecol)=std(time2fall_end);
        P_notfall_matrix_iter(ksrow,kdecol)=notfalliter/(falliter+notfalliter);
        fP=700;hP=figure(fP);
        surf(PARAMRANGE2,PARAMRANGE1,P_notfall_matrix_iter)%changed x<-->y in surf
         set(gca,'YTick',[0:tickparam1:PARAMRANGE1(end)])
         set(gca,'XTick',[0:tickparam2:PARAMRANGE2(end)])  
         zlim([0 1])
         zlabel('P')
         view([0 90])
       
         if ispc set(gcf,'position',[616 516 560 420]);end
        
         fsz=15;set(gca,'fontsize', fsz,'fontweight','bold')
         set(gcf,'Name',[ 'Fall or not   AP Sway   cond=' num2str( fix(gravity/(9.8*1.8))+1) 'g'  ' APswgoal=' num2str(APswgoal) ' iter=' num2str(iterruns)])
          xlabel('K_{ED}')
          ylabel('K_{S}^{FB}')

          view([0 90])
          caxis([0 1])
          colormap jet   
         
        
        t_fall_matrix1=t_fall_matrix/max(max(t_fall_matrix));
        fT=500;hT=figure(fT);
        surf(PARAMRANGE2,PARAMRANGE1,t_fall_matrix)%changed x<-->y in surf
        set(gca,'YTick',[0:tickparam1:PARAMRANGE1(end)])
        set(gca,'XTick',[0:tickparam2:PARAMRANGE2(end)])  
        zlim([0 20])
        zlabel('T')
        if ispc 
            set(gcf,'position',[13 12 560 420])
        else
            set(gcf,'position',[1296 554 560 420]);
        end
      fsz=15;set(gca,'fontsize', fsz,'fontweight','bold')
         set(gcf,'Name',[ 'T to fall   AP Sway   cond=' num2str( fix(gravity/(9.8*1.8))+1) 'g' ' APswgoal=' num2str(APswgoal) ' iter=' num2str(iterruns)])
          xlabel('K_{ED}')
          ylabel('K_{S}^{FB}')

          view([0 90])
          caxis([0 20])
          colormap jet   
        
         disp('#################################################')
    end




    end %END OF KED SEARCH
end%END OF KS SEARCH

endtime=strrep(datestr(now),':','-')


saveas(gcf,[num2str( fix(gravity/(9.8*1.2))+1) 'g ' get(gcf,'name') '  a2a=' stra2a 'cm' '  ' endtime],'fig')





