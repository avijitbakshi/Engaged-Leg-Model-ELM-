


close all

%lt is +y     and     rt is -y
clc;clear all
global ktm mset negmset mrandix velnow4beta velpast4beta mold mnew cold cnew ktmgt4 swcori disply disply2 centerborder pauseordisplay_interval generaltempstorageshelf betaloadedengleg disenglegload c_lt_begin c_lt_end c_rt_begin c_rt_end cend flagcom flagactiveleg alphatauml slopebeta gravity goal APmusclegoalbound MLmusclebound correcting vMLtrigthresh runtillthisktm runtillthisktmATedge randloopdirectionchooseswitchgoingUP randloopdirectionchooseswitchgoingDN rando randnu switch_x_or_z_tork swmlV lots_of_noise amp amplow Ataumid edgeloopcontrol midVcontrol APswgoal fdeltaktm APedgestiffens MLedgestiffens APvthresh4goal brake Vbrakethresh






%PARAM SETTING / DIALS
slopetaurange=0:50:1000;slopebetarange=0:10;

APmusclegoalbound=0.04;%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
MLmusclebound=0.01;%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
APgoalthresh=APmusclegoalbound;%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Vbrakethresh=0.1/4;


PARAMRANGE1=slopetaurange;

%SWITCHES
swcori=0;swm=1;
randinitials=1;
MLedgestiffens=1;




%DISPLAY SETTINGS
disply2=0;
plotrealtime=0;
if length(slopetaurange)==1 plotrealtime=1;end
byby=5;
sync=1;syncdisp=500;forcedisplaylim=2000;
shortperiod=20;
tauxz_vs_swV_testing=0;
markersizenum=8;
linewdth=2.5;



%SIMULATION SETTINGS
tstop=20;
centerborder =0.01/4;0.004/2;
lots_of_noise=0;





%CONSTANTS
gravity=9.8*1;%1.8;
mub=75 ;wtub=mub*gravity;
y=[0 1 0];ank2ank=.27;anklt=(ank2ank/2).*y;ankrt=(ank2ank/2).*(-y);








%INITILIZERS
somecounter=0;shitfallen=0;
ksrow=0;
if length(slopetaurange)>1 
    donematrix=zeros(length(slopetaurange),length(slopebetarange));
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
for alphatauml=slopetaurange
    ksrow=ksrow+1;kdecol=0;
    for slopebeta=slopebetarange
                    kdecol=kdecol+1;
                    
                    
     if length(slopetaurange)>1 iterend=10;else iterend=1;end

    countiter=0;clear time2fall_end 
    falliter=0;notfalliter=0;
    [slopebeta alphatauml]
for iterrun=1:iterend
    iterrun               
                    
                    shitfallen=0;onceplot=0;
                    c=0;

                    somecounter=somecounter+1;
                    clear sr sv
                    
                    
                    
                    
                    
                    if randinitials==1
                        xlean0=(rand-0.5)/10; 
                        ylean0=(2*rand-1)*.02;
                        z0=sqrt(hcomub^2-xlean0^2-ylean0^2);
                        rcomub0=[xlean0 ylean0 z0];
                        v0(1)=-sign(xlean0)*abs(v0(1));
                        if  alphatauml==450 slopebeta==4
                            xlean0=(2*rand-1)*.015;
                            z0=sqrt(hcomub^2-xlean0^2-ylean0^2);
                            rcomub0=[xlean0 ylean0 z0];
                        end
                    end

                    rcomo=rcomub0;
                    vo=v0;
                    ktm=1;
                    htt(ktm)=norm(rcomo);httlt(ktm)=norm(rcomo-anklt);httrt(ktm)=norm(rcomo-ankrt);
                    for i=1:3 sr(ktm,i)=rcomo(i);sv(ktm,i)=vo(i);end

                    goal=1;
                    correcting=0;
                    ncorrkting=0;
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



            if sync==1 pauseordisplay_interval=syncdisp;end; %26 May 2014 moved up , commented below
            if ktm==shortperiod pauseordisplay_interval=shortperiod;end
            if length(slopetaurange)==1 pauseordisplay_interval=shortperiod;end

            if ~lots_of_noise
            rando=rand(1,10);randnu=randn(1,10);
            else
            rando=rand(1,16);randnu=randn(1,16);
            end
brake=0;
            % RUNGE KUTTA Multivariable integration (r and l for com ; v and k for velocity)
            disply=0;pauseordisplay1=20;if sync==1 pauseordisplay1=syncdisp;end;if ktm>forcedisplaylim pauseordisplay1=1; end
            if mod(ktm,pauseordisplay1)==0 disply=1;end ;disply=0;
            a=accln_Post_rotation_Coriolis_off(rcomo,vo ,ank2ank,swm,c,mub,flag0,1);
            disply=0;
            k1 = h *a;
            l1 = h *vo;
            k2 = h *accln_Post_rotation_Coriolis_off( rcomo+ l1./2, vo + k1./2,ank2ank,swm,c,mub,flag0,2);
            l2 = h *(vo + k1./2 );
            k3 = h *accln_Post_rotation_Coriolis_off( rcomo+ l2./2, vo + k2./2,ank2ank,swm,c,mub,flag0,3);
            l3 = h *(vo + k2./2 );

 
            dispNplot= mod(ktm,pauseordisplay_interval)==0; 

            if dispNplot disply2=0;else disply2=0; end
            k4 = h *accln_Post_rotation_Coriolis_off( rcomo+ l3, vo + k3,ank2ank,swm,c,mub,flag0,4);
            disply2=0;
            l4 = h *(vo + k3 );
            vn = vo + 1/6 *(k1 + 2*k2 + 2*k3 +  k4); 
            rn = rcomo  + 1/6 *(l1 + 2*l2 + 2*l3 +  l4); 







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





            % Controlling EDGE Loops
            meanedgedelktm=10;
            if rn(1)>APmusclegoalbound & rcomo(1)<APmusclegoalbound
                randloopdirectionchooseswitchgoingUP=rand>=0.5;%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                edgedelktm=fix(abs(rand*20-10+meanedgedelktm));%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                runtillthisktmATedge=ktm+edgedelktm;
            elseif rn(1)<-APmusclegoalbound & rcomo(1)>-APmusclegoalbound
                randloopdirectionchooseswitchgoingDN=rand>0.8;%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                edgedelktm=fix(abs(rand*20-10+meanedgedelktm));%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                runtillthisktmATedge=ktm+edgedelktm;
            end

            edgerang=[0 0 1];
            if randloopdirectionchooseswitchgoingUP==1&randloopdirectionchooseswitchgoingDN==1 edgerang=[1 1 0];end
            if randloopdirectionchooseswitchgoingUP==1&randloopdirectionchooseswitchgoingDN==0 edgerang=[1 0 0];end
            if randloopdirectionchooseswitchgoingUP==0&randloopdirectionchooseswitchgoingDN==1 edgerang=[0 1 0];end
            if randloopdirectionchooseswitchgoingUP==0&randloopdirectionchooseswitchgoingDN==0 edgerang=[0 0 0];end

            if randloopdirectionchooseswitchgoingUP==1 uparedgerang=[0 1 0];else  uparedgerang=[0 0 1];end








            %Below for Velocity Threshold Corrections Control when calculating lateral muscle torques
            rap=rcomo(1);
            rml=rcomo(2);
            vap=vo(1);
            vml=vo(2);
            usetogetsmallvthresh=10;
            usetogetbigvthresh=10;
            if ~correcting

                    if rml<0 & vap>0 % right up
                        vMLtrigthresh= rand/usetogetsmallvthresh;   %NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                        if vml<-vMLtrigthresh correcting=1;runtillthisktm=ifcorrectdothis(ktm);end
           
                    elseif rml<0 & vap<0 % right down
                        vMLtrigthresh= rand/usetogetbigvthresh;  %NOISE SOURCE  %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                        if vml>vMLtrigthresh correcting=1;runtillthisktm=ifcorrectdothis(ktm);end
            

                    elseif rml>0 & vap>0 % left up
                        vMLtrigthresh= rand/usetogetbigvthresh;  %NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                        if vml<-vMLtrigthresh correcting=1;runtillthisktm=ifcorrectdothis(ktm);end
            

                    elseif rml>0 & vap<0 % left down
                        vMLtrigthresh= rand/usetogetsmallvthresh; %NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                        if vml>vMLtrigthresh correcting=1;runtillthisktm=ifcorrectdothis(ktm);end
           
                    end

                    ncorrkting=ncorrkting+correcting;
            end

             if ktm>runtillthisktm correcting=0;end

            





















            % DATA STORAGE
            
                sr(ktm,:)=rcomo(:);
                sv(ktm,:)=vo(:);
         
br(ktm)=brake;


















            % DISPLAY AND PLOTS
            if ktm>forcedisplaylim pauseordisplay_interval=1;end
           
            if dispNplot 
                pod=pauseordisplay_interval-1;
            
             mod9=mod(somecounter,byby^2);fix9=fix(somecounter/byby^2);
            if somecounter/byby^2-fix(somecounter/byby^2)~=0 figno=fix9+1;splotno=mod9;
            else figno=fix9;splotno=byby^2;end


            if plotrealtime
                
                
                
                  
                         
                    
                   
                    
                    
                    figure(figno);

                    if length(slopetaurange)==1 extratit=['   ' 'ktm:' num2str(ktm)];
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
                                                                    
                    if length(slopetaurange)==1 
                    
                    plot(sr(lsr,2),sr(lsr,1),'bo','MarkerSize',2)
                    end
                                            
                     if ~onceplot  
                         plot([-0.1 0.1],[0 0],'color',mc('orange'));plot([0 0],[-0.1 0.1],'color',mc('orange'))
                         plot(rcomub0(2),rcomub0(1),'bo','MarkerSize',markersizenum,'LineWidth', linewdth)
                         
                         title(['     Post    ' 'K_S' '= ' num2str(alphatauml) '   ' 'K_{ED}' '= ' num2str(slopebeta)]);
                         ylabel('ap (mt)');xlabel('ml (mt)');
                         axis([-0.1,0.1, -0.1, 0.1]*(1+.05));%1=x=ap  2=y=ml
                        set(gca,'XDir','rev');
                                            if length(slopetaurange)==1
                            plot([0+centerborder 0+centerborder], [-0.1 0.1], 'y')
                            plot([0-centerborder 0-centerborder], [-0.1 0.1], 'y')
                            hv(1,APmusclegoalbound);
                            hv(1,-APmusclegoalbound);
                            hv(2,MLmusclebound);
                            hv(2,-MLmusclebound);
                         end
                                                       onceplot=1;
                     end
                     
     
            end

           
            
            end























            outofbound=0.1;
            if abs(sr(ktm-1,1))>outofbound | abs(sr(ktm-1,2))>outofbound 
                shitfallen=1;tfall=ktm*dt;


             ktm=4000;t=tstop; 

            
            end 

                                                  




















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
                      plot(sr(end,2),sr(end,1),'r*','MarkerSize',markersizenum,'LineWidth', linewdth) 
                  else 
                      plot(sr(end,2),sr(end,1),'g*','MarkerSize',markersizenum,'LineWidth', linewdth)
                  end
            plot(sr(:,2),sr(:,1),'k-')
            end
            
            
            
            
end %END ITERRUN























    if length(slopetaurange)>1
        
         donematrix(ksrow,kdecol)=0.1;
         fD=1200;hD=figure(fD);
         surf(slopebetarange,slopetaurange,donematrix)%changed x<-->y in surf
         set(gca,'YTick',[0:100:1000])
         set(gca,'XTick',[0:10])  
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
        surf(slopebetarange,slopetaurange,P_notfall_matrix_iter)%changed x<-->y in surf
         set(gca,'YTick',[0:100:1000])
         set(gca,'XTick',[0:10])  
         zlim([0 1.5])
         zlabel('P')
         if ispc set(gcf,'position',[616 516 560 420]);end
         
         
        
        t_fall_matrix1=t_fall_matrix/max(max(t_fall_matrix));
        fT=500;hT=figure(fT);
        surf(slopebetarange,slopetaurange,t_fall_matrix)%changed x<-->y in surf
        set(gca,'YTick',[0:100:1000])
        set(gca,'XTick',[0:10])  
        zlim([0 20])
        zlabel('T')
        if ispc 

set(gcf,'position',[24 1 560 420])
        end
       
        
         disp('#################################################')
    end
%     pause



    end %END OF KED SEARCH
end%END OF KS SEARCH




















































if disply2
    for j=1:2
    figure(98000+j)
    for ii=1:12 subplot(4,3,ii); hv(1,0);end
    end 
end
    
    















    
                




if length(PARAMRANGE1)>1
pathfig='C:\Users\abakshi\Desktop\simulation to bci comp\figs\';
pathfig='C:\Users\abakshi\Desktop\Simul Matlab\CurrentFigs2save\'
saveas(hP,[pathfig 'P.fig'])
saveas(hT,[pathfig 'T.fig'])
hvar=figure(1300);
surf(slopebetarange,slopetaurange,var_t_fall_matrix)%changed x<-->y in surf
 set(gca,'YTick',[0:100:1000])
 set(gca,'XTick',[0:10])  
 zlim([0 20])
 zlabel('Variance')
  if ispc set(gcf,'position',[1963 507 560 420]);end
saveas(hvar,[pathfig 'Var.fig'])
end


if length(PARAMRANGE1)==1
figure(54684)

hold on
plot(sv(:,2),'r*')
plot(sr(:,2),'go')
 grid minor
 hv(1,MLmusclebound)
 hv(1,-MLmusclebound)
hv(1,Vbrakethresh)
 hv(1,-Vbrakethresh)
 ibr1=find(br==1)';
 plot(ibr1,sv(ibr1,2),'mo')
 plot(ibr1,sr(ibr1,2),'co')
end