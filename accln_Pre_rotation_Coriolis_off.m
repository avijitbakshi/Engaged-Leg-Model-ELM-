
function accl=accln_for_no_coriolis_PRE_TRACE(rcomub,v ,ank2ank,swmuscl,c,mub,flag,RKcall);

global ktm mset negmset mrandix velnow4beta velpast4beta mold mnew cold cnew ktmgt4 swcori disply disply2 centerborder pauseordisplay2 generaltempstorageshelf betaloadedengleg disenglegload c_lt_begin c_lt_end c_rt_begin c_rt_end cend flagcom flagactiveleg alphatauml slopebeta gravity goal APmusclegoalbound MLmusclebound correcting vMLtrigthresh runtillthisktm runtillthisktmATedge randloopdirectionchooseswitchgoingUP randloopdirectionchooseswitchgoingDN rando randnu switch_x_or_z_tork swmlV lots_of_noise amp amplow Ataumid edgeloopcontrol midVcontrol APswgoal fdeltaktm APedgestiffens MLedgestiffens APvthresh4goal brake
clear tauwt taudiseng tautot






%%%SWITCHES INITIALIZATIONS DEFINITIONS SETTINGS AND PARAMS
%SWITCHES
%PARAM SETTING / DIALS
%MATLAB SETTINGS
%CONSTANTS
%INITILIZERS
%VARIABLE VALUES






%SWITCHES
tangent_dirn_switch=1;
radial_dirn_switch=1;
noiseonoffswitch=1;%NOISE SOURCE AAND SWITCH




 
ibyi=1;
COM_dependent_activeleg=1;
noiseintauxz=1;





% %PARAM SETTING / DIALS
if switch_x_or_z_tork
    if noiseintauxz %NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    ampmean= 170; amp= ampmean +  rando(2)*40-20;
    faclow=(rando(3)*100+1);

    amplow=ampmean/faclow;

    else
        amp=170;amplow=amp/3.5;%JIM MIRACLE SHOW VALUES 23MAY2014  

    end
                % ctz=1;
end


 if ~lots_of_noise 
     for j=1:6 addn(j)=0;end % muscle has noice via randn, but not add noise to multiplying factors
 else %add noise to multiplying factors
     for j=1:6 rand_1to1(j)=(rando(10+j) -0.5)/0.5;end
     f=5;
     j=1;addn(j)=rand_1to1(j)*20*f; 
     j=2;addn(j)=rand_1to1(j)*20*f; 
     j=3;addn(j)=rand_1to1(j)*1*f; 
     j=4;addn(j)=rand_1to1(j)*0.2*f; 
     j=5;addn(j)=rand_1to1(j)*10*f; 
     j=6;addn(j)=rand_1to1(j)*10*f; 
 end
facclnoiseap=80+addn(1);facclnoiseml=80+addn(2);facclnoisez=0+addn(3);%simulation surf with values =10%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
fconstbiasnoise=0.1+addn(4);%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
ftaumusnoisemlaffectap=25*2+addn(5);ftaumusnoiseapaffectml=25*1+addn(6);%NOISE SOURCE
ftaumusnoise=[ftaumusnoiseapaffectml ftaumusnoisemlaffectap fconstbiasnoise];%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


rap=rcomub(1);rml=rcomub(2);
vap=v(1);
alphatauap=0;
KsEdgeEnhancingfactor=1;Vbrakethresh=0.1;
increaseKsfactoratEdgeTo=2;

 if abs(rap)>APmusclegoalbound & APedgestiffens==1 & abs(vap) > Vbrakethresh  % AP edge regions

     if rap>0 & vap>Vbrakethresh
         KsEdgeEnhancingfactor=increaseKsfactoratEdgeTo;
         brake=1;
     elseif rap<0 & vap<-Vbrakethresh
         KsEdgeEnhancingfactor=increaseKsfactoratEdgeTo;
         brake=1;
     end
 end
 if abs(rml)>MLmusclebound & MLedgestiffens==1
     alphatauap=50;
 end


disengamplifypercentfac=0/100;%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@











%DISPLAY SETTINGS
nr=4;nc=3;l=80;%j=2;






%CONSTANTS
x=[1 0 0];y=[0 1 0];z=[0 0 1];
w=mub*gravity;
omega=pi/3;
mft=0.6;
anklt=(ank2ank/2).*y;ankrt=(ank2ank/2).*(-y);
fn=w.*z;







%INITILIZERS

taumuscle=[0 0 0];taumusclecyl=[0 0 0];
accllft=0;acclrt=0;acclcentposn=0;
taumusgoal=[0 0 0];
tauxz=[0 0 0];
taumuslateralboundV=[0 0 0];





















%VARIABLE VALUES
hcomub=norm(rcomub);%used for I=MR^2 and centriaccln during flag=0 ie Cylinder state
ncapub=rcomub./hcomub;

%compute coriolis
fcori=(-2*mub*ibyi*swcori).*cross(omega.*z,v);

%compute disengaged variables 
if RKcall == 1 %No change in beta for Rkcall>1 RKcall is RungeKutta call
          if COM_dependent_activeleg==1
              flagactiveleg=flag;
              if flagactiveleg~=0 beta_com_side_leg=betaf(abs(rcomub(2)));
                  betaloadedengleg=beta_com_side_leg; 
                  disenglegload=(1-betaloadedengleg).*fn;
              end%%ACTIVATE IF DISENG/ENG WRT COM
              flagcom=flag;
          else
              legcutoff=55/100;flagcom=flag;

                    if flagcom==1          disp(['_______________    LEFT COM    ' num2str(ktm) '   ________________________']);
                    elseif  flagcom==-1    disp(['________________   RIGT COM  ' num2str(ktm) '   _________________________']);
                    else                   disp(['_________________  CENTER COM  '   num2str(ktm) ' ________________________']);
                    end

                    beta_com_side_leg=betaf(abs(rcomub(2))); 
                    if  beta_com_side_leg>legcutoff flagactiveleg=flagcom; betaloadedengleg=beta_com_side_leg; disenglegload=(1-betaloadedengleg).*fn;
                    elseif beta_com_side_leg < 1-legcutoff  flagactiveleg=-flagcom; betaloadedengleg=1-beta_com_side_leg; disenglegload=(1-betaloadedengleg).*fn;
                    else flagactiveleg=0;%%No eng/diseng concept
                    end
          end  
end
    


    
   
    
    
    
    
    
    
    









    
if flagactiveleg == 1 %left
   
                    if disply==1 disp(['***************  LEFT side eng    ' num2str(ktm) '   **************************']);% disp(['  LEFT side eng'  ]);
                    end
                %     RKcall
                    rcomlt= rcomub-anklt;rcomfromengagedank=rcomlt;hcomlt=norm(rcomlt);ncaplt= rcomlt./hcomlt; %vector from lt leg ankle to com
                                         pos_of_diseng_fn=ank2ank.*(-y)+dot(rcomub,x).*x; % in fn not w is put bcos mass is cancelled in next eqn
                    taudiseng=cross(pos_of_diseng_fn,disenglegload)*(1+disengamplifypercentfac);
                    tauwt= cross(rcomfromengagedank,w.*(-z));
                    taucori=cross(rcomlt,fcori);

                    if swmuscl taumuscle=muscle2(rcomub,v,flagactiveleg,alphatauap,ftaumusnoise,noiseonoffswitch,KsEdgeEnhancingfactor);%swmuscl.* taumusub;
                    end
                    if APswgoal taumusgoal=muscle2goal(rcomub,v);end

                    if swmlV taumuslateralboundV=musclelateralboundV(rcomub,v);
                    end
                    if switch_x_or_z_tork if v(1)<0 tzC=amp;else tzC=amplow;end; tauxz=-sign(v(1))*tzC*norm(v)*x;end;
                                  tautot=tauwt  + swcori*taucori + swmuscl*taumuscle + taudiseng  +APswgoal*taumusgoal + swmlV*taumuslateralboundV +switch_x_or_z_tork*tauxz;

                    if disply2
                               disp(['***@@@@@@@@@@@@@@@@@@  LEFT side eng    ' num2str(ktm) '   @@@@@@@@@@@@@@@@@@@@@@******']);
                      

                               for j=1:2
                               figure(98000+j)
                               if j==1 set(gcf,'position',[576 524 1463 462]);elseif j==2 set(gcf,'position',[578 -1 1463 462]);end
                               clr=[0 1 0];
                               subplot(nr,nc,1);plot(ktm,taudiseng(j),'*','color',clr);hold on;title('taudiseng');
                               subplot(nr,nc,2);plot(ktm,tauwt(j),'*','color',clr);hold on;title('tauwt');
                               subplot(nr,nc,4);plot(ktm,taucori(j),'*','color',clr);hold on;title('taucori');ylim([-l l])
                                       subplot(nr,nc,9);plot(ktm,taumuscle(j),'*','color',clr);hold on;title('taumuscle');ylim([-l l])
                                     subplot(nr,nc,12);plot(ktm,tautot(j),'*','color',clr);hold on;title('tautot');ylim([-l l])
                               subplot(nr,nc,3);plot(ktm,taudiseng(j)+tauwt(j),'*','color',clr);hold on;title('taudiseng+tauwt');ylim([-l l])
                                 subplot(nr,nc,7);plot(ktm,taumusgoal(j),'*','color',clr);hold on;title('taumusgoal');ylim([-l l])
                                       subplot(nr,nc,11);plot(ktm,taumuslateralboundV(j),'*','color',clr);hold on;title('taumuslateralboundV');ylim([-l l])
                               end
                    end
           
                     
                    alpha=( tautot )./(mub*hcomlt^2) ;%mass taken off below from muscle factor also
                    vtan=v-dot(v,ncaplt).*ncaplt;
                    accllft=cross(alpha,rcomlt)*tangent_dirn_switch+(dot(vtan,vtan)/hcomlt).*(-ncaplt)*radial_dirn_switch;
    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
elseif flagactiveleg == -1 %right
    
                   if disply==1 disp(['###############   RIGT side eng  ' num2str(ktm) '   ##########################']); %disp([ '  RIGT side eng ']);
                   end
                %    RKcall
                    rcomrt= rcomub-ankrt;rcomfromengagedank=rcomrt;hcomrt=norm(rcomrt);ncaprt= rcomrt./hcomrt; %vector from lt leg ankle to com
                                                            pos_of_diseng_fn=ank2ank.*y + dot(rcomub,x).*x;  % in fn not w is put bcos mass is cancelled in next eqn
                     taudiseng=cross(pos_of_diseng_fn,disenglegload)*(1+disengamplifypercentfac);
                     tauwt=cross(rcomfromengagedank,w.*(-z));
                     taucori=cross(rcomrt,fcori);

                     if swmuscl taumuscle=muscle2(rcomub,v,flagactiveleg,alphatauap,ftaumusnoise,noiseonoffswitch,KsEdgeEnhancingfactor);%swmuscl.* taumusub;
                    end
                     if APswgoal taumusgoal=muscle2goal(rcomub,v);end
                     if swmlV taumuslateralboundV=musclelateralboundV(rcomub,v);
                     end
                             if switch_x_or_z_tork if v(1)>0 tzC=amp;else tzC=amplow;end; tauxz=-sign(v(1))*tzC*norm(v)*x;%disp([tzC tauxz]);
                     end;
                                    tautot=tauwt  + swcori*taucori + swmuscl*taumuscle + taudiseng  +APswgoal*taumusgoal + swmlV*taumuslateralboundV +switch_x_or_z_tork*tauxz;
                     
                     if disply2
                                 disp(['###&&&&&&&&&&&&   RIGT side eng  ' num2str(ktm) '   &&&&&&&&&&&&&&&&&#########']);
                               disp(flagactiveleg)
                              
                               for j=1:2
                               figure(98000+j)
                               if j==1 set(gcf,'position',[576 524 1463 462]);elseif j==2 set(gcf,'position',[578 -1 1463 462]);end
                               clr=[1 0 0];
                               subplot(nr,nc,1);plot(ktm,taudiseng(j),'*','color',clr);hold on;title('taudiseng');
                               subplot(nr,nc,2);plot(ktm,tauwt(j),'*','color',clr);hold on;title('tauwt');
                               subplot(nr,nc,4);plot(ktm,taucori(j),'*','color',clr);hold on;title('taucori');ylim([-l l])
        %                        subplot(nr,nc,10);plot(ktm,taucori(j),'*','color',clr);hold on;title('taucori');
        %                        subplot(nr,nc,10);plot(ktm,tauxz(j),'o','color',clr);hold on;title('*taucori otauxz');ylim([-l l])
                               subplot(nr,nc,9);plot(ktm,taumuscle(j),'*','color',clr);hold on;title('taumuscle');ylim([-l l])
        %                        subplot(nr,nc,5);plot(ktm,tauxz(j),'*','color',clr);hold on;title('tauxz');ylim([-l l]);ylim([-l l])
                               subplot(nr,nc,12);plot(ktm,tautot(j),'*','color',clr);hold on;title('tautot');ylim([-l l])
                               subplot(nr,nc,3);plot(ktm,taudiseng(j)+tauwt(j),'*','color',clr);hold on;title('taudiseng+tauwt');ylim([-l l])
        %                        subplot(nr,nc,6);plot(ktm,taucori(j)+tauxz(j),'*','color',clr);hold on;title('taucori+tauxz');ylim([-l l])
                               subplot(nr,nc,7);plot(ktm,taumusgoal(j),'*','color',clr);hold on;title('taumusgoal');ylim([-l l]*2)
        %                        subplot(nr,nc,8);plot(ktm,taumuslateralbound(j),'*','color',clr);hold on;title('taumuslateralbound');ylim([-l l])
                        %                        subplot(nr,nc,11);plot(ktm,taumuslateralboundintermittent(j),'*','color',clr);hold on;title('taumuslateralboundintermittent');ylim([-l l])
                               subplot(nr,nc,11);plot(ktm,taumuslateralboundV(j),'*','color',clr);hold on;title('taumuslateralboundV');ylim([-l l])
                               end
                     end
                    alpha=( tautot )./(mub*hcomrt^2) ;%mass taken off below from muscle factor also
                    vtan=v-dot(v,ncaprt).*ncaprt;
                    acclrt=cross(alpha,rcomrt)*tangent_dirn_switch+(dot(vtan,vtan)/hcomrt).*(-ncaprt)*radial_dirn_switch; 








                    
                    
                    
                    
                    
                    
                    
                    
                    




elseif flagactiveleg == 0 
    
                    if disply==1 disp(['oooooooooooooooo   CENTER  '   num2str(ktm) '   ooooooooooooooooooooooooooo']); %disp([ '  CENTER  ' ]);
                    end
                                                                                         tauwt=cross(rcomub,w.*(-z));
                     taucori=cross(rcomub,fcori);
                      if swmuscl taumusclecyl=muscle2(rcomub,v,flagactiveleg,alphatauap,ftaumusnoise,noiseonoffswitch,KsEdgeEnhancingfactor);end
                      if switch_x_or_z_tork if v(1)>0 tzC=amp;else tzC=amp;end; tauxz=-sign(v(1))*tzC*norm(v)*x;%disp([tzC tauxz]);
                      end;
                     tautot=tauwt + swcori.*taucori + swmuscl*taumusclecyl + switch_x_or_z_tork*tauxz;
                    if disply2
                        disp(['o$$$$$$$$$   CENTER  '   num2str(ktm) '   $$$$$$$$$$oooooooooo']);
                                            for j=1:2
                       figure(98000+j)
                       clr=[0 0 1];
%                        subplot(nr,nc,1);plot(ktm,taudiseng(j),'*','color',clr);hold on;title('taudiseng');
                       subplot(nr,nc,2);plot(ktm,tauwt(j),'*','color',clr);hold on;title('tauwt');
                       subplot(nr,nc,4);plot(ktm,taucori(j),'*','color',clr);hold on;title('taucori');ylim([-l l])
%                        subplot(nr,nc,10);plot(ktm,taucori(j),'*','color',clr);hold on;title('taucori');
%                        subplot(nr,nc,10);plot(ktm,tauxz(j),'o','color',clr);hold on;title('*taucori otauxz');ylim([-l l])
                       subplot(nr,nc,9);plot(ktm,taumusclecyl(j),'*','color',clr);hold on;title('taumuscle');ylim([-l l])
%                        subplot(nr,nc,5);plot(ktm,tauxz(j),'*','color',clr);hold on;title('tauxz');ylim([-l l])
                       subplot(nr,nc,12);plot(ktm,tautot(j),'*','color',clr);hold on;title('tautot');ylim([-l l])
%                         subplot(nr,nc,6);plot(ktm,taucori(j)+tauxz(j),'*','color',clr);hold on;title('taucori+tauxz');ylim([-l l])
                       end
                    end
                    
                    alpha=(tautot)./(mub*hcomub^2);
                    vtan=v-dot(v,ncapub).*ncapub; 
                    acclcylncrossalpha=cross(alpha,rcomub);
                    acclcylncentri=(dot(vtan,vtan)/hcomub).*(-ncapub);
                    acclcyln=acclcylncrossalpha + acclcylncentri;
                    acclcylncori=[0 0 0];
                    acclcentposn=acclcyln+acclcylncori;
%                     pause
end % end of flagactiveleg











accl=accllft+acclrt+acclcentposn;
percentnoise=[facclnoiseap*randnu(1) facclnoiseml*randnu(2) facclnoisez*randnu(3)];%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
acclnoise=(percentnoise./100).*accl;%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
accl=accl+acclnoise*noiseonoffswitch;











if disply==1
%     display(['ktm =' num2str(ktm)]);
    if flagactiveleg==0

    else

            betadiseng=(1-betaloadedengleg);
            display(['betadiseng        =  ' num2str(betadiseng)])
            display(['disenglegload     =  ' num2str(disenglegload)])
            display(['pos_of_diseng_fn  =  ' num2str(pos_of_diseng_fn)])
            display(['rcom              =  ' num2str(rcomub)])
            display(['rcomfromengagedank=  ' num2str(rcomfromengagedank)])
%             w.*(-z)
            disp('Taus: ---');
            display(['taudiseng= ' num2str( taudiseng)])
            display(['tauwt    = ' num2str( tauwt)])
            tauengmus=taumusub;
            display(['tauengmus= ' num2str( tauengmus)])
            display(['taucori  = ' num2str(  taucori)])
            disp('__________+R(-)/-L(+)   +U/-D');
%             disp('         +R(-)/-L(+)   +U/-D');   
            display(['tautot   =' num2str(tautot)])
            disp('       ')
    end
end
















return
end % end of acclnquick
































































%--------------------------------------------------------------------------
%--------------SUB FUNCTION---------------------------------------------

function  taumusub2=muscle2(rcomu,v ,flagactivelegx,alphatauap,ftaumusnoise,noiseonoffswitch,KsEdgeEnhancingfactor);
global ktm mset negmset mrandix velnow4beta velpast4beta mold mnew cold cnew ktmgt4 swcori disply disply2 centerborder pauseordisplay2 generaltempstorageshelf betaloadedengleg disenglegload c_lt_begin c_lt_end c_rt_begin c_rt_end cend flagcom flagactiveleg alphatauml slopebeta gravity goal APmusclegoalbound MLmusclebound correcting vMLtrigthresh runtillthisktm runtillthisktmATedge randloopdirectionchooseswitchgoingUP randloopdirectionchooseswitchgoingDN rando randnu switch_x_or_z_tork swmlV lots_of_noise amp amplow Ataumid edgeloopcontrol midVcontrol APswgoal fdeltaktm APedgestiffens MLedgestiffens APvthresh4goal brake
UPFRONTSIGN_AFTER_TAUENGCODE_UNDSTOOD=[1 -1];%[LEFT RIGHT] HERE



sign_tauplate2flap_lt_rt=[1 -1];
                                if flagactivelegx==1 %left +y
                                    tauplate2ftml=-alphatauml*KsEdgeEnhancingfactor*rcomu(1);% 1  proj = ap here (diff from expt 2)
                                     tauplate2ftap=UPFRONTSIGN_AFTER_TAUENGCODE_UNDSTOOD(1)*alphatauap*abs(rcomu(2));%put on 25May2104

                                elseif flagactivelegx==-1
                                    tauplate2ftml=-alphatauml*KsEdgeEnhancingfactor*rcomu(1);% 1  proj = ap here (diff from expt 2)
                                      tauplate2ftap=UPFRONTSIGN_AFTER_TAUENGCODE_UNDSTOOD(2)*alphatauap*abs(rcomu(2));%put on 25May2104
                                elseif flagactivelegx==0
                                tauplate2ftml=-2*alphatauml*rcomu(1);

                                tauplate2ftap=0*alphatauap*rcomu(1);
                                end
                                
               
tauplate2ft=[tauplate2ftap tauplate2ftml 0];% ap is x in old coord sys
percentnoiseintau=[ftaumusnoise(1)*randnu(4) abs(ftaumusnoise(2)*randnu(5)) randnu(6)];%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
fconstbias=ftaumusnoise(3);constbiasnoise=fconstbias.*[randnu(7) randnu(8) randnu(9)];%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

taumusnoise=(percentnoiseintau./100).*tauplate2ft;%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


taumusub2=tauplate2ft+(taumusnoise+constbiasnoise)*noiseonoffswitch;%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



return
end
































%--------------------------------------------------------------------------
%--------------SUB FUNCTION---------------------------------------------

function betaengagedlegload=betaf(absyproj);
%%%%Gives load under engaged leg
global ktm mset negmset mrandix velnow4beta velpast4beta mold mnew cold cnew ktmgt4 swcori disply disply2 centerborder pauseordisplay2 generaltempstorageshelf betaloadedengleg disenglegload c_lt_begin c_lt_end c_rt_begin c_rt_end cend flagcom flagactiveleg alphatauml slopebeta gravity goal APmusclegoalbound MLmusclebound correcting vMLtrigthresh runtillthisktm runtillthisktmATedge randloopdirectionchooseswitchgoingUP randloopdirectionchooseswitchgoingDN rando randnu switch_x_or_z_tork swmlV lots_of_noise amp amplow Ataumid edgeloopcontrol midVcontrol APswgoal fdeltaktm APedgestiffens MLedgestiffens APvthresh4goal brake
if swcori ==0 | ktmgt4==0

betaengagedlegload=slopebeta*absyproj+0.5;
elseif swcori==1 & ktmgt4==1

betaengagedlegload=slopebeta*absyproj+0.5;

   
end%swc
return
end





































function  taumusgoal=muscle2goal(rcomu,v);
global ktm mset negmset mrandix velnow4beta velpast4beta mold mnew cold cnew ktmgt4 swcori disply disply2 centerborder pauseordisplay2 generaltempstorageshelf betaloadedengleg disenglegload c_lt_begin c_lt_end c_rt_begin c_rt_end cend flagcom flagactiveleg alphatauml slopebeta gravity goal APmusclegoalbound MLmusclebound correcting vMLtrigthresh runtillthisktm runtillthisktmATedge randloopdirectionchooseswitchgoingUP randloopdirectionchooseswitchgoingDN rando randnu switch_x_or_z_tork swmlV lots_of_noise amp amplow Ataumid edgeloopcontrol midVcontrol APswgoal fdeltaktm APedgestiffens MLedgestiffens APvthresh4goal brake

if abs(rcomu(1))<APmusclegoalbound & abs(v(1))<APvthresh4goal
    taumusgoal=[0 goal*30 0];

else taumusgoal=[0 0 0];
end 

end









































function  taumuslateralboundV=musclelateralboundV(rcomu,v);
global ktm mset negmset mrandix velnow4beta velpast4beta mold mnew cold cnew ktmgt4 swcori disply disply2 centerborder pauseordisplay2 generaltempstorageshelf betaloadedengleg disenglegload c_lt_begin c_lt_end c_rt_begin c_rt_end cend flagcom flagactiveleg alphatauml slopebeta gravity goal APmusclegoalbound MLmusclebound correcting vMLtrigthresh runtillthisktm runtillthisktmATedge randloopdirectionchooseswitchgoingUP randloopdirectionchooseswitchgoingDN rando randnu switch_x_or_z_tork swmlV lots_of_noise amp amplow Ataumid edgeloopcontrol midVcontrol APswgoal fdeltaktm APedgestiffens MLedgestiffens APvthresh4goal brake





    
    
Ataumid=30/3;
Vartaumid=10;
fac4lesstauif2mid=2;
Atauedge=3;
Vartauedge=1; 





 
vap=v(1);
vml=v(2);
rap=rcomu(1);
rml=rcomu(2);
taux=[0 0 0]; 
if midVcontrol
   if correcting
   
           if abs(rap)<APmusclegoalbound % center region
               taux=[0 0 0]; 
               Atau=Ataumid+ Vartaumid*((rando(4)-1/2)/0.5);%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
               Atau=abs(Atau);

               
               
               
               if rml<0 % right side   
                   if vml<0 %going more -ve => away to more right
                       taux=-Atau*heaviside(runtillthisktm-ktm); taux=[taux 0 0]; 

                 elseif vml>0 
                      
                   end
               elseif rml>0 %left side
                   if vml>0 %going more +ve => away to left
                      taux=Atau/fac4lesstauif2mid*heaviside(runtillthisktm-ktm); taux=[taux 0 0];
                   elseif vml<0 
                      
                   end
               end
               

               
               
           end
   end
end

   
   
  if edgeloopcontrol 
           if abs(rap)>APmusclegoalbound % AP edge regions
               taux=[0 0 0]; 
               
               Atau=(Atauedge+Vartauedge*((rando(5)-1/2)/0.5));%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                           if vap>0&rap>0&randloopdirectionchooseswitchgoingUP==1 % up
                               taux=[-Atau*heaviside(runtillthisktmATedge-ktm) 0 0];

                           elseif vap<0&rap<0&randloopdirectionchooseswitchgoingDN==1 %down
                                taux=[Atau*heaviside(runtillthisktmATedge-ktm) 0 0];  
                           elseif vap>0&rap<0&randloopdirectionchooseswitchgoingDN==0&randloopdirectionchooseswitchgoingUP==0 % =0 => lsst top edge was normal=> went rightward loop for ging up=> now in down edge i send it to big leftward loop.
                                taux=[-Atau*heaviside(runtillthisktmATedge-ktm) 0 0]; 
                           end
                           

               end

    end
               
              
               

               
    taumuslateralboundV=taux;      
end 