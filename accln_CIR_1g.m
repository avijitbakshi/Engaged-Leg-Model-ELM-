function accl=accln__CIRCULAR_selectKS(rcomub,v ,ank2ank,swmuscl,c,mub,flag,RKcall);%sw=switch
global ktm mset negmset mrandix velnow4beta velpast4beta mold mnew cold cnew ktmgt4 swcori disply centerborder pauseordisplay2 generaltempstorageshelf betaloadedengleg disenglegload c_lt_begin c_lt_end c_rt_begin c_rt_end cend flagcom flagactiveleg alphatauml slopebeta gravity alphatauRL vtannedswitch quadrantsign
clear tauwt taudiseng tautot 
noiseonoffswitch=0;
facclnoiseap=80;facclnoiseml=80;facclnoisez=0;%simulation surf with values =10
fconstbiasnoise=0.1;
ftaumusnoisemlaffectap=20*1;ftaumusnoiseapaffectml=20*1;ftaumusnoise=[ftaumusnoiseapaffectml ftaumusnoisemlaffectap fconstbiasnoise];
disengamplifypercentfac=0/100;
x=[1 0 0];y=[0 1 0];z=[0 0 1];taumusub=[0 0 0];newmus=1;
alphatauap=14;75;60;
w=mub*gravity;omega=pi/3;mft=0.6;anklt=(ank2ank/2).*y;ankrt=(ank2ank/2).*(-y);
hcomub=norm(rcomub);ncapub=rcomub./hcomub;costh=rcomub(3)./hcomub;sinth=sqrt(rcomub(1)^2+rcomub(2)^2)/hcomub;sinphi= rcomub(2)/(hcomub*sinth);cosphi=rcomub(1)/(hcomub*sinth);
ncapltleg=ncapub;ncaprtleg=ncapub; % just change this defn incase you want diff from parallel legs paradigm 
normaltorcomcap=[costh*cosphi costh*sinphi -sinth];

ibyi=1;fcori=(-2*mub*ibyi*swcori).*cross(omega.*z,v);
fn=w.*z;
accllft=0;acclrt=0;acclcentposn=0;
KsEdgeEnhancingfactor=1;
quadrantsignselectforRL=-1;


COM_dependent_activeleg=1;
if RKcall == 1 %No change in beta for Rkcall>1
          if COM_dependent_activeleg==1
              flagactiveleg=flag;if flagactiveleg~=0 beta_com_side_leg=betaf(abs(rcomub(2)));betaloadedengleg=beta_com_side_leg; disenglegload=(1-betaloadedengleg).*fn;end%%ACTIVATE IF DISENG/ENG WRT COM
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
    
 
amp=50;ctz=1;x_or_z_tork=1;
    if swcori==0 x_or_z_tork=0;end
if flagactiveleg == 1 %left
                    if disply==1 disp(['***************  LEFT side eng    ' num2str(ktm) '   **************************']);% disp(['  LEFT side eng'  ]);
                    end
                               rcomlt= rcomub-anklt;rcomfromengagedank=rcomlt;hcomlt=norm(rcomlt);ncaplt= rcomlt./hcomlt; %vector from lt leg ankle to com
                     if newmus==1 
                    if quadrantsign==quadrantsignselectforRL
                    taumusub2 =muscle2RL(rcomub,v,flagactiveleg,alphatauap,ftaumusnoise,noiseonoffswitch);taumusub=taumusub2;
                    else
                    taumusub2 =muscle2(rcomub,v,flagactiveleg,alphatauap,ftaumusnoise,noiseonoffswitch,KsEdgeEnhancingfactor);taumusub=taumusub2;
                    end
                    end
                        pos_of_diseng_fn=ank2ank.*(-y)+dot(rcomub,x).*x; % in fn not w is put bcos mass is cancelled in next eqn
                    taudiseng=cross(pos_of_diseng_fn,disenglegload)*(1+disengamplifypercentfac);
                    tauwt= cross(rcomfromengagedank,w.*(-z));
                    taucori=cross(rcomlt,fcori);
                    tautot=tauwt  + swcori*taucori + swmuscl* taumusub + taudiseng ;
                     if x_or_z_tork if v(1)<0 tzC=amp*ctz;else tzC=ctz;end; tauz=tzC*norm(v)*x;tautot=tautot+tauz;end
                    alpha=( tautot )./(mub*hcomlt^2) ;%mass taken off below from muscle factor also
                    if vtannedswitch==1
                    vtan=v-dot(v,ncaplt).*ncaplt;
                    else
                     vtan=v;
                    end
                    accllft=cross(alpha,rcomlt)+(dot(vtan,vtan)/hcomlt).*(-ncaplt);
    
elseif flagactiveleg == -1 %right
                   if disply==1 disp(['###############   RIGT side eng  ' num2str(ktm) '   ##########################']); %disp([ '  RIGT side eng ']);
                   end
                                 rcomrt= rcomub-ankrt;rcomfromengagedank=rcomrt;hcomrt=norm(rcomrt);ncaprt= rcomrt./hcomrt; %vector from lt leg ankle to com
                        if newmus==1 
                        if quadrantsign==quadrantsignselectforRL
                    taumusub2 =muscle2RL(rcomub,v,flagactiveleg,alphatauap,ftaumusnoise,noiseonoffswitch);taumusub=taumusub2;
                    else
                    taumusub2 =muscle2(rcomub,v,flagactiveleg,alphatauap,ftaumusnoise,noiseonoffswitch,KsEdgeEnhancingfactor);taumusub=taumusub2;
                    end
                    end
                     pos_of_diseng_fn=ank2ank.*y + dot(rcomub,x).*x;  % in fn not w is put bcos mass is cancelled in next eqn
                     taudiseng=cross(pos_of_diseng_fn,disenglegload)*(1+disengamplifypercentfac);
                     tauwt=cross(rcomfromengagedank,w.*(-z));
                     taucori=cross(rcomrt,fcori);
                     tautot=tauwt  + swcori.*taucori + swmuscl.* taumusub +  taudiseng;
                      if x_or_z_tork if v(1)>0 tzC=amp*ctz;else tzC=ctz;end; tauz=-tzC*norm(v)*x;tautot=tautot+tauz;end
                    alpha=( tautot )./(mub*hcomrt^2) ;%mass taken off below from muscle factor also
                     if vtannedswitch==1
                     vtan=v-dot(v,ncaprt).*ncaprt;
                     else
                     vtan=v;
                     end
                    acclrt=cross(alpha,rcomrt)+(dot(vtan,vtan)/hcomrt).*(-ncaprt) ;
elseif flagactiveleg == 0 
    disp('center 0')
                    if disply==1 disp(['oooooooooooooooo   CENTER  '   num2str(ktm) '   ooooooooooooooooooooooooooo']); %disp([ '  CENTER  ' ]);
                    end
                                    if newmus==1 
                    if quadrantsign==quadrantsignselectforRL 
                    taumusub2 =muscle2RL(rcomub,v,flagactiveleg,alphatauap,ftaumusnoise,noiseonoffswitch);taumusub=taumusub2;
                    else
                    taumusub2 =muscle2(rcomub,v,flagactiveleg,alphatauap,ftaumusnoise,noiseonoffswitch,KsEdgeEnhancingfactor);taumusub=taumusub2;
                    end
                    end
                                    tauwt=cross(rcomub,w.*(-z));
                     taucori=cross(rcomub,fcori);
                     tautot=tauwt + swcori.*taucori + swmuscl.*taumusub;
                     alpha=(tautot)./(mub*hcomub^2);
                     
 if vtannedswitch==1
     vtan=v-dot(v,ncapub).*ncapub; 
 else
     vtan=v;
 end
                    acclcyln=cross(alpha,rcomub)+(dot(vtan,vtan)/hcomub).*(-ncapub);
                    acclcentposn=acclcyln+fcori./mub;
end % end of flagactiveleg

accl=accllft+acclrt+acclcentposn;
percentnoise=[facclnoiseap*randn facclnoiseml*randn facclnoisez*randn];
acclnoise=(percentnoise./100).*accl;
accl=accl+acclnoise*noiseonoffswitch;


if disply==1

    if flagactiveleg==0

    else

            betadiseng=(1-betaloadedengleg);
            display(['betadiseng        =  ' num2str(betadiseng)])
            display(['disenglegload     =  ' num2str(disenglegload)])
            display(['pos_of_diseng_fn  =  ' num2str(pos_of_diseng_fn)])
            display(['rcom              =  ' num2str(rcomub)])
            display(['rcomfromengagedank=  ' num2str(rcomfromengagedank)])

            disp('Taus: ---');
            display(['taudiseng= ' num2str( taudiseng)])
            display(['tauwt    = ' num2str( tauwt)])
            tauengmus=taumusub;
            display(['tauengmus= ' num2str( tauengmus)])
            display(['taucori  = ' num2str(  taucori)])
            disp('__________+R(-)/-L(+)   +U/-D');
  
            display(['tautot   =' num2str(tautot)])
            disp('       ')
    end
end
sx=0;mllim=40;aplim=30;aplimengwt=130;
subpltixrow=5;subpltixcol=2;tlim=2000;

generaltempstorageshelf(1)=tautot(1);
return
end % end of acclnquick







%--------------------------------------------------------------------------
%--------------SUB FUNCTION---------------------------------------------

function betaengagedlegload=betaf(absyproj);%coming in in meters %%DEFINED ON OCT 1,08
%%%%Gives load under engaged leg
global ktm mset negmset mrandix velnow4beta velpast4beta mold mnew cold cnew ktmgt4 swcori disply centerborder pauseordisplay2 generaltempstorageshelf betaloadedengleg disenglegload c_lt_begin c_lt_end c_rt_begin c_rt_end cend flagcom flagactiveleg alphatauml slopebeta gravity alphatauRL vtannedswitch quadrantsign
if swcori ==0 | ktmgt4==0

betaengagedlegload=slopebeta*absyproj+0.5;
elseif swcori==1 & ktmgt4==1

betaengagedlegload=slopebeta*absyproj+0.5;

   
end%swc
return
end



















function  taumusub2=muscle2RL(rcomu,v ,flagactivelegx,alphatauap,ftaumusnoise,noiseonoffswitch);


global ktm mset negmset mrandix velnow4beta velpast4beta mold mnew cold cnew ktmgt4 swcori disply centerborder pauseordisplay2 generaltempstorageshelf betaloadedengleg disenglegload c_lt_begin c_lt_end c_rt_begin c_rt_end cend flagcom flagactiveleg alphatauml slopebeta gravity alphatauRL vtannedswitch quadrantsign

localswitchRLmover=1;
localswitchFBhelper=1;

UPFRONTSIGN_AFTER_TAUENGCODE_UNDSTOOD=[1 -1];%[LEFT RIGHT] HERE
constantAPtorqueRL=0;
sign_tauplate2flap_lt_rt=[1 -1];
extra=0;thresh2=0.03;
                                if flagactivelegx==1 %left +y  
                                    tauplate2ftml=-alphatauml*rcomu(1)*localswitchRLmover;% 1  proj = ap here (diff from expt 2)
                                    if abs(rcomu(2))>thresh2  constantAPtorqueRL=extra;else constantAPtorqueRL=0;end;%AVIJIT Oct 4 2015  PUT BECOZ NEEDED DRIVE , TO SHOW OSCILL TRACES, ELSE WAS JUST STAYING ERECT AT CENTER OF ML DIMN , ie ML=0
                                     tauplate2ftap=UPFRONTSIGN_AFTER_TAUENGCODE_UNDSTOOD(1)*alphatauRL*abs(rcomu(2))+constantAPtorqueRL;

                                elseif flagactivelegx==-1
                                    tauplate2ftml=-alphatauml*rcomu(1)*localswitchRLmover;% 1  proj = ap here (diff from expt 2)
                                    if abs(rcomu(2))>thresh2 constantAPtorqueRL=extra;else constantAPtorqueRL=0;end;%AVIJIT Oct 4 2015  THIS AND ABOVE SWITH ON ONLY FOR PARAMS WHERE NEED RT RL SWAY
                                    tauplate2ftap=UPFRONTSIGN_AFTER_TAUENGCODE_UNDSTOOD(2)*alphatauRL*abs(rcomu(2))-constantAPtorqueRL;
                                elseif flagactivelegx==0
                                tauplate2ftml=-2*alphatauml*rcomu(1);
                                tauplate2ftap=0*alphatauap*rcomu(1);
                                end
                                
                                

                                thresh1=0.04; pc_KSFB_increase=100;
                                thresh1=0.08; pc_KSFB_increase=100;
                                thresh1=0.03; pc_KSFB_increase=100;

thresh1=0; pc_KSFB_increase=0;%we choose 0 0
                                if abs(rcomu(1))>thresh1 
                                    tauplate2ftml=-(1+pc_KSFB_increase/100)*alphatauml*rcomu(1)*localswitchFBhelper;
                                end;

tauplate2ft=[tauplate2ftap*1 tauplate2ftml 0];% ap is x in old coord sys
percentnoiseintau=[ftaumusnoise(1)*randn ftaumusnoise(2)*randn randn];
fconstbias=ftaumusnoise(3);constbiasnoise=fconstbias.*[randn randn randn];

taumusnoise=(percentnoiseintau./100).*tauplate2ft;


taumusub2=tauplate2ft+(taumusnoise+constbiasnoise)*noiseonoffswitch;



   

return
end














function  taumusub2=muscle2(rcomu,v ,flagactivelegx,alphatauap,ftaumusnoise,noiseonoffswitch,KsEdgeEnhancingfactor);
  global ktm mset negmset mrandix velnow4beta velpast4beta mold mnew cold cnew ktmgt4 swcori disply         centerborder pauseordisplay2 generaltempstorageshelf betaloadedengleg disenglegload c_lt_begin c_lt_end c_rt_begin c_rt_end cend flagcom flagactiveleg alphatauml slopebeta gravity alphatauRL vtannedswitch

UPFRONTSIGN_AFTER_TAUENGCODE_UNDSTOOD=[1 -1];%[LEFT RIGHT] HERE

localswitchFBmover=1;
localswitchRLhelper=1;
sign_tauplate2flap_lt_rt=[1 -1];
                                if flagactivelegx==1 %left +y
                                    tauplate2ftml=-alphatauml*KsEdgeEnhancingfactor*rcomu(1)*localswitchFBmover;% 1  proj = ap here (diff from expt 2)
                                     tauplate2ftap=UPFRONTSIGN_AFTER_TAUENGCODE_UNDSTOOD(1)*alphatauap*abs(rcomu(2))*localswitchRLhelper;%put on 25May2104

                                elseif flagactivelegx==-1
                                      tauplate2ftap=UPFRONTSIGN_AFTER_TAUENGCODE_UNDSTOOD(2)*alphatauap*abs(rcomu(2))*localswitchRLhelper;%put on 25May2104
                                elseif flagactivelegx==0
                                tauplate2ftml=-2*alphatauml*rcomu(1);

                                tauplate2ftap=0*alphatauap*rcomu(1);
                                end
                                
          
tauplate2ft=[tauplate2ftap*1 tauplate2ftml 0];% ap is x in old coord sys
percentnoiseintau=[ftaumusnoise(1)*randn ftaumusnoise(2)*randn randn];
fconstbias=ftaumusnoise(3);constbiasnoise=fconstbias.*[randn randn randn];
taumusnoise=(percentnoiseintau./100).*tauplate2ft;%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


taumusub2=tauplate2ft+(taumusnoise+constbiasnoise)*noiseonoffswitch;%NOISE SOURCE %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



return
end








