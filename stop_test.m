clear all; clc; close all;
set(gcf,'color','w')                         
IS=0.25;                                
wlen=1100;                              
inc=440;                                
filedir=[];                              
filename='NM_TT_2_01.wav';              
fle=[filedir filename]                      
SNR=5;                                 
[xx,fs]=audioread(fle);                     
y=enframe(x,wnd,inc)';                    
fn=size(y,2);                             
frameTime=frame2time(fn, wlen, inc, fs);
h=waitbar(0,'Running...');                  
for k=1 : fn
    v=y(:,k);                           
    imf=emd(v);                         
    L=size(imf,1);                       
    Etg=zeros(1,wlen);
    for i=1 : L                           
        Etg=Etg+steager(imf(i,:));
    end
    Tg(k,:)=Etg;
    Tgf(k)=mean(Etg);                   
    waitbar(k/fn,h)                     

    set(h,'name',['¶Ëµã¼ì²â - ' sprintf('%2.1f',k/fn*100) '%'])
end
close(h)                                
df=fs/wlen;                             
fx1=fix(250/df)+1; fx2=fix(3500/df)+1;  
km=floor(wlen/8);                       
K=0.5;                                 
for i=1:fn
    A=abs(fft(y(:,i)));                    
    E=zeros(wlen/2+1,1);            
    E(fx1+1:fx2-1)=A(fx1+1:fx2-1);         
    E=E.*E;                            
    P1=E/sum(E);                       
    index=find(P1>=0.9);                 
    if ~isempty(index), E(index)=0; end 
    for m=1:km                         
        Eb(m)=sum(E(4*m-3:4*m));
    end
    prob=(Eb+K)/sum(Eb+K);              
    Hb(i) = -sum(prob.*log(prob+eps));      
end
Tgfm1=multimidfilter(Tgf,10);              
Tgfm=movmedian(Tgfm1,20);
Enm=multimidfilter(Hb,10);                
% Enm=movmedian(Enm1,20);
Mtg=max(Tgfm);                          
Tmth=mean(Tgfm(1:NIS));
T1=0.00005;
T2=0.0002;
Me=min(Enm);                            
eth=mean(Enm(1:NIS));
Det=eth-Me;
T3=0.905*Det+Me;
T4=0.5*Det+Me;                          
[voiceseg,vsl,SF,NF]=vad_param2D_revr(Tgfm,Enm,T1,T2,T3,T4);
%% ×÷Í¼
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1), pos(2)-200,pos(3),pos(4)+150]) 
subplot 311; plot(time,x,'k');
title('Spectrum Proportion');
ylabel('Percentage'); axis([0 max(time) -1.2*max(x) 1.2*max(x)]);

subplot 312; plot(frameTime,Tgfm,'k'); 
title('Marginal Spectrum TEO'); ylim([0 1.2*Mtg]);
ylabel('Amplitude'); axis([0 max(time) 0 1.2*max(Tgfm)])
line([0 max(time)], [T1 T1],'color','k','lineStyle','--');
line([0 max(time)], [T2 T2], 'color','k','lineStyle','-');
subplot 313; plot(frameTime,Enm,'k');
title('Marginal Spectrum MFCC'); ylim([0 1.2*max(Enm)]);
xlabel('Frames'); ylabel('Average MFCC'); axis([0 max(time) 0 1.2*max(Enm)])
line([0 max(time)], [T3 T3], 'color','k','lineStyle','--');
line([0 max(time)], [T4 T4], 'color','k','lineStyle','-');
nx=0;nxm=0;
for k=1 : vsl
    nx1=voiceseg(k).begin; nx2=voiceseg(k).end;
    nx3=((nx2-nx1-1)*inc+wlen)*1/44100;                
    nx=nx+nx3;                                     
    if nx3>nxm                                      
        nxm=nx3;
    end    
    fprintf('%4d   %4d   %4d   %4d\n',k,nx1,nx2,nx3);
    figure(1); subplot 311; 
    line([frameTime(nx1) frameTime(nx1)],[-1.2*max(x) 1.2*max(x)],'color','k','lineStyle','-');
    line([frameTime(nx2) frameTime(nx2)],[-1.2*max(x) 1.2*max(x)],'color','k','lineStyle','--');
    subplot 312; 
    line([frameTime(nx1) frameTime(nx1)],[0 1.2*Mtg],'color','k','lineStyle','-');
    line([frameTime(nx2) frameTime(nx2)],[0 1.2*Mtg],'color','k','lineStyle','--');
    subplot 313;
    line([frameTime(nx1) frameTime(nx1)],[0 1.2*max(Enm)],'color','k','lineStyle','-');
    line([frameTime(nx2) frameTime(nx2)],[0 1.2*max(Enm)],'color','k','lineStyle','--');
end

t=k-1;                 
spro=nx/max(time);     
maxn=nxm;           