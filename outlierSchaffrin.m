function [Mout]=outlierSchaffrin(tobs,yobs,Pd_obs,algorithm)
%   given time series, output the outliers.
% Using Dr. Burkhard Schaffrin's algorithm for detection of a single outlier.
% given: tobs: time epoch vector in year, n by 1
%        yobs: elevation vector, n by 1
%        Pd_obs: the diagnal components of the weight matrix, n by 1, assume there are no off-diagnal components.
%        if Pd_obs=[], use equal weight;
%        algorithm: 1, or 2. If algorith==1, the fitting model is y=a+bt; if 2, model is y=a; 

%Output Mout: logical vector, 0 not outlier, 1 outlier
%  Usage: [Mout]=outlierSchaffrin(tobs,yobs,[],2);
%       similar to Mout=isoutlier(yobs);
%       Code written by Chunli Dai, March 26 2021.

    Mout=false(size(yobs));
    ni=length(yobs);
    AMa=ones(ni,1);
    yr=365.25;
    t=tobs*yr;
    yobs_org=yobs;tobs_org=tobs;
    flagplot=1;
    
    if isempty(Pd_obs) %assume equal weight
        Pd_obs=ones(size(yobs));
    end
    
%     P=diag(Pd_obs);
%     wd=zeros(size(yobs)); wd(:)=Pd_obs.^0.5; %use wd to use the faster matlab algorithm \.

    tm=mean(tobs);
    if (algorithm == 1 ) %linear trend
        mp=2; 
        model='y=a+bt';
    elseif (algorithm == 2) % constant
        mp=1;
        model='y=a';
    end

    % original GMM that does not include an outlier.
    idout=[];idoutpre=[];
    idkp=1:length(yobs);
    
    count=0;flagcond=0;
    
    while 1 %iteration until no more outlier is detected.
    idkp= idkp(~ismember(idkp,idout));
    P=diag(Pd_obs(idkp)); 
    yobs=yobs_org(idkp);tobs=tobs_org(idkp);
    wd=zeros(size(yobs)); wd(:)=Pd_obs(idkp).^0.5;
    
    ny=length(yobs);
    nr=ny-mp;%%number of redundancy.
    if nr<=1
        Mout=false(size(yobs_org)); return;    %unable to detect outliers
    end
    
    %search for each observation to see if it is an outlier
    id=[];
    for j=1:ny
    AM=[];
    AM=zeros(length(yobs),mp);  %
    AM(:,1)=AMa(idkp,:);
    if (algorithm == 1 ) %linear trend            
    AM(:,2)=tobs-tm;
    elseif (algorithm == 2) % constant
        %do nothing
    end

    cdA=cond(AM'*P*AM);
    if(cdA>1e5) 
        %display(['condtion number of AM*P*AM is large: ',num2str(cdA),';isel:',num2str(isel)]);
        display(['condtion number of AM*P*AM is large: ',num2str(cdA)]);
        flagcond=1;
%         break
    end
%     est=inv(AM'*P*AM)*AM'*P*yobs;
    est=(wd.*AM)\(wd.*yobs);
    %             fit=est(1) + est(2)*(epoch(idkp)-tm);
    fit=AM*est;
    etilde=yobs-AM*est;
    
    %see /Users/chunlidai/Downloads/651Adjustment computations2/Final08_2_outlier.m
    Omiga=etilde'*P*etilde;
    Qe=inv(P)-AM*inv(AM'*P*AM)*AM';
    Rj=etilde(j)^2/Qe(j,j); %equation 10.106 in adjustments
    Tj=Rj*(nr-1)/(Omiga-Rj);

    %F-distribution https://www.mathworks.com/help/stats/finv.html
    v1=1; %number of parameters
    v2=ny-mp-1; %
%     alpha=0.05; % significance level of 5%.
    alpha=1-0.9973; % significance level.
    xf=finv(1-alpha,v1,v2); %%You would observe values greater than xf only 0.27% of the time by chance.
    %Hypothesis:Test H0: j is not an outlier; against Ha: j is an outlier.
    if Tj<=xf %accept H0
    else % reject H0
        %j is an outlier
        id=[id(:);j];
    end
    
    end
    idout=[idoutpre;idkp(id)];idsign=id;
    count=count+1;
    
    if flagplot==1 && 1% plotting
    tkp=t(idkp);[ts,idsort]=sort(tkp);

    figure % (1)
    set(gcf,'Color','white')
    set(gca,'FontSize', 18);
    set(gcf, 'PaperPosition', [0.25 2.5 4 3]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
    hold all
    plot(t(idkp(idsort)),yobs_org(idkp(idsort)),'b>-','MarkerSize',12,'linewidth',4)
    Mt1=Pd_obs<1; %less weight;
    plot(t(Mt1),yobs_org(Mt1),'m+','MarkerSize',18,'linewidth',4)
    plot(t(idout),yobs_org(idout),'ks','MarkerSize',12,'linewidth',4)
    plot(t(idkp(idsort)),fit(idsort),'g*-','MarkerSize',12,'linewidth',4)
    datetick('x','mm/yy')
    box on
    ylabel('DEM time series') 
    title(['Outliers detection model: ',model,', iter=',num2str(count)])
    end % plotting
    
    
    if isempty(idsign); break;   end %If no outliers detected,  stop the iteration
    if count>=9 %99
    %                 fprintf(['Too many iterations ',count,']);
        break
    end
    idoutpre=idout;    
    end %while 1
    
    Mout(idout)=1;
    
    if 0 %Final plot for all outliers
%       Mout=isoutlier(yobs);
%       t=tobs*yr;   yobs_org=yobs;tobs_org=tobs;
    [ts,idsort]=sort(t);
    fit=nanmedian(yobs_org)*ones(size(yobs_org));
    multi=3;
    stdi=mad(yobs_org,1); %median absolute deviation: median(abs(X ? median(X))).
    figure % (1)
    set(gcf,'Color','white')
    set(gca,'FontSize', 18);
    set(gcf, 'PaperPosition', [0.25 2.5 4 3]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
    hold all
    plot(t((idsort)),yobs_org((idsort)),'b>-','MarkerSize',12,'linewidth',4)
    plot(t(Mout),yobs_org(Mout),'ks','MarkerSize',12,'linewidth',4)
    if 1 
    plot(t((idsort)),fit(idsort),'g*-','MarkerSize',12,'linewidth',4)
    plot(t((idsort)),fit(idsort)-multi*stdi,'r-','linewidth',4)
    plot(t((idsort)),fit(idsort)+multi*stdi,'r-','linewidth',4)   
    end
    datetick('x','mm/yy')
    box on
    ylabel('DEM time series') 
    title(['Outliers detection'])
    end

return
end