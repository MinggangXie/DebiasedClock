clear;close all
FontSize=15;
lambda = 1:5:200; % Poisson parameter
lambda = round(logspace(0,3.5));
index_plot=round(linspace(1,length(lambda)*0.9,5));
% Define range of values
dx=1;

% Parameters
% mu = 0; % Gaussian mean
% sigma = 0.2; % Gaussian standard deviation

b=3;
sigma_SR=[0.02 0.08 0.2];
systematic_error_probability_density_function={'Gaussian distribution'};%'Uniform distribution','Gaussian distribution'
% Calculate the probability density function (pdf) of the sum

colors='rgbmck';
set(gcf,'position',[566   93  635  869]);
axes1=axes('position',[0.1    0.74    0.87    0.25]);
cd '../../'
comparison_accuracy_and_Poisson_error(FontSize);
text(-0.585,0.45e7,'(a)','FontSize',round(FontSize*1.1),'FontName','Times New Roman')
cd 'Xie2023EPSL_code\GUI'

axes2=axes('position',[0.1    0.41    0.87    0.25]);
axes3=axes('position',[0.1    0.08    0.87    0.25]);
i_str=0;
for i_se=1:length(systematic_error_probability_density_function)
    for i=1:length(sigma_SR)
        for j=1:length(lambda)
            [x,pdf]=pdf_of_sum(lambda(j),sigma_SR(i),b,systematic_error_probability_density_function{i_se});
            
            pdf_poiss=poisspdf(x, lambda(j));
            pdf_poiss = pdf_poiss / sum(pdf_poiss);
            c=1/max(pdf_poiss);
            cdf_poiss=cumsum(pdf_poiss*dx);
            % Normalize the cdf
            cdf_poiss = cdf_poiss/max(cdf_poiss);
            cdf=cumsum(pdf*dx);
            % Normalize the cdf
            cdf = cdf / max(cdf);
            
            
            
            
            if sum(j==index_plot)>0%j==length(lambda)
%                 if i==1
               axes(axes2)
%                 end
                % Plot the probability density distribution
                %         figure(1)
                if i_se==1&&i==1
                    semilogx(x, pdf_poiss*c,'k-.','LineWidth',3.5);hold on
                    plot([1,1]*lambda(j),[0 1.2],'--','LineWidth',1,'color',[0.5 0.5 0.5])
                end
                %         plot(x, cdf_poiss,'k--')
                if i_se==1
                    plot(x, pdf*c,[colors(i) '-'],'LineWidth',1.5)
                else
                    plot(x, pdf*c,[colors(i+length(sigma_SR)) '-'],'LineWidth',1.5)
                end
                %         plot(x, cdf,'g--')
                
                %         title('Probability Density Distribution of observed crater number')
                drawnow
                ylim([0 1.2])
                %         hold off
            end
            confidence_interval=0.6827;
            %         j
            %         pdf
            [lower_bound_1sigma(j),upper_bound_1sigma(j)]=get_uncertainty(confidence_interval,x,pdf);
            [lower_bound_poiss_1sigma(j),upper_bound_poiss_1sigma(j)]=get_uncertainty(confidence_interval,x,pdf_poiss);
        end
        xlim([0 3000])
        % figure;
        % errorbar(lambda,lambda,lambda-lower_bound_1sigma,upper_bound_1sigma-lambda,'ro');hold on
        % errorbar(lambda,lambda,lambda-lower_bound_poiss_1sigma,upper_bound_poiss_1sigma-lambda,'ko');
        
        axes(axes3)
        sigma_neg=(lambda-lower_bound_1sigma)./lambda;
        sigma_pos=(upper_bound_1sigma-lambda)./lambda;
        if i_se==1
            loglog(lambda,sigma_neg,[colors(i) 'o'],lambda,sigma_pos,[colors(i) 'd'],'markersize',FontSize*0.5);hold on
        else
            loglog(lambda,sigma_neg,[colors(i+length(sigma_SR)) 'o'],lambda,sigma_pos,[colors(i+length(sigma_SR)) 'd'],'markersize',FontSize*0.5);hold on
        end
        % loglog(lambda,lambda-lower_bound_poiss_1sigma,'ro',lambda,upper_bound_poiss_1sigma-lambda,'rd');hold on
        if i_se==1&&i==1
            loglog(lambda,lambda.^0.5./lambda,'k-.','LineWidth',1.5)
        end
        % Set up fittype and options.
        ft = fittype( 'power2' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        % opts.StartPoint = [0.332941131511425 1.00418441372602];
        
        if 0
        % Fit model to data.
        [fitresult, gof] = fit( lambda', sigma_pos', ft, opts );
        if i_se==1
            loglog(lambda,fitresult(lambda),[colors(i) '-'],'LineWidth',1.5)
        else
            loglog(lambda,fitresult(lambda),[colors(i+2) '-'],'LineWidth',1.5)
        end
        
        ylim([min([sigma_neg,sigma_pos]) max([sigma_neg,sigma_pos])])
        end
        xlim([0 3000])
        ylim([1e-2 2])
        i_str=i_str+1;
        str_txt{i_str}=sprintf('%s, \\sigma_{S}=%.2f',systematic_error_probability_density_function{i_se},sigma_SR(i));
    end
end
axes(axes2)
set(axes2,'FontSize',round(FontSize*0.65),'FontName','Times New Roman','TickLength',[0.015 0.012],'linewidth',0.8,...
    'xminortick','on','yminortick','on','box','on');hold on
xlabel('Observed number of craters,\itn_{\rmobs}','FontSize',FontSize,'FontName','Times New Roman')
ylabel('Probability Density','FontSize',FontSize,'FontName','Times New Roman')
text(0.5,1.2,'(b)','FontSize',round(FontSize*1.1),'FontName','Times New Roman')

axes(axes3)
set(axes3,'FontSize',round(FontSize*0.65),'FontName','Times New Roman','TickLength',[0.015 0.012],'linewidth',0.8,...
    'xminortick','on','yminortick','on','box','on');hold on
xlabel('Expected actual number of craters','FontSize',FontSize,'FontName','Times New Roman')
ylabel('1\sigma relative uncertainty','FontSize',FontSize,'FontName','Times New Roman')
text(0.5,2.2,'(c)','FontSize',round(FontSize*1.1),'FontName','Times New Roman')

for i=1:length(str_txt)
    if length(str_txt)==3
        text(1.3,0.05*1.6^-(i-1),str_txt{i},'FontSize',round(FontSize*0.75),'FontName','Times New Roman','color',colors(i));
    elseif length(str_txt)==6
        text(1.3,0.16*1.6^-(i-1),str_txt{i},'FontSize',round(FontSize*0.75),'FontName','Times New Roman','color',colors(i));
    end
end   


print('-djpeg90','-r400',['../../Figure/probability_density_function_of_systematic_error_in_number_b' num2str(b) '.jpg'])







