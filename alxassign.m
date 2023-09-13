Filepath1=input('Enter the file path to the excel workbook where your wind data (direction and speed data) is stored.\n');
Sheetnumber1=input('Enter the sheet number conntaining the hourly windspeed data, make sure only wind speed data is stored here\n');
cell11=input('Enter the cell number of the first wind speed data entry in exxcell sheet e.g A1\n');
cell12=input('Enter the cell number of the last wind speed data enetry in excell sheet e.g C32\n');
Sheetnumber2=input('enter the sheet number conntaining the corresponding hourly wind direction data\n');
cell21=input('Enter the cell number of the first wind direction data entry in exxcell sheet e.g A1\n');
cell22=input('Enter the cell number of the last wind direction data enetry in excell sheet e.g C32\n');
Sheetnumber3=input('enter the sheet number conntaining the Time(month and year) vector whose elements are the month and year in cosecutive columns say A(month) and B(year) corresponding to each wind speed entry\n');
cell31=input('Enterthe cell number of the first entry in months column vector data entry in excell sheet e.g A1\n');
cell32=input('Enter the cell number of the last entry in months column vector in excell sheet e.g B32\n');
c=input('enter number of months in a year\n');
b=input('enter threshold percentage\n');
density=input('input air density\n');

%filename=input('Enter the File Name:Make sure this file has been created in your desktop\n');
Range1=sprintf('%s:%s', cell11, cell12);
Range2=sprintf('%s:%s', cell21, cell22);
Range3 = sprintf('%s:%s', cell31, cell32);
WSD0 = xlsread(Filepath1,Sheetnumber1,Range1)
DirectionData = xlsread(Filepath1,Sheetnumber2,Range2)
Time_data = xlsread(Filepath1,Sheetnumber3,Range3);
Month_data=Time_data(:,1);
YearsData=Time_data(:,2);
%pointer matrices, that allows for the conservation of the arrangement of
%the data in the matrices during dletion, insertion and so on
a=25;
Days_pointers= cumsum(ones(size(WSD0,1),24),1);
Hours_pointers= repmat(1:24,size(WSD0,1),1);

%making column vectors from the pointer matrices and wind speed matrix that were made
A= reshape(WSD0.',[],1);
B=Month_data(1:numel(A));
Days_pointersf= reshape(Days_pointers.',[],1);
Hours_pointersf= reshape(Hours_pointers.',[],1);

 %start of data filtration process
 %C1=[A B];
% First phase of filtration (Checking for the percentage contribtion of the invalid cells)
 S = numel(A(A<0))+numel(A(A>25));
 a1=a*24;
 N =  numel(A);
 percentage_S = (S/N)*100;
 if percentage_S<b
     fprintf('Data is good proceed with analysis %g ', percentage_S)
 else 
     return
 end

 %second phase of filtration
nrows= numel(A);
ncols=1;
for i=1:nrows
    for j=1:ncols
        if A(i,j)<0
            A(i,j)=0.1;
        elseif A(i,j)>25
            A(i,j)=25;
        elseif A(i,j)==0
            A(i,j)=0.5;
        else
        end
    end 
end

A = [A;0];
%First NaN
idx = find(isnan(A),1);
%First number after NaN
k = find(~isnan(A(idx+1:end)),1)+idx;
while ~isempty(k)
    if k-idx>=6
        A(idx:k-1)=[];
        Days_pointersf(idx:k-1)=[]; %WDD must be the same dimensions as TDD 
        Hours_pointersf(idx:k-1)=[];
        YearsData(idx:k-1)=[];
        B(idx:k-1)=[];
        k=idx;
    end
    idx = find(isnan(A(k+1:end)),1)+k;
    k = find(~isnan(A(idx+1:end)),1)+idx;
end
%Remove the number
A = A(1:end-1);
A = fillmissing(A,'movmean',10); %fill the missing data
%the third stage of data filtration
for i=1:1:size(WSD0,1)
    selectedRowsLogicalIndices = Days_pointersf == i;
    idxcount=numel(Days_pointersf(Days_pointersf==i));
    if idxcount<=12
       A(selectedRowsLogicalIndices)=[];
       Days_pointersf(selectedRowsLogicalIndices)=[];
       Hours_pointersf(selectedRowsLogicalIndices)=[];
       YearsData(selectedRowsLogicalIndices)=[];
       B(selectedRowsLogicalIndices)=[];
    end
    if i==size(WSD0,1)
        A;
    end
end
%Deal with the zeros
nrows= numel(A);
for i=1:nrows
    for j=1:ncols
        if A(i,j)==0.1
            %A(i,j)=mean(A)*rand(1);
            A(i,j)=NaN;
        
        end
    end 
end
A = fillmissing(A,'movmean',10);
Clean_data=[A  Hours_pointersf Days_pointersf B];
%Hourly mean and stdv category
%D=A(116:end);
percentage_S %output remeber to call it at the end
for k=1:1:24
    selectedRowsLogicalIndices = Hours_pointersf == k;
    subsetAdata = A(selectedRowsLogicalIndices);
    Hourly_mean(k,1)=mean(subsetAdata);
    Hourly_std(k,1)=std(subsetAdata);
    if k==24
        display(Hourly_mean)
        display(Hourly_std)
        Hours_of_day=[1:1:24];
        %turbulence intensity evolution over the day
        Turbulence_intensity_daily=(Hourly_std)./(Hourly_mean)
        %Codes for exporting the results to excell worksheets
        Hrs_of_the_day=["00:00","01:00","02:00","03:00","04:00","05:00","06:00","07:00","08:00","09:00","10:00","11:00","12:00","13:00","14:00","15:00","16:00","17:00","18:00","19:00","20:00","21:00","22:00","23:00"]'
        T1=table(Hrs_of_the_day,Hourly_mean);
        T2=table(Hrs_of_the_day,Hourly_std);
        T3=table(Hrs_of_the_day,Turbulence_intensity_daily);
        filename = 'Wind Analysis Ouput.xlsx';
        writetable(T1,filename,'Sheet',1,'Range','A1')
        writetable(T2,filename,'Sheet',2,'Range','A1')
        writetable(T3,filename,'Sheet',3,'Range','A1')
        %plotting mean wind speed varriation during the day
        figure1=figure();
        plot(Hours_of_day,Hourly_mean)
        xlabel('Hours of the day hrs');
        ylabel('Wind speed m/s');
        set(gca,'XTick',Hours_of_day)
        set(gca,'XTickLabel',str2mat('00:00','01:00','02:00','03:00','04:00','05:00','06:00','07:00','08:00','09:00','10:00','11:00','12:00','13:00','14:00','15:00','16:00','17:00','18:00','19:00','20:00','21:00','22:00','23:00'))
        grid on
        grid minor
        %turbulence intensity plot
        figure3=figure();
        plot(Hours_of_day,Turbulence_intensity_daily)
        xlabel('Hours of the day hrs');
        ylabel('Turbulence intensity');
        set(gca,'XTick',Hours_of_day)
        set(gca,'XTickLabel',str2mat('00:00','01:00','02:00','03:00','04:00','05:00','06:00','07:00','08:00','09:00','10:00','11:00','12:00','13:00','14:00','15:00','16:00','17:00','18:00','19:00','20:00','21:00','22:00','23:00'))
        grid on
        grid minor
    end
end

%Annual Weubull parameters are obtained via wblfit() function. The
%finctinin executes Maximum likelihood method and it is actually inbiult in
%matlab
[Overal_parameters]= wblfit(A);
figure6=figure();
histogram(A,'Normalization','probability')
hold on
plot(sort(A),wblpdf(sort(A),Overal_parameters(1),Overal_parameters(2)))
legend('Observed Samples','Estimated Distribution')
xlabel('Mean wind speed (m/s)')
ylabel('Annual probabilty density')
hold off

%plot cumulative distribution function
figure7=figure();
plot(sort(A),wblcdf(sort(A),Overal_parameters(1),Overal_parameters(2)))
hold on
cdfplot(A)
legend('Estimated Distribution','Observed Samples')
grid off
xlabel('Mean wind speed m/s')
ylabel('Probability')
title('')
%Error analysis using the RMSE technique and the R-squared techniques.
%Error approximation by making use of cummulative distribution
%determine the emprical cumulative distribution
Empc=ecdf(A);
%weibull cumulative
Estc=wblcdf(A,Overal_parameters(1),Overal_parameters(2));
%group the Estc
figure
Hs=histogram(Estc,numel(Empc));
binEdges = Hs.BinEdges;
x = binEdges(1:end-1) + Hs.BinWidth/2;
fitlm(Empc,x')

% This following division is the loop that regrouped data into their
% respective month of the year. Monthly; mean wind speed, standard
% deviation, weibull parameters, turbulence intensity, power density, most
% probable wind speed and wind speed carrying maximum energy are determined
% in the loop. You should note that, corresponding pictorial reperesentation
% of data is determined result from the loop
for i=1:1:c
     selectedRowsLogicalIndices = B == i;
     if numel(B(B==i))>a1&numel(B(B==i))>a1;
         subsetAdata = A(selectedRowsLogicalIndices);
         Monthly_mean(i,1)=mean(subsetAdata);
         Monthly_standard_deviation(i,1)=std(subsetAdata);
         Monthly_parameters(:,i)= wblfit(subsetAdata);
         if i==c
             display(Monthly_mean)
             display(Monthly_standard_deviation)
             display(Monthly_parameters)
             display(Overal_parameters)
            
             %turbulence intensity evolution over the year
             Turbulence_intensity_overtheyear=(Monthly_standard_deviation)./(Monthly_mean)
             
             % Exporting output data to Excell
             Months_of_the_yr=["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sept","Oct","Nov","Dec"]'
             T4=table(Months_of_the_yr,Monthly_mean);
             T5=table(Months_of_the_yr,Monthly_standard_deviation);
             T6=table(Months_of_the_yr,Turbulence_intensity_overtheyear);
             T61=table(Months_of_the_yr,(Monthly_parameters)');
             filename = 'Wind Analysis Ouput.xlsx';
             writetable(T4,filename,'Sheet',4,'Range','A1')
             writetable(T5,filename,'Sheet',5,'Range','A1')
             writetable(T6,filename,'Sheet',6,'Range','A1')
             writetable(T61,filename,'Sheet',6,'Range','A15')
             
             %mean wind speed over the months pf the year plot
             month_of_theyear=[1:1:12];
             figure2=figure();
             bar(month_of_theyear,Monthly_mean)
             set(gca,'XTick',month_of_theyear)
             set(gca,'XTickLabel',str2mat('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sept','Oct','Nov','Dec'))
             xlabel('Months of the year');
             ylabel('Mean wind speed m/s');
             figure4=figure();
             bar(month_of_theyear,Turbulence_intensity_overtheyear)
             set(gca,'XTick',month_of_theyear)
             set(gca,'XTickLabel',str2mat('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sept','Oct','Nov','Dec'))
             xlabel('Months of the year');
             ylabel('Turbulence intensity');
             %the plot of wind speed varriation with height from the
             %exrapolated wind speed, using decadal mean as the constant
             %speed
             H11=10:1:100;
             %z_0=exp((mean(A)*log(z_y)-v_y*log(10))/(mean(A)-v_y));
             alph_a=0.3%(0.37-0.088*log(mean(A)))/(1-0.088*log(z_0/10));
             for r=1:1:numel(H11)
                 Vz(r,:)=mean(A)*(H11(r)/10)^(alph_a);
                 if r==numel(H11)
                     display(Vz);
                     figure5=figure();
                     plot(Vz,H11)
                     ylabel('Height (m)')
                     xlabel('Wind speed (m/s)')
                 end
             end
             
             for n=1:c
                 Power_density_monthly(n,:)=0.5*density*((Monthly_parameters(1,n)).^3)*gamma((Monthly_parameters(2,n)+3)/Monthly_parameters(2,n));
                
                 if n==c
                     display(Power_density_monthly)
                     Power_density_overall=0.5*density*((Overal_parameters(1)).^3)*gamma((Overal_parameters(2)+3)/Overal_parameters(2));
                     figure8=figure();
                     bar(month_of_theyear,Power_density_monthly)
                     set(gca,'XTick',month_of_theyear)
                     set(gca,'XTickLabel',str2mat('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sept','Oct','Nov','Dec'))
                     xlabel('Months of the year');
                     ylabel('Wind power density (W/m^2)');
                     Power_density=[Power_density_monthly;Power_density_overall];
                     Months_of_the_yr=["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sept","Oct","Nov","Dec","Annual wind power density"]'
                     T9=table(Months_of_the_yr,Power_density);
                     filename = 'Wind Analysis Ouput.xlsx';
                     writetable(T9,filename,'Sheet',9,'Range','A1')
                     
                 end
             end         
         end
     else
        F(i,:)=find(A ==i);
        if i==c
            display(F)
        end

     end 
end

%Power dnsity extrapolation with respect to height. The equation executed
%was derived fromthe formula for wind speed extrapolation (refer to the report)
H1=10:1:100; %matrix of heights
%Extrapolations of Weibull parameters
for r=1:1:numel(H11)
    for i=1:numel(Monthly_parameters(1,:))
    c_z(r,i)=Monthly_parameters(1,i)*(H11(r)/10)^(alph_a);
    k_z(r,i)=Monthly_parameters(2,i)/(1-0.088*log(H11(r)/10));
    c_zoverall(r,1)=Overal_parameters(1)*(H11(r)/10)^(alph_a);
    k_zoverall(r,1)=Overal_parameters(2)/(1-0.088*log(H11(r)/10));
    if i==12
        display(c_z);
        display(k_z);
        display(c_zoverall);
        display(k_zoverall);
    end
    end
end
%power extrapolation
for j=1:1:12
    Power_density_monthlyextr(:,j)=0.5*density*((c_z(:,j)).^3).*gamma((k_z(:,j)+3)./k_z(:,j));
    if j==12
       Power_density_overallextr(:,1)=0.5*density*((c_zoverall(:,1)).^3).*gamma((k_zoverall(:,1)+3)./k_zoverall(:,1));
       display( Power_density_monthlyextr);
       figure();
       plot(Power_density_monthlyextr,H1)
       ylabel('Height (m)')
       xlabel('Power density W/m^2')
       legend ('Jan','Feb','Mar','Apr','May','June','Jul','Aug','Sept','Oct','Nov','Dec');
       figure9=figure();
       plot(Power_density_overallextr,H1)
       ylabel('Height (m)')
       xlabel('Power density W/m^2')
       end
end


%WSD0=input('Enter wind wind data(Multidimentional matrix)\n');
%TE=input('Enter corresponding time data for the time span data was taken in months\n');%time_entry
thresh=6;
%First we convert the data which was in multidimensional matrix to cmoulunm
%vector
%colunm vector is easier to work with as opposed to multidimentional matrix
TDD = reshape(DirectionData.',[],1);% Transposing the matrix "DirectionData" row-wise into a Column Matrix
WSD1= reshape(WSD0.',[],1);% Transposing the matrix "wind speed data (WSD)" row-wise into a Column vector
%converting the wind speed matrix to have the same dimensions as the TDD
WSD=WSD1(1:numel(TDD));
%TDD is the trasposed direction matrix which is basically 1D vector
% The fisrt step of data filtering, checking for the negatives, NaN and
% exegerated values within the matrix.
nrows= numel(TDD);
ncols=1;
for i=1:nrows
    for j=1:ncols
        if TDD(i,j)<0 %checking for negative entries
            TDD(i,j)=0.005; %fixing the negative entries to a number close to zero
        elseif TDD(i,j)>360
            TDD(i,j)=360;
        else
        end
    end 
end
for i=1:nrows
    for j=1:ncols
        if WSD(i,j)<0
            WSD(i,j)=0.005;
        elseif WSD(i,j)>25
            WSD(i,j)=25;
        else
        end
    end  
end
%Adding a number to the end to account fro the perriferal data
TDD = [TDD;0];
WSD=[WSD;0];
%First NaN
idx = find(isnan(TDD),1);
%First number after NaN
k = find(~isnan(TDD(idx+1:end)),1)+idx;
while ~isempty(k)
    if k-idx>=thresh
        TDD(idx:k-1)=[];
        WSD(idx:k-1)=[]; %WDD must be the same dimensions as TDD 
        k=idx;
    end
    idx = find(isnan(TDD(k+1:end)),1)+k;
    k = find(~isnan(TDD(idx+1:end)),1)+idx;
end
%Remove the number
TDD = TDD(1:end-1);
disp(TDD);


idx0 = find(isnan(WSD),1);
j = find(~isnan(WSD(idx0+1:end)),1)+idx0;
while ~isempty(j)
    if j-idx0>=thresh
        WSD(idx0:j-1)=[];
        TDD(idx0:j-1)=[];
        %WDD must be the same dimensions as TDD 
        j=idx0;
    end
    idx0 = find(isnan(WSD(j+1:end)),1)+j;
    j = find(~isnan(WSD(idx0+1:end)),1)+idx0;
end
%Remove the number
WSD = WSD(1:end-1);
disp(WSD);
% fill in the NaNs that do not meet the the deletion criteria using moving
% average technique using the boaering ten numbers
F = fillmissing(TDD,'movmedian',10);
F1 = fillmissing(WSD,'movmedian',10);
FF=[F1 F];
%triming the time data to match the dimension of the F matrix, this is
%important for monthwise data exyraction
%TE1=TE(1:numel(F));
%plot the wind rose chart by calling the wind_rose() function which has
%been saved as a function file.

%The windrose of the whole time sampling period
figure10=figure
      wind_rose(F,F1)%calling the wind_rose() function
disp('Press any key to continue, you view your graphs before pressing any key')  % Press a key here.You can see the message 'Paused: Press any key' in        % the lower left corner of MATLAB window.
pause;
Reportte %report generation
