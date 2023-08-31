import mlreportgen.report.* 
import mlreportgen.dom.* 
rpt = Report('Wind Energy Assessment','html'); 

tp = TitlePage; 
tp.Title = 'WIND POWER POENTIAL REPORT'; 
tp.Author = 'Steven and Alvin Maingi'; 
add(rpt,tp); 
add(rpt,TableOfContents); 

ch1 = Chapter; 
ch1.Title = 'INTRODUCTION'; 
sec1 = Section; 
sec1.Title = 'Background'; 
para = Paragraph(['Wind assessment is a crucial process in evaluating the feasibility and potential of harnessing wind energy resources. It involves a comprehensive analysis of various factors, including wind direction, temporal and spatial wind speed variation, turbulence intensity, and wind power density. By examining the prevailing wind direction, which indicates the predominant path from which the wind blows, experts can determine the most suitable locations for wind farms and optimize turbine placements. Temporal and spatial wind speed variations, as well as turbulence intensity, are also assessed to understand the consistency, intensity, and stability of wind resources over time and across different geographical areas. Furthermore, wind power density is a key metric that quantifies the energy available in the wind and helps estimate the potential electricity generation capacity. In this assessment, wind direction, spatial and temporal wind speed variation, turbulence intensity, and wind power density will be evaluated to provide comprehensive insights for efficient wind energy utilization.']); 
add(sec1,para) 
add(ch1,sec1)
add(rpt,ch1)

ch2 = Chapter(); 
ch2.Title = sprintf('THE ASSESSMENT RESULTS'); 
sec1 = Section; 
sec1.Title = 'Background';
para=Paragraph(['In this section, a short report on every varriable mentioned previously will be presented. It should be noted that the  iterpretation of the results is left to the user. The application only gives an organized preview and visualization of the data but does not interpret the results'])
add(sec1,para) 
add(ch2,sec1)
sec2 = Section; 
sec2.Title = 'Temporal wind speed variation';
para=Paragraph(['Figure 1 dsiplays the wind speed behaviour of the area over the hours of the day. From the figure, we can deduce the number of hours that wind energy generation is suitable. The hours with wind speed graeter than 3 m/s whch is the cut-in wind speed. On the oher hand Figure 2 diplsys the wind speed eveolution over the month of the year. From the bar graph, wwe can duduced the behaviour of the wind regime over the entire year. '])
add(sec2,para) 
add(ch2,sec2)
fig = Figure(figure1);
fig.Snapshot.Height = '4in'; 
fig.Snapshot.Width = '6in';
fig.Snapshot.Caption = sprintf('Diurnal wind speed variation'); 
add(ch2,fig);
fig2 = Figure(figure2);
fig2.Snapshot.Height = '4in'; 
fig2.Snapshot.Width = '6in';
fig2.Snapshot.Caption = sprintf('Wind speed variation over the months of the year');
add(ch2,fig2);

sec3 = Section; 
sec3.Title = 'Turbulence Intensity';
para=Paragraph(['Figure 3 illustrates the turbulence intensity (TI) evolution over the hours of the day. If the turbulence intensity exceeds 0.25, then the win regime is considered as turbulent and wind turbines installed in the area at that height must take into account the turbulence intensity. Figure 4 likewise illustrates the turbulence intensity over the year. Judgment criteria applied in diurnal TI is applicable as well in monthly turbulence intensities.'])
add(sec3,para)
add(ch2,sec3)
fig3 = Figure(figure3);
fig4 = Figure(figure4);
fig3.Snapshot.Height = '4in'; 
fig3.Snapshot.Width = '6in';
fig4.Snapshot.Height = '4in'; 
fig4.Snapshot.Width = '6in';
fig3.Snapshot.Caption = sprintf('Diurnal wind turbulence intensity variation'); 
fig4.Snapshot.Caption = sprintf('Wind turbulence intensity variation over the months of the year');
add(ch2,fig3);
add(ch2,fig4);

sec4 = Section; 
sec4.Title = 'Wind speed variation with height';
para=Paragraph(['Usually, wind speed data used in analysis is at 10 m. In most cases, most wind turbines are placed at heights greater than 10 m, hence it is important to evaluate how wind speeds will vary with respect to height. Therefore, figure 5 has been plotted to reveal the expected wind speeds various heights up to 100 m.'])
add(sec4,para)
add(ch2,sec4)
fig5 = Figure(figure5);
fig5.Snapshot.Height = '4in'; 
fig5.Snapshot.Width = '6in';
fig5.Snapshot.Caption = sprintf('Mean wind speed variation with height'); 
add(ch2,fig5);

sec5 = Section; 
sec5.Title = 'Turbulence Intensity';
para=Paragraph(['Figure 6 illustrates the actual (histogram) and estimated (line curve) probability distributions compared in the same axes. The function probability distribution function (pdf) is two parameter Weibull distribution based on the fact that it is the most applied distribution. If the Weibull pdf (which will always be the case) fits the actual data well such that the R2 value is greater than or equal to 0.7, then the graph in figure 6 can be used to determine the probability of occurrence of a particular wind chosen on the x axis. Figure 7 is the cumulative distribution plot that can be used to evaluate the probability of a given wind speed read from the x axis not exceeding a particular speed., especially the cut-in wind speed.'])
add(sec5,para)
add(ch2,sec5)
fig6 = Figure(figure6);
fig7 = Figure(figure7);
fig6.Snapshot.Height = '4in'; 
fig6.Snapshot.Width = '6in';
fig7.Snapshot.Height = '4in'; 
fig7.Snapshot.Width = '6in';
fig6.Snapshot.Caption = sprintf('Wind speed probability density deistribution'); 
fig7.Snapshot.Caption = sprintf('Wind speed cumulative probability density distribution');
add(ch2,fig6);
add(ch2,fig7);

sec6 = Section; 
sec6.Title = 'Wind power density';
para=Paragraph(['Wind power potential is quantified by power density. Mean wind power density has been summarized in table 1. Monthly wind power densities are summarized in figure 8, while figure 9 show the extrapolated wind power density at various heights. Any wind power that exceeds 200 W/m2 is viable for large scale wind power investment.'])
add(sec6,para)
add(ch2,sec6)
fig8 = Figure(figure8);
fig9 = Figure(figure9);
fig8.Snapshot.Height = '4in'; 
fig8.Snapshot.Width = '6in';
fig9.Snapshot.Height = '4in'; 
fig9.Snapshot.Width = '6in';
fig8.Snapshot.Caption = sprintf('Monthly mean wind power density'); 
fig9.Snapshot.Caption = sprintf('Extrapolated mean wind powerv density at various heights');
add(ch2,fig8);
add(ch2,fig9);
tbl = Table(T9); 
tbl.Style = {... 
RowSep('solid','black','1px'),... 
ColSep('solid','black','1px'),}; 
tbl.Border = 'double'; 
tbl.TableEntriesStyle = {HAlign('center')}; 
add(ch2,tbl); 
%add(rpt,ch2); 

sec7 = Section; 
sec7.Title = 'Wind Direction';
para=Paragraph(['Figure 10 shows the predominant wind direction which dictates the orientation of wind turbines.'])
add(sec7,para)
add(ch2,sec7)
fig10 = Figure(figure10);
fig10.Snapshot.Height = '4in'; 
fig10.Snapshot.Width = '6in';
fig10.Snapshot.Caption = sprintf('Wind diretion of the wind regime'); 
add(ch2,fig10);
add(rpt,ch2); 

delete(gcf) 
close(rpt)
rptview(rpt)
