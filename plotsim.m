

%File is a tab-delimited table of simulation results
%Row 1 contains headers
%Column 1 contains time values
%Generate an example in LTspice by ru55555555nning a simulation
% and choosing File -> Export data as text

close all;
%Column 1 contains time values
%Generate an example in LTspice by running a simulation
% and choosing File -> Export data as text

%Input file name
datafile = "simdata.txt";

%Get simulation data.
simfile = importdata(datafile,'\t',1);

simdata = simfile.data;

%Loop over columns of data
figure;
names={'Magnitude (dB)', 'Phase (deg)'};
subplot(size(simdata,2),1,2);
semilogx(simdata(:,1),20*log10(simdata(:,2))); %Add the data
ylabel(names{1});
xlabel('Frequency (Hz)');
yticks([-100 -80 -60 -40 -20 0]);

hold on
subplot(size(simdata,2),1,3);
semilogx(simdata(:,1),simdata(:,3)); %Add the data
ylabel(names{2});
xlabel('Frequency (Hz)');
yticks([-90 -60 -30 0 30 60 90]);
