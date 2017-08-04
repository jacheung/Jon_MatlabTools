depth=T.depth;
dp=num2str(depth);
figure;T.plot_spike_raster(0,'BehavTrialNum');
print(gcf,'-dpng',['Z:\Users\Jon\Presentations\LabMeetings\160805\' dp '_' T.mouseName '_Raster']);
