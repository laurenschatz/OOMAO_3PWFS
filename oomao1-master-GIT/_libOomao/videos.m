data=fitsread('/Users/laurenschatz/Documents/exao0Pull/cam4p_tMultiplier2_Gain_0.8.fits');

figure; imagesc(data(:,:,1)); axis equal






%%
utilities.videoMaker(data)




