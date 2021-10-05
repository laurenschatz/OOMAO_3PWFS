%% Plot

x=[0,2,4,6,8,10,12];

%% 3PWFS Gain 1.8
y3RI=[19.203, 19.201, 19.186, 19.12, 18.632, 15.906, 8.764];

y3SM=[19.11, 19.109, 19.093, 19.017, 18.518, 15.737, 8.66];

%%  3PWFS Gain 2.25

y3SM2=[19.81, 19.797, 19.733, 19.425, 17.759, 11.746, 6.622];
y3RI2=[19.912, 19.907, 19.835, 19.633, 17.91, 12.067, 6.851];

%% 4PWFS Gain 1.8

y4RI=[ 20.894, 20.903, 20.894, 20.811,20.361, 17.726, 10.09];
y4SM=[ 20.693, 20.693, 20.679, 20.602, 20.08, 17.489, 9.65];

%%

figure; plot(x,y3SM); hold on;  
plot(x,y3RI);
plot(x, y4SM);
plot(x, y4RI);
plot(x, y3RI2);
plot(x, y3SM2);
l=legend('3PWFS SM loop gain 1.8', '3PWFS RI loop gain 1.8', '4PWFS SM loop gain 1.8', '4PWFS RI gain 1.8', '3PWFS RI loop gain 2.25', '3PWFS SM loop gain 1.8');
l.FontSize=20;
ylabel('Strehl Ratio', 'FontSize', 20)
xlabel('Guide Star Magnitude', 'FontSize', 20)




