% Plotting a square


x = linspace(0,1,100);

plot(x,x.*0,'k-','LineWidth',3)
hold on
plot(x,x.*0+1,'k-','LineWidth',3)
plot(x.*0,x,'k-','LineWidth',3)
plot(x.*0+1,x,'k-','LineWidth',3)
plot(x,x.*0+0.5,'k-','LineWidth',1)
plot(x.*0+0.5,x,'k-','LineWidth',1)
axis([0 1 0 1])
axis equal
box off

text(0.75,0.75,'1','FontSize',30)
text(0.25,0.75,'2','FontSize',30)
text(0.25,0.25,'3','FontSize',30)
text(0.75,0.25,'4','FontSize',30)

xlabel('x');
ylabel('y');
exportfig('init.eps')

hold off