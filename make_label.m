function [] = make_label(x_text,y_text)

set(gca,'XGrid','off','YGrid','off','XMinorTick','off','YMinorTick','off','fontsize', 24,'fontname','Times new roman');

xlabel(x_text,'FontSize',24,'fontname','Times new roman');

ylabel(y_text,'FontSize',24,'fontname','Times new roman');

box on;

grid
