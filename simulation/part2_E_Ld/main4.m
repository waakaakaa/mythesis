subplot(2,1,1);
plot(Ld_plot, Ec(:,1), '-o');hold on;
plot(Ld_plot, Ec(:,2), '-o');hold on;
plot(Ld_plot, Ec(:,3), '-o');hold on;
subplot(2,1,2);
plot(Ld_plot, Ehh(:,1), 'r-o');hold on;
plot(Ld_plot, Ehh(:,2), 'r-o');hold on;
plot(Ld_plot, Ehh(:,3), 'r-o');hold on;

plot(Ld_plot, Elh(:,1), 'g-o');hold on;
plot(Ld_plot, Elh(:,2), 'g-o');hold on;