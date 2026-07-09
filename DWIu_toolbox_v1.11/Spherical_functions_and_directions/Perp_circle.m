function Dir_per = Perp_circle(gnz, nodir)

if(abs(gnz(1))==1);
    Dir_per=Perp_circle_y(gnz(1),gnz(2),gnz(3),nodir);
else
    Dir_per=Perp_circle_x(gnz(1),gnz(2),gnz(3),nodir);
end