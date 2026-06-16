function K = sim_zt(AllD,DT,f)

K1111_i=Kijkl(AllD,DT,f,1,1,1,1);%1
K2222_i=Kijkl(AllD,DT,f,2,2,2,2);%2
K3333_i=Kijkl(AllD,DT,f,3,3,3,3);%3
K1112_i=Kijkl(AllD,DT,f,1,1,1,2);%4
K1113_i=Kijkl(AllD,DT,f,1,1,1,3);%5
K1222_i=Kijkl(AllD,DT,f,1,2,2,2);%6
K2223_i=Kijkl(AllD,DT,f,2,2,2,3);%7
K1333_i=Kijkl(AllD,DT,f,1,3,3,3);%8
K2333_i=Kijkl(AllD,DT,f,2,3,3,3);%9
K1122_i=Kijkl(AllD,DT,f,1,1,2,2);%10
K1133_i=Kijkl(AllD,DT,f,1,1,3,3);%11
K2233_i=Kijkl(AllD,DT,f,2,2,3,3);%12
K1123_i=Kijkl(AllD,DT,f,1,1,2,3);%13
K1322_i=Kijkl(AllD,DT,f,1,3,2,2);% new 14
K1233_i=Kijkl(AllD,DT,f,1,2,3,3);%15
K1212_i=Kijkl(AllD,DT,f,1,2,1,2);% new 16
K1313_i=Kijkl(AllD,DT,f,1,3,1,3);% new 17
K2323_i=Kijkl(AllD,DT,f,2,3,2,3);% new 18
K1223_i=Kijkl(AllD,DT,f,1,2,2,3);% new 19 (old 14)
K1323_i=Kijkl(AllD,DT,f,1,3,2,3);% new 20
K1213_i=Kijkl(AllD,DT,f,1,2,1,3);% new 21


K =[K1111_i K2222_i K3333_i K1112_i K1113_i ...
    K1222_i K2223_i K1333_i K2333_i K1122_i ...
    K1133_i K2233_i K1123_i  K1322_i ...
    K1233_i K1212_i K1313_i K2323_i ... 
    K1223_i K1323_i K1213_i...
];