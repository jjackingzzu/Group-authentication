% clc;clear all;
% K=5;
%% 
Constellation_user1=[-1,1];
Constellation_user2=[-1i,1i];
ii=0;
for i1=1:2
    for i2=1:2
        ii=ii+1;
        Constellation_table_2(ii,1:2)=[Constellation_user1(i1),Constellation_user2(i2)];
                Constellation_table_2(ii,3)=sum(Constellation_table_2(ii,1:2));
                Constellation_index_2(ii,1:2)=[i1,i2]-1;
    end
end

%% K=3
Constellation_user1=[-2-1i,1i];
Constellation_user2=[-1,1-2*1i];
Constellation_user3=[1,1+2*1i];
ii=0;
for i1=1:2
    for i2=1:2
        for i3=1:2
                ii=ii+1;
                Constellation_table_3(ii,1:3)=[Constellation_user1(i1),Constellation_user2(i2),...
                    Constellation_user3(i3)];
                Constellation_table_3(ii,4)=sum(Constellation_table_3(ii,1:3));
                Constellation_index_3(ii,1:3)=[i1,i2,i3]-1;
        end
    end
end

%% K=4
Constellation_user1=[-1,-2-1i*2];
Constellation_user2=[-2,2-1i];
Constellation_user3=[1+2*1i,-2*1i];
Constellation_user4=[2*1i,2+1i];
ii=0;
for i1=1:2
    for i2=1:2
        for i3=1:2
            for i4=1:2
                ii=ii+1;
                Constellation_table_4(ii,1:4)=[Constellation_user1(i1),Constellation_user2(i2),...
                    Constellation_user3(i3),Constellation_user4(i4)];
                Constellation_table_4(ii,5)=sum(Constellation_table_4(ii,1:4));
                Constellation_index_4(ii,1:4)=[i1,i2,i3,i4]-1;
            end
        end
    end
end
clear i1 i2 i3 i4 ii;
clear Constellation_user4 Constellation_user3 Constellation_user2 Constellation_user1;

%% K=5
Constellation_user1=[-3,-1i];
Constellation_user2=[-1-2*1i,1+2*1i];
Constellation_user3=[-1+2*1i,1-2*1i];
Constellation_user4=[-1i,3*1i];
Constellation_user5=[-1i,3];
ii=0;
for i1=1:2
    for i2=1:2
        for i3=1:2
            for i4=1:2
                for i5=1:2
                ii=ii+1;
                Constellation_table_5(ii,1:5)=[Constellation_user1(i1),Constellation_user2(i2),...
                    Constellation_user3(i3),Constellation_user4(i4) ,Constellation_user5(i5)];
                Constellation_table_5(ii,6)=sum(Constellation_table_5(ii,1:5));
                Constellation_index_5(ii,1:5)=[i1,i2,i3,i4,i5]-1;
                end
            end
        end
    end
end
clear i i1 i2 i3 i4 i5 ii Constellation_user5 Constellation_user4;
clear Constellation_user3 Constellation_user2 Constellation_user1;
if K==2 
    Constellation_table=Constellation_table_2;
    Constellation_index=Constellation_index_2;
end
if K==3 
    Constellation_table=Constellation_table_3;
    Constellation_index=Constellation_index_3;
end
if K==4 
    Constellation_table=Constellation_table_4;
    Constellation_index=Constellation_index_4;
end
if K==5 
    Constellation_table=Constellation_table_5;
    Constellation_index=Constellation_index_5;
end
clear Constellation_table_5 Constellation_table_4 Constellation_table_3 Constellation_table_2;
clear Constellation_index_2 Constellation_index_3 Constellation_index_4 Constellation_index_5;