    Program JiangYuChaZhi
    Implicit None

    Integer:: ZhanN,Hang,Lie,XN,YN,Cell,ChaXY(2),X0,Y0,YanXY(2),Mini(1),ios,honghao,daihao,honghao0,bianhao
    Character:: Name
    Character(len=30):: dat(5)
    Integer,Allocatable:: ZhanXY(:,:),DEM(:,:),YiJi(:),LiuJu(:,:)
    Real,Allocatable:: RainXY(:,:),RainT(:),JieDian(:),TongJi(:)
    Real(kind=8),Allocatable:: long(:)
    InTeger:: i,j,m,n,k,Shu
    Real:: Chang,Fenmu,step
    Real(kind=8):: JuLi
    character(len=8)::mingzi1,mingzi2

    open(111,file='elevation.asc')  !高程数据
    open(112,file='秃尾河站点计算.csv')  !站点位置
    Open(222,file='distance.asc')   !汇流距离
    open(666,file='秃尾河雨洪筛选-聂启阳.csv')

    read(111,*)Name,Lie
    read(111,*)Name,Hang
    read(111,*)Name,XN
    read(111,*)Name,YN
    read(111,*)Name,Cell
    read(111,*)

    read(112,*)ZhanN
    Shu=100 !汇流距离分层数
    Allocate(ZhanXY(ZhanN,3),DEM(Lie,Hang),RainXY(Lie,Hang),LiuJu(Lie,Hang),RainT(ZhanN),Long(ZhanN+1),YiJi&
        &(ZhanN),JieDian(0:Shu),TongJi(Shu))

    read(111,*)DEM

    DO i=1,ZhanN
        read(112,*)Name,ZhanXY(i,1),ZhanXY(i,2),ZhanXY(i,3)
    ENDDO
    Close(111)
    Close(112)
    open(111,file='elevation.asc')


    DO i=1,5
        read(111,'(A30)')dat(i)
    ENDDO
    read(111,*)
    close(111)

    ios=0
    honghao0=0
    bianhao=0
    DO  
        read(666,*，iostat=ios)honghao,daihao,daihao,daihao,daihao,daihao,RainT
		if(ios/=0)exit
        if(honghao==honghao0)then
            bianhao=bianhao+1
        else
            bianhao=1
            honghao0=honghao
        endif
        write(mingzi1,'(i8)')honghao
        write(mingzi2,'(i4)')bianhao
        if(all(rainT<=0.0))cycle

        open(113,file=trim(adjustl(mingzi1))//trim(adjustl(mingzi2))//'.asc')
        DO i=1,5
            write(113,'(A30)')dat(i)
        ENDDO
        write(113,*)'NODATA_value   -99.00'

        !RainT=(/0.0,-999.0,0.0,0.683,2.917,0.0,-999.0,3.3,2.667,8.9,7.903/) !测试数据，降雨后期直接输给插值函数
        RainXY=-99

        X0=XN+cell/2
        Y0=YN+cell/2
        step=cell/2.0
        DO j=1,hang
            write(*,*)j,'/',hang
            DO i=1,lie  !遍历每一个点进行插值操作
                if(DEM(i,j)>=0)then
                    RainXY(i,j)=0.0
                    ChaXY(1)=X0+(i-1)*cell          !插值点X坐标
                    ChaXY(2)=Y0+(hang-j)*cell       !插值点Y坐标
                    YiJi=0                          !记录一级插值站
                    DO n=1,hang
                        DO m=1,lie
                            if(DEM(m,n)>=0)then
                                YanXY(1)=X0+(m-1)*cell          !验证点X坐标
                                YanXY(2)=Y0+(hang-n)*cell       !验证点Y坐标
                                Long=9999999
                                DO k=1,ZhanN
                                    if(RainT(k)<-1.0)then
                                        cycle      !若降雨数值为负数，则为缺测点，不计入筛选范围
                                    endif
                                    long(k)=JuLi(YanXY(1),YanXY(2),DEM(m,n),ZhanXY(k,1),ZhanXY(k,2),ZhanXY(k,3))
                                ENDDO
                                long(ZhanN+1)=JuLi(YanXY(1),YanXY(2),DEM(m,n),ChaXY(1),ChaXY(2),DEM(i,j))
                                Mini=Minloc(long)
                                if(Mini(1)==ZhanN+1)then
                                    long(ZhanN+1)=9999999
                                    MiNi=Minloc(long)
                                    YiJi(MiNi(1))=1
                                endif
                            endif
                        ENDDO
                    ENDDO
                    n=0
                    Fenmu=0.0
                    DO k=1,ZhanN
                        if(YiJi(k)==1)then
                            n=n+1
                            Chang=JuLi(ChaXY(1),ChaXY(2),DEM(i,j),ZhanXY(k,1),ZhanXY(k,2),ZhanXY(k,3))
                            if(Chang<step)then
                                n=0
                                RainXY(i,j)=RainT(k)
                                exit
                            else
                                RainXY(i,j)=RainXY(i,j)+(RainT(k)*10000)/(Chang**2)
                                Fenmu=Fenmu+10000/(Chang**2)
                            endif
                        endif
                    ENDDO
                    if(n>=1)then
                        RainXY(i,j)=RainXY(i,j)/Fenmu
                    endif
                endif
            ENDDO
        ENDDO

        DO j=1,hang
            write(113,'(<lie>f12.4)')(RainXY(i,j),i=1,lie)
        ENDDO
        close(113)
    ENDDO
    
    !统计空间降雨的汇流距离分布
    !Do i=1,6
    !    read(222,*)
    !END DO
    !read(222,*)LiuJu
    !close(222)
    !step=maxval(LiuJu)/Shu
    !JieDian(0)=0.0
    !Do i=1,Shu
    !    JieDian(i)=JieDian(i-1)+step
    !ENDDO
    !TongJi=0.0
    !DO j=1,hang
    !    DO i=1,lie
    !        if(DEM(i,j)>=0)then
    !            DO k=1,Shu
    !                if(LiuJu(i,j)>=JieDian(k-1).and.LiuJu(i,j)<JieDian(k))then
    !                    TongJi(k)=TongJi(k)+RainXY(i,j)
    !                endif
    !            ENDDO
    !        endif
    !    ENDDO
    !ENDDO
    !Open(222,file='统计.csv')
    !DO i=1,Shu
    !    write(222,*)JieDian(i),',',TongJi(i)
    !ENDDO
    !close(222)

    Deallocate(ZhanXY,DEM,RainXY,RainT,Long,YiJi,LiuJu,TongJi)
    
    EndProgram

    Real(kind=8) Function JuLi(X0,Y0,H0,X1,Y1,H1)
    Implicit None

    Integer,Intent(in):: X0,Y0,H0,X1,Y1,H1
    real(kind=8):: Long

    if((abs(X0-X1)<=1).and.(abs(Y0-Y1)<=1))then
        Long=0.0
    else
        long=((X0-X1)/100.0)**2
        long=long+((Y0-Y1)/100.0)**2
        long=long+((H0-H1)/100.0)**2
        long=100.0*sqrt(long)
    endif
    Juli=long

    End Function Juli

