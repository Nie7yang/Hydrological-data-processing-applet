    !******************************************************************
    !************************  长安大学 聂启阳  ************************
    !************************    2019.01.06    ************************
    !******************************************************************
    Module ran_mod
    !生成正态分布的模块
    Implicit None
    ! ran 返回回0-1之间的均匀随机数
    ! norma 返回正态分布值
    contains
    function ran()   ! ran 返回回0-1之间的均匀随机数
    implicit none
    integer , save :: flag = 0
    double precision :: ran
    if(flag==0) then
        call random_seed()
        flag = 1
    endif
    call random_number(ran)     
    end function ran

    function normal(mean,sigma)  ! norma 返回正态分布值
    implicit none
    integer :: flag
    double precision, parameter :: pi = 3.141592653589793239
    double precision :: u1, u2, y1, y2, normal, mean, sigma
    save flag
    data flag /0/
    u1 = ran(); u2 = ran()
    if (flag.eq.0) then
        y1 = sqrt(-2.0d0*log(u1))*cos(2.0d0*pi*u2)
        normal = mean + sigma*y1
        flag = 1
    else
        y2 = sqrt(-2.0d0*log(u1))*sin(2.0d0*pi*u2)
        normal = mean + sigma*y2
        flag = 0
    endif
    end function normal
    End Module ran_mod

    program MiDuQuXian
    use ran_mod
    Implicit None

    INTEGER hang,lie,i,j,k
    INTEGER(8) LiuYuN,GC
    REAL(8) MianJi
    REAL(4) GCC
    REAL(4),allocatable :: DEM(:,:),PoChang(:,:)
    REAL(4),allocatable :: HWMiDu(:)
    INTEGER,allocatable :: ACC(:,:),DIR(:,:),YuZhi(:),HeChang(:)
    CHARACTER DATE

    Integer , parameter :: N = 40
    Real( Kind = 8 ) :: a( N )

    OPEN (111,FILE='dem.txt')
    OPEN (222,FILE='dir.txt')
    OPEN (333,FILE='acc.txt')
    OPEN (444,FILE='输出.csv')

    READ (111,*) DATE,lie	!读取总列数
    READ (111,*) DATE,hang	!读取总行数
    READ (111,*)
    READ (111,*)
    READ (111,*) DATE,GC
    READ (111,*)

    DO i=1,6
        READ (222,*)
        READ (333,*)
    ENDDO

    allocate (DEM(hang,lie),DIR(hang,lie),ACC(hang,lie),PoChang(hang,lie),YuZHi(50),HeChang(50),HWMiDu(50))

    DO I=1,hang
        READ (111,*) (DEM(i,j),J=1,lie)
        READ (222,*) (DIR(i,j),J=1,lie)
        READ (333,*) (ACC(i,j),J=1,lie)
    ENDDO

    LiuYuN=0
    DO j=1,lie
        DO i=1,hang
            if(DEM(i,j)>0)then
                LiuYuN=LiuYuN+1
            endif
        ENDDO
    ENDDO

    !MianJi=LiuYuN*GC**2        !流域面积

    Do i = 1 , N                !前一部分阈值用正态分布保证可能的均匀性，后一部分根据经验赋予主要的阈值
        a( i ) = normal( 2500.0D0 , 1000.0D0 )
        if(a(i)<0)then
            a(i)=-a(i)
        endif
        YuZhi(i)=int(a(i))
    End Do
    YuZhi(41:50)=(/5,10,20,50,100,200,400,600,800,1000/)   
    !DO i=1,50
    !    YuZhi(i)=5+(LiuYuN/100-5)*(i-1)/50
    !ENDDO

    HeChang=0
    DO j=1,lie
        DO i=1,hang
            DO k=1,50
                if(ACC(i,j)>=YuZhi(k))then
                    HeChang(k)=HeChang(k)+1
                !else
                !    exit       ! 因阈值未按照小到大顺序，因此不可提前退出
                endif
            ENDDO
        ENDDO
    ENDDO

    GCC=real(GC)/1000.0
    HwMiDu=real(HeChang)/real(LiuYuN*GCC)

    DO i=1,50
        write(444,*)YuZhi(i),',',HWMiDu(i)
    ENDDO
    write(*,*)'结束'
    close(111)
    close(222)
    close(333)
    close(444)
    DEALLOCATE (DEM,DIR,ACC,PoChang,YuZhi,HeChang,HWMiDu)
    end program 
