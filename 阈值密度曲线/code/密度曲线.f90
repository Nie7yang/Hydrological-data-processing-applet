    !******************************************************************
    !************************  长安大学 聂启阳  ************************
    !************************    2019.01.06    ************************
    !******************************************************************
    MODULE HeWangMiDu
    REAL(4),allocatable :: HWMiDu(:)
    ENDMODULE HeWangMiDu

    program MiDuQuXian
    
    USE HeWangMiDu
    Implicit None

    INTEGER hang,lie,i,j,k,N
    INTEGER(8) LiuYuN,GC
    REAL(8) MianJi
    REAL(4) GCC
    REAL(4),allocatable :: DEM(:,:),PoChang(:,:)
    REAL(4),allocatable :: HeChang(:)
    INTEGER,allocatable :: ACC(:,:),DIR(:,:),YuZhi(:)
    CHARACTER DATE

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
    
    N=700  !等距离点数

    allocate (DEM(hang,lie),DIR(hang,lie),ACC(hang,lie),PoChang(hang,lie),YuZHi(N),HeChang(N),HWMiDu(N))

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

    
    YuZhi=(/(i,i=10,7000,10)/)

    HeChang=0.0
    DO j=1,lie
        DO i=1,hang
            DO k=1,N
                if(ACC(i,j)>=YuZhi(k))then
                    if(ANY((/1,4,16,64/)==DIR(i,j)))then
                        HeChang(k)=HeChang(k)+1.0
                    else
                        HeChang(k)=HeChang(k)+sqrt(2.0)
                    endif
                else
                    exit
                endif
            ENDDO
        ENDDO
    ENDDO

    GCC=real(GC)/1000.0
    HwMiDu=real(HeChang)/real(LiuYuN*GCC)  !河网密度单位为：Km/(Km^2)

    !DO i=1,N
    !    write(444,*)YuZhi(i),',',HWMiDu(i)
    !ENDDO
    
    
    call JunZhiBianDian(N)
    
    write(*,*)'结束'
    close(111)
    close(222)
    close(333)
    close(444)
    DEALLOCATE (DEM,DIR,ACC,PoChang,YuZhi,HeChang,HWMiDu)
    end program
