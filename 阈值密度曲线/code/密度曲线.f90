    !******************************************************************
    !************************  ������ѧ ������  ************************
    !************************    2019.01.06    ************************
    !******************************************************************

    program MiDuQuXian
    Implicit None

    INTEGER hang,lie,i,j,k
    INTEGER(8) LiuYuN,GC
    REAL(8) MianJi
    REAL(4) GCC
    REAL(4),allocatable :: DEM(:,:),PoChang(:,:)
    REAL(4),allocatable :: HWMiDu(:),HeChang(:)
    INTEGER,allocatable :: ACC(:,:),DIR(:,:),YuZhi(:)
    CHARACTER DATE

    OPEN (111,FILE='dem.txt')
    OPEN (222,FILE='dir.txt')
    OPEN (333,FILE='acc.txt')
    OPEN (444,FILE='���.csv')

    READ (111,*) DATE,lie	!��ȡ������
    READ (111,*) DATE,hang	!��ȡ������
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

    !MianJi=LiuYuN*GC**2        !�������

    YuZhi(1:5)=(/10,20,40,60,80/)
    YuZhi(6:45)=(/(i,i=100,4000,100)/)
    YuZhi(46:50)=(/(i,i=5000,9000,1000)/)

    HeChang=0.0
    DO j=1,lie
        DO i=1,hang
            DO k=1,50
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
    HwMiDu=real(HeChang)/real(LiuYuN*GCC)  !�����ܶȵ�λΪ��Km/(Km^2)

    DO i=1,50
        write(444,*)YuZhi(i),',',HWMiDu(i)
    ENDDO
    write(*,*)'����'
    close(111)
    close(222)
    close(333)
    close(444)
    DEALLOCATE (DEM,DIR,ACC,PoChang,YuZhi,HeChang,HWMiDu)
    end program 
