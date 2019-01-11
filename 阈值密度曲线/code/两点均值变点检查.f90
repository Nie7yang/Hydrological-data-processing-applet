    subroutine JunZhiBianDian_2(N)
    use HeWangMiDu
    implicit none

    INTEGER N,i,K,G,M
    REAL(8) JunZhi,S,Xi1,Xi2,Xi3,Si1,Si2,Si3
    REAL(8),ALLOCATABLE:: Si(:,:)

    OPEN (123,FILE='均值变点检查.txt')
    
    ALLOCATE(Si(2:N-1,3:N))
    JunZhi=SUM(HWMiDu)/REAl(N)
    DO i=1,N
        S=S+(HWMiDu(i)-JunZhi)**2
    ENDDO
    Si=-9999
    DO K=2,N-1
        DO G=K+1,N
            Xi1=SUM(HWMiDu(1:K-1))/real(K-1)
            Xi2=SUM(HWMiDu(K:G-1))/real(G-K)
            Xi3=SUM(HWMiDu(G:N))/real(N-G+1)
            Si1=0.0
            Si2=0.0
            si3=0.0
            DO M=1,K-1
                Si1=Si1+(HWMiDu(M)-Xi1)**2
            ENDDO
            DO M=K,G-1
                Si2=Si2+(HWMiDu(M)-Xi2)**2
            ENDDO
            DO M=G,N
                Si3=Si3+(HWMiDu(M)-Xi3)**2
            ENDDO
            Si(K,G)=Si1+Si2+Si3
        ENDDO       
    ENDDO

    !DO K=2,N-1
    !    DO G=K+1,N
    !        write(123,*)K,',',G,',',Si(K,G)
    !    ENDDO
    !ENDDO
    
    do K = 2,N-1
        write(123,'(<N-2>F14.4)') (Si(K,G), G = 3,N)
    enddo
    
    DEALLOCATE(Si)
    close(123)
    
    return
    end subroutine JunZhiBianDian_2