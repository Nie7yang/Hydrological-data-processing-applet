    subroutine JunZhiBianDian(N)
    use HeWangMiDu
    implicit none

    INTEGER N,i,j
    REAL(8) JunZhi,S,Xi1,Xi2,XiZ,XiY,Si1,Si2
    REAL(8) Si(N-1)

    OPEN (123,FILE='均值变点检查.csv')
    S=0.0
    JunZhi=SUM(HWMiDu)/REAl(N)
    DO i=1,N
        S=S+(HWMiDu(i)-JunZhi)**2
    ENDDO
    Si=0.0
    DO i=2,N
        Xi1=SUM(HWMiDu(1:i-1))
        XiZ=Xi1/real(i-1)
        Xi2=SUM(HWMiDu(i:N))
        XiY=Xi2/real(N+1-i)
        Si1=0.0
        Si2=0.0
        DO j=1,i-1
            Si1=Si1+(HWMiDu(j)-XiZ)**2
        ENDDO
        DO j=i,N
            Si2=Si2+(HWMiDu(j)-XiY)**2
        ENDDO
        Si(i-1)=Si1+si2
    ENDDO
    !Si=-(Si-S)

    DO i=1,N-1
        write(123,*)i,',',Si(i),',',(i+1)*10,',',HWMiDu(i+1)
    ENDDO
    close(123)
    
    return
    end subroutine JunZhiBianDian