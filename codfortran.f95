! Material Suplementar para “Evoluc¸ao temporal da func¸ ˜ ao de distribuic¸ ˜ ao de sistemas ˜ unidimensionais”
! Codigo em FORTRAN para resolver a equac¸ ´ ao de Vlasov unidimensional ˜
! USAMOS O M´ETODO LIE-SPLITTING PARA SEPARAR A EQUAC¸˜AO DE VLASOV EM DUAS EQUAC¸˜OES DE TRANSPORTE.
! RESOLVEMOS A EQUAC¸˜AO DE TRANSPORTE USANDO O M´ETODO PFC DE SEGUNDA ORDEM COM LIMITADOR DE
! INCLINAC¸˜AO MC (monotonized central-difference) LIMITER.

PROGRAM vlasov1d
    USE OMP_LIB

    IMPLICIT NONE


INTEGER :: I,J,K,dK
REAL*8, DIMENSION(:), ALLOCATABLE :: X, Y, entropia, massa, ek, ep
REAL*8, DIMENSION(:,:), ALLOCATABLE :: F, f1, den, pot
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: FF ! FD
REAL*8, PARAMETER :: pi = 4.0d0*ATAN(1.0d0)
REAL*8, PARAMETER :: Lx = 8.0d0, Ly = 8.0d0 ! DETERMINA O DOM´INIO NO ESPAC¸O DE FASES
INTEGER, PARAMETER :: Nt = 500, Nx = 800, Ny = 200 ! n´umero de n´ıveis (Nt) e n´umero de celulas (Nx,Ny)
REAL*8, PARAMETER :: dt = 0.002, dX = Lx/Nx, dY = Ly/Ny ! tamanhos de passo
INTEGER, PARAMETER :: Nn = 10, dNt = Nt/Nn
! Nn = n´umero de macro-n´ıveis para imprimir a FD, dNt = n´umero de n´ıveis
! dentro de um macro-n´ıvel
WRITE(*,FMT='(2X,"__________________________________________________________________")')
WRITE(*,FMT='(8X," EQUACAO DE VLASOV 1D com Splitting e PFC")')
WRITE(*,FMT='(2X,"__________________________________________________________________")')
WRITE(*,FMT='(8X,"Intervalo de tempo para imprimir os dados: ",2X,F10.5)')dNt*dt
WRITE(*,FMT='(8X,"Tempo Final: ",2X,F9.3)')Nt*dt
WRITE(*,FMT='(8X,"dX: ",2X,F8.5,8X,"Courant number x: ",2X,F10.5)')dX,dt/dX
WRITE(*,FMT='(8X,"dY: ",2X,F8.5,8X,"Courant number y: ",2X,F10.5)')dY,dt/dY
WRITE(*,*)


ALLOCATE(X(Nx),Y(Ny),F(Nx,Ny),f1(Nx,Ny))
ALLOCATE(FF(0:Nn,Nx,Ny))
ALLOCATE(entropia(0:Nt),massa(0:Nt),ek(0:Nt),ep(0:Nt))
ALLOCATE(den(0:Nt,Nx),pot(0:Nt,Nx))

forall(I = 1:Nx) X(I) = (-Lx/2 + dX/2) + (I-1)*dX
forall(J = 1:Ny) Y(J) = (-Ly/2 + dY/2) + (J-1)*dY


DO I = 1,Nx; DO J = 1,Ny
FF(0,I,J) = finicial(X(I),Y(J)) ! condic¸˜ao inicial da FD
END DO; END DO


F(:,:) = FF(0,:,:)

FORALL (I=1:Nx) den(0,I) = SUM(F(I,:))*dY
massa(0) = sum(den(0,:))*dX
CALL potencial_n(den(0,:),pot(0,:))
ep(0) = 0.5*dot_product(pot(0,:),den(0,:))*dX
CALL energia_k(F,Y,dX,dY,ek(0))
CALL entropy(F,dX,dY,entropia(0))

DO K=0,Nt-1
    CALL lie_splitting(F,Y,dt,dX,dY,den(K+1,:),pot(K+1,:))
    massa(K+1) = sum(den(K+1,:))*dX
    CALL entropy(F,dX,dY,entropia(K+1))
    ep(K+1) = 0.5*dot_product(pot(K+1,:),den(K+1,:))*dX
    CALL energia_k(F,Y,dX,dY,ek(K+1))
    IF(MOD(K+1,dNT)==0) FF((K+1)/dNT,:,:) = F(:,:)
    IF(MOD(K+1,50)==0) THEN
    WRITE(*,FMT='(8X,"ITERACAO: ",2X,I9,"TEMPO: ",2X,F10.5,5X,5X,"massa: ",2X,F10.7)')&
    &K+1, (K+1)*dt,massa(K+1)
    END IF
END DO


WRITE(*,*)
WRITE(*,FMT='(8X,"TEMPOS EM QUE A FD SERA IMPRIMIDA: ")')
DO K = 1,Nn; WRITE(*,FMT='(8X,F8.3)')K*dNt*dt; END DO

! IMPRESS˜AO DE DADOS
DO K=0,Nt
    WRITE(11,FMT='(5(2X,F24.16))')K*dt,massa(K),massa(K)/massa(0)-1.0d0,entropia(K),&
    &(entropia(K)-entropia(0))/ABS(entropia(0))
END DO
DO I=1,Nx
    DO J=1,Ny
        WRITE(12,FMT='(13(2X,F16.10))')X(I),Y(J),FF(0,I,J),FF(1,I,J),FF(2,I,J),FF(3,I,J),FF(4,I,J),&
        & FF(5,I,J),FF(6,I,J),FF(7,I,J),FF(8,I,J),FF(9,I,J),FF(10,I,J)
    END DO
    WRITE(12,*)
END DO


DO K=0,Nt
    WRITE(14,FMT = '(4(2X,F24.16))')K*dt,ek(K),ep(K),ek(K)+ep(K) !energias do sistema
END DO
dk = 4 ! tamanho de passo para n˜ao sobrecarregar a impress˜ao quando Nt
!seja grande
DO K=0,Nt,dK
    DO I = 1,Nx,8
        WRITE(15,FMT='(4(2X,F16.10))')K*dt,X(I),den(K,I),pot(K,I)
    END DO
    WRITE(15,*)""
END DO

DO I = 1,Nx,20 ! impress˜ao da fd em func¸˜ao da energia orbital
    DO J = 1,Ny,20
        WRITE(16,FMT='(12(2X,F16.10))')0.5*Y(J)**2+pot(0,I),FF(0,I,J),0.5*Y(J)**2+pot(dNt,I),&
        & FF(1,I,J),0.5*Y(J)**2+pot(2*dNt,I),FF(2,I,J),0.5*Y(J)**2+pot(3*dNt,I),FF(3,I,J),&
        & 0.5*Y(J)**2+pot(4*dNt,I),FF(4,I,J),0.5*Y(J)**2+pot(Nn*dNt,I),FF(Nn,I,J)
    END DO
END DO

DO I = 1,Nx,2
    WRITE(17,FMT='(12(2X,F16.10))')X(I),den(0,I),den(dNt,I),den(2*dNt,I),den(3*dNt,I),den(4*dNt,I), &
    & den(5*dNt,I),den(6*dNt,I),den(7*dNt,I),den(8*dNt,I),den(9*dNt,I),den(10*dNt,I)
    WRITE(18,FMT='(12(2X,F16.10))')X(I),pot(0,I),pot(dNt,I),pot(2*dNt,I),pot(3*dNt,I),&
    & pot(4*dNt,I),pot(5*dNt,I),pot(6*dNt,I),pot(7*dNt,I), pot(8*dNt,I),pot(9*dNt,I),pot(10*dNt,I)
END DO


!=======================================================================
CONTAINS
FUNCTION finicial(x1,y1)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: x1,y1
    REAL*8 :: finicial
    finicial = senoidal(x1,0.0d0,2.0d0)*senoidal(y1,0.0d0,1.0d0) ! g_I
    !finicial = quadrada(x1,0.0d0,2.0d0)*quadrada(y1,0.0d0,1.0d0) ! g_II
    !finicial = 0.5*(senoidal(x1,1.0d0,1.0d0)*senoidal(y1,0.0d0,1.0d0)+&
    ! &senoidal(x1,-1.0d0,1.0d0)*senoidal(y1,0.0d0,1.0d0)) ! g_III
    !finicial = disco_senoidal(x1,y1,1.0d0) ! g_IV
END FUNCTION finicial

FUNCTION quadrada(x1,centerq,wideq) !normalizada
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: x1, centerq, wideq; REAL*8 :: quadrada
    IF(ABS(x1-centerq)<=wideq/2) THEN
    quadrada = 1.0d0/wideq
    ELSE; quadrada = 0.0d0; END IF
END FUNCTION quadrada

FUNCTION senoidal(x1,centerq,wideq) !normalizada
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: x1, centerq, wideq; REAL*8 :: senoidal
    IF(ABS(x1-centerq)<=wideq/2) THEN
    senoidal = Cos(pi*(x1-centerq)/wideq)*pi/(2*wideq)
    ! o fator pi/(2*wideq) ´e para normalizar a func¸˜ao
    ELSE; senoidal = 0.0d0; END IF
END FUNCTION senoidal

FUNCTION disco_senoidal(x1,y1,rq) !normalizada
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: x1, y1, rq; REAL*8 :: disco_senoidal
    IF(norma(x1,y1)<=rq) THEN
    disco_senoidal = Cos(pi*norma(x1,y1)/(2*rq))*pi/(4*(pi-2)*rq**2)
    ELSE; disco_senoidal = 0.0d0; END IF
END FUNCTION disco_senoidal

FUNCTION norma(x1,y1)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: x1,y1; REAL*8 :: norma
    norma = SQRT(x1**2+y1**2)
END FUNCTION norma

FUNCTION kernelc(Iq) !kernel de tipo gravitacional
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Iq; REAL*8 :: kernelc
    REAL*8, PARAMETER :: epsq = 2.0d-1
    IF(Iq == 0) THEN
    kernelc = -2*LOG((dX+2*epsq)/(2*epsq))
    ELSE
    kernelc = -LOG((ABS(Iq)*dX + dX/2 + epsq)/(ABS(Iq)*dX - dX/2 + epsq))
    END IF
END FUNCTION kernelc

FUNCTION slope(x1,y1,z1)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: x1, y1, z1; REAL*8 :: slope
    REAL*8 :: delta_m, delta_p, delta_c
    delta_p = z1-y1
    delta_m = y1-x1
    delta_c = z1-x1
    slope = minmod(delta_c/2,minmod(2*delta_m,2*delta_p)) ! MC limiter
END FUNCTION slope

FUNCTION minmod(x1,y1)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: x1, y1; REAL*8 :: minmod
    IF(x1*y1>1.0d-20) THEN
    minmod = min(ABS(x1),ABS(y1))*x1/ABS(x1)
    ELSE
    minmod = 0.0d0
    END IF
END FUNCTION minmod

SUBROUTINE entropy(Fq,dXq,dYq,entropiaq) ! Funcional H
    REAL*8, DIMENSION(:,:), INTENT(IN) :: Fq
    REAL*8, INTENT (IN):: dXq, dYq
    REAL*8, INTENT (OUT):: entropiaq
    INTEGER :: N1q, N2q, Iq, Jq
    REAL*8 :: aux01q
    N1q = SIZE(Fq,1)
    N2q = SIZE(Fq,2)
    aux01q = 0.0d0
    DO Iq=1,N1q; DO Jq=1,N2q
    IF (Fq(Iq,Jq)>0.0d0) aux01q = aux01q - Fq(Iq,Jq)*LOG(Fq(Iq,Jq))
    END DO; END DO
    entropiaq = aux01q*dXq*dYq
END SUBROUTINE entropy

SUBROUTINE energia_k(Fq,Yq,dXq,dYq,ekq) ! energia cin´etica
    REAL*8, DIMENSION(:,:), INTENT(IN) :: Fq
    REAL*8, DIMENSION(:), INTENT(IN) :: Yq
    REAL*8, INTENT (IN):: dXq, dYq
    REAL*8, INTENT (OUT):: ekq
    INTEGER :: Nq, Jq
    REAL*8 :: aux01q
    Nq = SIZE(Fq,2)
    aux01q = 0.0d0
    DO Jq = 1,Nq
    aux01q = aux01q + (0.5*Yq(Jq)**2)*SUM(Fq(:,Jq))*dXq
    END DO
    ekq = aux01q*dYq
END SUBROUTINE energia_k

SUBROUTINE potencial_n(denq,potnq)
    REAL*8, DIMENSION(:), INTENT(IN) :: denq
    REAL*8, DIMENSION(:), INTENT(OUT) :: potnq
    INTEGER :: Nq, Iq, Jq
    REAL*8 :: aux01q
    Nq = SIZE(denq)
    DO Iq = 1,Nq
    aux01q = 0.0d0

    DO Jq = 1,Nq
    aux01q = aux01q + kernelc(Iq-Jq)*denq(Jq)
    END DO
    potnq(Iq) = aux01q
    END DO
END SUBROUTINE potencial_n

SUBROUTINE campo_n(potnq,dXq,camponq)
    REAL*8, DIMENSION(:), INTENT(IN) :: potnq
    REAL*8, INTENT(IN) :: dXq
    REAL*8, DIMENSION(:), INTENT(OUT) :: camponq
    INTEGER :: Nq, Iq
    Nq = SIZE(potnq)
    camponq(1) = -(potnq(2)-potnq(Nq))/(2*dXq)
    camponq(Nq) = -(potnq(1)-potnq(Nq-1))/(2*dXq)
    forall (Iq=2:Nq-1) camponq(Iq) = -(potnq(Iq+1)-potnq(Iq-1))/(2*dXq)
END SUBROUTINE campo_n

SUBROUTINE PFC_2_MC(sq,Fq,F2q) ! pfc de segunda ordem com condic¸ao de contorno PERI´ODICA
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: sq
    REAL*8, DIMENSION(:), INTENT(IN) :: Fq
    REAL*8, DIMENSION(:), INTENT(OUT) :: F2q
    REAL*8, DIMENSION(:), ALLOCATABLE :: phip
    INTEGER :: Nq, Iq, Jq, Jm1, J2q, Jp1, Kq, K2q
    REAL*8 :: alfaq, aux01q
    Nq = SIZE(Fq)
    ALLOCATE(phip(Nq))
    if (ABS(sq)<1.0d-12) then ! para um n´umero de Courant muito pequeno
    FORALL(Iq=1:Nq) F2q(Iq) = Fq(Iq)
    else
    DO Iq=1,Nq
    Jq = floor(Iq+1-sq)
    if (sq>0) then
    alfaq = Jq-Iq+sq
    J2q = Jq; IF(J2q<1) J2q = J2q + Nq
    Jm1 = J2q-1; IF(Jm1<1) Jm1 = Jm1 + Nq
    Jp1 = J2q+1; IF(Jp1>Nq) Jp1 = Jp1 - Nq
    aux01q = 0.0d0
    DO Kq=Jq+1,Iq; K2q = Kq; IF(K2q<1) K2q = K2q + Nq; aux01q = aux01q + Fq(K2q); END DO
    phip(Iq) = alfaq*(Fq(J2q) + 0.5*slope(Fq(Jm1),Fq(J2q),Fq(Jp1))*(1-alfaq)) + aux01q
    else
    alfaq = Jq-Iq-1+sq
    J2q = Jq; IF(J2q>Nq) J2q = J2q - Nq
    Jm1 = J2q-1; IF(Jm1<1) Jm1 = Jm1 + Nq
    Jp1 = J2q+1; IF(Jp1>Nq) Jp1 = Jp1 - Nq
    aux01q = 0.0d0
    DO Kq=Iq+1,Jq-1; K2q = Kq; IF(K2q>Nq) K2q = K2q - Nq; aux01q = aux01q + Fq(K2q); END DO
    phip(Iq) = alfaq*(Fq(J2q) - 0.5*slope(Fq(Jm1),Fq(J2q),Fq(Jp1))*(1+alfaq)) - aux01q
    end if
    END DO
    F2q(1) = Fq(1) + phip(Nq) - phip(1) ! v´alido somente para c. de contorno PERI´ODICA
    FORALL(Iq=2:Nq) F2q(Iq) = Fq(Iq) + phip(Iq-1) - phip(Iq)
    end if

    deallocate(phip)
END SUBROUTINE PFC_2_MC


SUBROUTINE lie_splitting(Fq,Yq,dtq,dXq,dYq,denq,potq)
    REAL*8, DIMENSION(:,:), INTENT(INOUT) :: Fq
    REAL*8, DIMENSION(:), INTENT(IN) :: Yq
    REAL*8, INTENT(IN) :: dtq, dXq, dYq
    REAL*8, DIMENSION(:), INTENT(OUT) :: denq, potq
    REAL*8, DIMENSION(:,:), allocatable :: F2q
    REAL*8, DIMENSION(:), allocatable :: campoq
    INTEGER :: N1q, N2q, Iq, Jq
    REAL*8 :: sq
    N1q = SIZE(Fq,1)
    N2q = SIZE(Fq,2)
    ALLOCATE(F2q(N1q,N2q),campoq(N1q))
    !$OMP PARALLEL DO
    DO Jq=1,N2q ! (I) transporta a FD somente na direc¸˜ao de x
    sq = Yq(Jq)*dtq/dXq
    CALL PFC_2_MC(sq,Fq(:,Jq),F2q(:,Jq))
    !CALL PFC_3(sq,Fq(:,Jq),F2q(:,Jq))
    END DO
    !$OMP END PARALLEL DO
    FORALL (Iq=1:N1q) denq(Iq) = SUM(Fq(Iq,:))*dYq ! densidade de massa
    CALL potencial_n(denq,potq)
    CALL campo_n(potq,dXq,campoq)
    !$OMP PARALLEL DO
    DO Iq=1,N1q ! (II) transporta a FD somente na direc¸˜ao de v
    sq = campoq(Iq)*dtq/dYq
    CALL PFC_2_MC(sq,F2q(Iq,:),Fq(Iq,:))
    !CALL PFC_3(sq,F2q(Iq,:),Fq(Iq,:))
    END DO
    !$OMP END PARALLEL DO
    DEALLOCATE(F2q,campoq)
END SUBROUTINE lie_splitting
END PROGRAM vlasov1d
