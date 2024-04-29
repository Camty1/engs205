        PROGRAM comparison
        IMPLICIT REAL*8 (A-H, O-Z)
        
        INTEGER, PARAMETER :: N = 80, k = 4, D = 1
        REAL*8, PARAMETER :: r_min = 0.1, r_max = 1, phi_min = 0
        REAL*8, PARAMETER :: phi_max = 1.570796326794897
        REAL*8, PARAMETER :: sigma = 1.0, I0 = 1.0
        REAL*8, PARAMETER :: theta = 0.0, rich = 0.1
        REAL*8, PARAMETER :: alpha = 1.0
        LOGICAL, PARAMETER :: TYPE_III = .FALSE.
        
        REAL*8, DIMENSION(N) :: r_arr, phi_arr
        REAL*8, DIMENSION(N**2) :: U, U_new, b_vec
        REAL*8, DIMENSION(N**2, 2*N+1) :: A_mat

        INTEGER :: i, j
        REAL*8 :: x, y, dr, dt
        CHARACTER(LEN = 50) :: filename

        CALL LINSPACE(r_min, r_max, N, r_arr)
        CALL LINSPACE(phi_min, phi_max, N, phi_arr)

        dr = r_arr(2)-r_arr(1)
        dt = rich * dr**2 / DBLE(D)

        WRITE(filename, '("output/r", I3.3, ".dat")') N
        OPEN(1, file=filename)

        WRITE(filename, '("output/phi", I3.3, ".dat")') N
        OPEN(2, file=filename)

        WRITE(filename, '("output/x", I3.3, ".dat")') N
        OPEN(3, file=filename)

        WRITE(filename, '("output/y", I3.3, ".dat")') N
        OPEN(4, file=filename)

        DO j=1,N
        DO i=1,N
                WRITE (1,'( *(g0, ",") )') r_arr(i)
                WRITE (2,'( *(g0, ",") )') phi_arr(j)

                CALL POLAR_2_CART(r_arr(i), phi_arr(j), x, y)

                WRITE (3,'( *(g0, ",") )') x
                WRITE (4,'( *(g0, ",") )') y

        END DO
        END DO

        CLOSE(1)
        CLOSE(2)
        CLOSE(3)
        CLOSE(4)

        IF (TYPE_III) THEN
                PRINT *, "TYPE III"
                CALL GEN_SS_EQNS_TYPE_III(N, alpha, r_arr, phi_arr, k, I0, sigma, A_mat, b_vec)
                WRITE(filename, '("output/u_ss_3_", I3.3, ".dat")') N

        ELSE 
                PRINT *, "TYPE I"
                CALL GEN_SS_EQNS_TYPE_I(N, r_arr, phi_arr, k, I0, sigma, A_mat, b_vec)
                
                WRITE(filename, '("output/u_ss_1_", I3.3, ".dat")') N

        END IF

        CALL DSOLVE(3, A_mat, b_vec, N**2, N, N**2, 2*N+1)

        OPEN(1, file=filename)

        DO i=1,N**2
                WRITE(1, '(*(g0, "," ))') b_vec(i)
        END DO

        CLOSE(1)

        END PROGRAM

        SUBROUTINE get_A(r, dr, A)
                REAL*8, INTENT(IN) :: r, dr
                REAL*8, INTENT(OUT) :: A

                A = (r + dr/2)/(r * dr**2)

        END SUBROUTINE

        SUBROUTINE get_B(r, dr, B)
                REAL*8, INTENT(IN) :: r, dr
                REAL*8, INTENT(OUT) :: B

                B = (r - dr/2)/(r * dr**2)

        END SUBROUTINE

        SUBROUTINE get_C(r, dr, dphi, C)
                REAL*8, INTENT(IN) :: r, dr, dphi
                REAL*8, INTENT(OUT) :: C

                C = -2 * (1 / dr**2 + 1 / (r * dphi)**2)

        END SUBROUTINE

        SUBROUTINE get_D(r, dphi, D)
                REAL*8, INTENT(IN) :: r, dphi
                REAL*8, INTENT(OUT) :: D

                D = 1 / (r * dphi) ** 2

        END SUBROUTINE

        SUBROUTINE BC(k, I0, sigma, r_max, phi, f)
                INTEGER, INTENT(IN) :: k
                REAL*8, INTENT(IN) :: I0, sigma, r_max, phi
                REAL*8, INTENT(OUT) :: f
                
                f = - I0 / sigma * r_max * COS(k * phi)
        END SUBROUTINE

        SUBROUTINE get_E(r, dphi, alpha, E)
                REAL*8, INTENT(IN) :: r, dphi, alpha
                REAL*8, INTENT(OUT) :: E

                E = -2*r*dphi*alpha
        END SUBROUTINE

        SUBROUTINE get_g(r, dphi, k, I0, sigma, g)
                INTEGER, INTENT(IN) :: k
                REAL*8, INTENT(IN) :: r, dphi, I0, sigma
                REAL*8, INTENT(OUT) :: g

                g = -2*k*r*I0*dphi/sigma
        END SUBROUTINE
         
        SUBROUTINE GEN_SS_EQNS_TYPE_I(N, r_arr, phi_arr, k, I0, sigma, A_mat, b_vec)
                INTEGER, INTENT(IN) :: N, k
                REAL*8, INTENT(IN) :: I0, sigma

                REAL*8, DIMENSION(N), INTENT(IN) :: r_arr, phi_arr
                REAL*8, DIMENSION(N**2, 2*N+1), INTENT(OUT) :: A_mat
                REAL*8, DIMENSION(N**2), INTENT(OUT) :: B_vec
                
                INTEGER :: i, j, row, col, new_col
                
                REAL*8 :: dr, dphi, A, B, C, D, f, r_max

                r_max = MAXVAL(r_arr)

                dr = r_arr(2) - r_arr(1)
                dphi = phi_arr(2) - phi_arr(1)

                A_mat = 0.0
                B_vec = 0.0

                DO j=1,N
                DO i=1,N

                CALL get_A(r_arr(i), dr, A)
                CALL get_B(r_arr(i), dr, B)
                CALL get_C(r_arr(i), dr, dphi, C)
                CALL get_D(r_arr(i), dphi, D)
                row = N*(j-1)+i
                IF (i == 1 .AND. j == 1) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = C

                        col = N*(j-1)+i+1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = A + B

                        col = N*j+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = 2 * D
                        B_vec(row) = -2 * A * dr * (1 + 1)

                ELSE IF (i == 1 .AND. 1 < j .AND. j < N) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = C

                        col = N*(j-1)+i+1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = A + B

                        col = N*(j-2)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = D

                        col = N*j+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = D
                        B_vec(row) = -2 * A * dr * 1

                ELSE IF (1 < i .AND. i < N .AND. j == 1) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = C

                        col = N*(j-1)+i+1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = A

                        col = N*(j-1)+i-1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = B

                        col = N*j+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = 2 * D

                        B_vec(row) = -2 * r_arr(i) * dphi * D * 1

                ELSE IF (1 < i .AND. i < N .AND. 1 < j .AND. j < N) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = C

                        col = N*(j-1)+i+1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = A

                        col = N*(j-1)+i-1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = B

                        col = N*(j-2)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = D

                        col = N*j+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = D

                ELSE IF (j == N) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = 1
                        B_vec(row) = -1

                ELSE IF (i == N) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        CALL BC(k, I0, sigma, r_max, phi_arr(j), f)
                        A_mat(row, new_col) = 1
                        B_vec(row) = f

                ELSE 
                        PRINT *, "You should not be here"

                END IF

                END DO
                END DO
        END SUBROUTINE

        SUBROUTINE GEN_SS_EQNS_TYPE_III(N, alpha, r_arr, phi_arr, k, I0, sigma, A_mat, b_vec)
                INTEGER, INTENT(IN) :: N, k
                REAL*8, INTENT(IN) :: alpha, I0, sigma

                REAL*8, DIMENSION(N), INTENT(IN) :: r_arr, phi_arr
                REAL*8, DIMENSION(N**2, 2*N+1), INTENT(OUT) :: A_mat
                REAL*8, DIMENSION(N**2), INTENT(OUT) :: B_vec
                
                INTEGER :: i, j, row, col, new_col
                
                REAL*8 :: dr, dphi, A, B, C, D, E, f, g, r_max

                r_max = MAXVAL(r_arr)

                dr = r_arr(2)-r_arr(1)
                dphi = phi_arr(2)-phi_arr(1)

                A_mat = 0.0
                B_vec = 0.0
                DO j=1,N
                DO i=1,N

                CALL get_A(r_arr(i), dr, A)
                CALL get_B(r_arr(i), dr, B)
                CALL get_C(r_arr(i), dr, dphi, C)
                CALL get_D(r_arr(i), dphi, D)

                row = N*(j-1)+i
                IF (i == 1 .AND. j == 1) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = C

                        col = N*(j-1)+i+1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = A + B

                        col = N*j+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = 2 * D
                              
                ELSE IF (i == 1 .AND. 1 < j .AND. j < N) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = C

                        col = N*(j-1)+i+1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = A + B

                        col = N*(j-2)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = D

                        col = N*j+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = D

                ELSE IF (i == 1 .AND. j == N) THEN
                        CALL get_E(r_arr(i), dphi, alpha, E)
                        CALL get_g(r_arr(i), dphi, k, I0, sigma, g)

                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = C + D * E

                        col = N*(j-1)+i+1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = A + B

                        col = N*(j-2)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = 2 * D

                        B_vec(row) = -g * D

                ELSE IF (1 < i .AND. i < N .AND. j == 1) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = C

                        col = N*(j-1)+i+1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = A

                        col = N*(j-1)+i-1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = B

                        col = N*j+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = 2 * D
                
                ELSE IF (1 < i .AND. i < N .AND. 1 < j .AND. j < N) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = C

                        col = N*(j-1)+i+1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = A

                        col = N*(j-1)+i-1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = B

                        col = N*(j-2)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = D

                        col = N*j+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = D

                ELSE IF (1 < i .AND. i < N .AND. j == N) THEN
                        CALL get_E(r_arr(i), dphi, alpha, E)
                        CALL get_g(r_arr(i), dphi, k, I0, sigma, g)

                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = C + D * E

                        col = N*(j-1)+i+1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = A

                        col = N*(j-1)+i-1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = B

                        col = N*(j-2)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = 2 * D

                        B_vec(row) = -g * D

                ELSE IF (i == N) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        CALL BC(k, I0, sigma, r_max, phi_arr(j), f)
                        A_mat(row, new_col) = 1
                        B_vec(row) = f

                ELSE 
                        PRINT *, "You should not be here"

                END IF
                END DO
                END DO

        END SUBROUTINE

        SUBROUTINE GET_LHS_TYPE_I(N, theta, dt, r_arr, phi_arr, k, I0, sigma, A_mat)
                INTEGER, INTENT(IN) :: N, k
                REAL*8, INTENT(IN) :: I0, sigma, theta, dt

                REAL*8, DIMENSION(N), INTENT(IN) :: r_arr, phi_arr
                REAL*8, DIMENSION(N**2, 2*N+1), INTENT(OUT) :: A_mat
                
                INTEGER :: i, j, row, col, new_col
                
                REAL*8 :: dr, dphi, A, B, C, D, f, r_max

                r_max = MAXVAL(r_arr)

                dr = r_arr(2) - r_arr(1)
                dphi = phi_arr(2) - phi_arr(1)

                A_mat = 0.0
                B_vec = 0.0

                DO j=1,N
                DO i=1,N

                CALL get_A(r_arr(i), dr, A)
                CALL get_B(r_arr(i), dr, B)
                CALL get_C(r_arr(i), dr, dphi, C)
                CALL get_D(r_arr(i), dphi, D)

                row = N*(j-1)+i

                IF (i == 1 .AND. j == 1) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = 1 - C * theta * dt

                        col = N*(j-1)+i+1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -(A + B) * theta * dt

                        col = N*j+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -2 * D * theta * dt
                              
                ELSE IF (i == 1 .AND. 1 < j .AND. j < N) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = 1 - C * theta * dt

                        col = N*(j-1)+i+1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -(A + B) * theta * dt

                        col = N*(j-2)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -D * theta * dt

                        col = N*j+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -D * theta * dt

                ELSE IF (1 < i .AND. i < N .AND. j == 1) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = 1 - C * theta * dt

                        col = N*(j-1)+i+1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -A * theta * dt

                        col = N*(j-1)+i-1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -B * theta * dt

                        col = N*j+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -2 * D * theta * dt
                
                ELSE IF (1 < i .AND. i < N .AND. 1 < j .AND. j < N) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = 1 - C * theta * dt

                        col = N*(j-1)+i+1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -A * theta * dt

                        col = N*(j-1)+i-1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -B * theta * dt

                        col = N*(j-2)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -D * theta * dt

                        col = N*j+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -D * theta * dt

                ELSE IF (i == N) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        CALL BC(k, I0, sigma, r_max, phi_arr(j), f)
                        A_mat(row, new_col) = 1

                ELSE IF (j == N) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = 1

                ELSE 
                        PRINT *, "You should not be here"

                END IF

                END DO
                END DO
        END SUBROUTINE

        SUBROUTINE GET_RHS_TYPE_I(N, theta, dt, r_arr, phi_arr, k, I0, sigma, U, U_new)
                INTEGER, INTENT(IN) :: N, k
                REAL*8, INTENT(IN) :: I0, sigma, theta, dt

                REAL*8, DIMENSION(N), INTENT(IN) :: r_arr, phi_arr
                REAL*8, DIMENSION(N**2), INTENT(IN) :: U
                REAL*8, DIMENSION(N**2), INTENT(OUT) :: U_new
                
                INTEGER :: i, j, idx
                REAL*8 :: dr, dphi, A, B, C, D, f, r_max

                dr = r_arr(2) - r_arr(1)
                dphi = phi_arr(2) - phi_arr(1)

                r_max = MAXVAL(r_arr)

                U_new = 0.0

                DO j=1,N
                DO i=1,N

                CALL get_A(r_arr(i), dr, A)
                CALL get_B(r_arr(i), dr, B)
                CALL get_C(r_arr(i), dr, dphi, C)
                CALL get_D(r_arr(i), dphi, D)
                
                idx = N*(j-1)+i

                IF (i == 1 .AND. j == 1) THEN
                        U_new(idx) = (1+(1-theta)*C*dt)*U(idx)
                        U_new(idx) = U_new(idx) + (1-theta)*(A+B)*dt*U(idx+1)
                        U_new(idx) = U_new(idx) + (1-theta)*2*D*dt*U(idx+N)

                ELSE IF (i == 1 .AND. 1 < j .AND. j < N) THEN
                        U_new(idx) = (1+(1-theta)*C*dt)*U(idx) 
                        U_new(idx) = U_new(idx) + (1-theta)*(A+B)*dt*U(idx+1) 
                        U_new(idx) = U_new(idx) + (1-theta)*D*dt*(U(idx+N)+U(idx-N))

                ELSE IF (1 < i .AND. i < N .AND. j == 1) THEN
                        U_new(idx) = (1+(1-theta)*C*dt)*U(idx) 
                        U_new(idx) = U_new(idx) + (1-theta)*A*dt*U(idx+1) 
                        U_new(idx) = U_new(idx) + (1-theta)*B*dt*U(idx-1) 
                        U_new(idx) = U_new(idx) + (1-theta)*2*D*dt*U(idx+N)
                
                ELSE IF (1 < i .AND. i < N .AND. 1 < j .AND. j < N) THEN
                        U_new(idx) = (1+(1-theta)*C*dt)*U(idx) 
                        U_new(idx) = U_new(idx) + (1-theta)*A*dt*U(idx+1) 
                        U_new(idx) = U_new(idx) + (1-theta)*B*dt*U(idx-1) 
                        U_new(idx) = U_new(idx) + (1-theta)*D*dt*U(idx+N)
                        U_new(idx) = U_new(idx) + (1-theta)*D*dt*U(idx-N)

                ELSE IF (j == N) THEN
                        U_new(idx) = 0

                ELSE IF (i == N) THEN
                        CALL BC(k, I0, sigma, r_max, phi_arr(j), f)
                        U_new(idx) = f

                ELSE 
                        PRINT *, "You should not be here"

                END IF

                END DO
                END DO

        END SUBROUTINE

        SUBROUTINE GET_LHS_TYPE_III(N, alpha, theta, dt, r_arr, phi_arr, k, I0, sigma, A_mat)
                INTEGER, INTENT(IN) :: N, k
                REAL*8, INTENT(IN) :: I0, sigma, alpha, theta, dt

                REAL*8, DIMENSION(N), INTENT(IN) :: r_arr, phi_arr
                REAL*8, DIMENSION(N**2, 2*N+1), INTENT(OUT) :: A_mat
                
                INTEGER :: i, j, row, col, new_col
                
                REAL*8 :: dr, dphi, A, B, C, D, E, f, g, r_max

                r_max = MAXVAL(r_arr)

                dr = r_arr(2) - r_arr(1)
                dphi = phi_arr(2) - phi_arr(1)

                A_mat = 0.0
                B_vec = 0.0

                DO j=1,N
                DO i=1,N

                CALL get_A(r_arr(i), dr, A)
                CALL get_B(r_arr(i), dr, B)
                CALL get_C(r_arr(i), dr, dphi, C)
                CALL get_D(r_arr(i), dphi, D)

                row = N*(j-1)+i

                IF (i == 1 .AND. j == 1) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = 1 - C * theta * dt

                        col = N*(j-1)+i+1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -(A + B) * theta * dt

                        col = N*j+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -2 * D * theta * dt
                              
                ELSE IF (i == 1 .AND. 1 < j .AND. j < N) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = 1 - C * theta * dt

                        col = N*(j-1)+i+1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -(A + B) * theta * dt

                        col = N*(j-2)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -D * theta * dt

                        col = N*j+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -D * theta * dt

                ELSE IF (i == 1 .AND. j == N) THEN
                        CALL get_E(r_arr(i), dphi, alpha, E)

                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = 1 - (C + D * E) * theta * dt

                        col = N*(j-1)+i+1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -(A + B) * theta * dt

                        col = N*(j-2)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -2 * D * theta * dt

                ELSE IF (1 < i .AND. i < N .AND. j == 1) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = 1 - C * theta * dt

                        col = N*(j-1)+i+1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -A * theta * dt

                        col = N*(j-1)+i-1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -B * theta * dt

                        col = N*j+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -2 * D * theta * dt
                
                ELSE IF (1 < i .AND. i < N .AND. 1 < j .AND. j < N) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = 1 - C * theta * dt

                        col = N*(j-1)+i+1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -A * theta * dt

                        col = N*(j-1)+i-1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -B * theta * dt

                        col = N*(j-2)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -D * theta * dt

                        col = N*j+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -D * theta * dt

                ELSE IF (1 < N .AND. i < N .AND. j == N) THEN
                        CALL get_E(r_arr(i), dphi, alpha, E)

                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = 1 - (C + D * E) * theta * dt

                        col = N*(j-1)+i+1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -A * theta * dt

                        col = N*(j-1)+i-1
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -B * theta * dt

                        col = N*(j-2)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        A_mat(row, new_col) = -2 * D * theta * dt

                ELSE IF (i == N) THEN
                        col = N*(j-1)+i
                        CALL MODE2_INDEX_MAP(row, col, N, new_col)

                        CALL BC(k, I0, sigma, r_max, phi_arr(j), f)
                        A_mat(row, new_col) = 1

                ELSE 
                        PRINT *, "You should not be here"

                END IF

                END DO
                END DO
        END SUBROUTINE

        SUBROUTINE GET_RHS_TYPE_III(N, alpha, theta, dt, r_arr, phi_arr, k, I0, sigma, U, U_new)
                INTEGER, INTENT(IN) :: N, k
                REAL*8, INTENT(IN) :: I0, sigma, alpha, theta, dt

                REAL*8, DIMENSION(N), INTENT(IN) :: r_arr, phi_arr
                REAL*8, DIMENSION(N**2), INTENT(IN) :: U
                REAL*8, DIMENSION(N**2), INTENT(OUT) :: U_new
                
                INTEGER :: i, j, row, col, new_col
                
                REAL*8 :: dr, dphi, A, B, C, D, E, f, g, r_max

                dr = r_arr(2) - r_arr(1)
                dphi = phi_arr(2) - phi_arr(1)

                r_max = MAXVAL(r_arr)

                U_new = 0.0

                DO j=1,N
                DO i=1,N

                CALL get_A(r_arr(i), dr, A)
                CALL get_B(r_arr(i), dr, B)
                CALL get_C(r_arr(i), dr, dphi, C)
                CALL get_D(r_arr(i), dphi, D)
                
                idx = N*(j-1)+i

                IF (i == 1 .AND. j == 1) THEN
                        U_new(idx) = (1+(1-theta)*C*dt)*U(idx)
                        U_new(idx) = U_new(idx) + (1-theta)*(A+B)*dt*U(idx+1)
                        U_new(idx) = U_new(idx) + (1-theta)*2*D*dt*U(idx+N)

                ELSE IF (i == 1 .AND. 1 < j .AND. j < N) THEN
                        U_new(idx) = (1+(1-theta)*C*dt)*U(idx) 
                        U_new(idx) = U_new(idx) + (1-theta)*(A+B)*dt*U(idx+1) 
                        U_new(idx) = U_new(idx) + (1-theta)*D*dt*(U(idx+N)+U(idx-N))

                ELSE IF (i == 1 .AND. j == N) THEN
                        CALL get_E(r_arr(i), dphi, alpha, E)
                        CALL get_g(r_arr(i), dphi, k, I0, sigma, g)

                        U_new(idx) = (1+(1 - theta)*(C+D*E)*dt)*U(idx)
                        U_new(idx) = U_new(idx) + (1-theta)*(A+B)*dt*U(idx+1) 
                        U_new(idx) = U_new(idx) + (1-theta)*2*D*dt*U(idx-N)
                        U_new(idx) = U_new(idx) + D*dt*g

                ELSE IF (1 < i .AND. i < N .AND. j == 1) THEN
                        U_new(idx) = (1+(1-theta)*C*dt)*U(idx) 
                        U_new(idx) = U_new(idx) + (1-theta)*A*dt*U(idx+1) 
                        U_new(idx) = U_new(idx) + (1-theta)*B*dt*U(idx-1) 
                        U_new(idx) = U_new(idx) + (1-theta)*2*D*dt*U(idx+N)
                
                ELSE IF (1 < i .AND. i < N .AND. 1 < j .AND. j < N) THEN
                        U_new(idx) = (1+(1-theta)*C*dt)*U(idx) 
                        U_new(idx) = U_new(idx) + (1-theta)*A*dt*U(idx+1) 
                        U_new(idx) = U_new(idx) + (1-theta)*B*dt*U(idx-1) 
                        U_new(idx) = U_new(idx) + (1-theta)*D*dt*U(idx+N)
                        U_new(idx) = U_new(idx) + (1-theta)*D*dt*U(idx-N)

                ELSE IF (1 < i.AND. i < N .AND. j == N) THEN
                        CALL get_E(r_arr(i), dphi, alpha, E)
                        CALL get_g(r_arr(i), dphi, k, I0, sigma, g)

                        U_new(idx) = (1+(1 - theta)*(C+D*E)*dt)*U(idx)
                        U_new(idx) = U_new(idx) + (1-theta)*A*dt*U(idx+1) 
                        U_new(idx) = U_new(idx) + (1-theta)*B*dt*U(idx-1) 
                        U_new(idx) = U_new(idx) + (1-theta)*2*D*dt*U(idx-N)
                        U_new(idx) = U_new(idx) + D*dt*g

                ELSE IF (i == N) THEN
                        CALL BC(k, I0, sigma, r_max, phi_arr(j), f)
                        U_new(idx) = f

                ELSE 
                        PRINT *, "You should not be here"

                END IF

                END DO
                END DO

        END SUBROUTINE
