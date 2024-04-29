        SUBROUTINE MODE2_INDEX_MAP(row, col, width, new_col)
                INTEGER, INTENT(IN) :: row, col, width
                INTEGER, INTENT(OUT) :: new_col

                new_col = (width + 1) + (col-row)
        END SUBROUTINE

        SUBROUTINE LINSPACE(x_min, x_max, num_points, arr)
                REAL*8, INTENT(IN) :: x_min, x_max
                INTEGER, INTENT(IN) :: num_points
                REAL*8, DIMENSION(num_points), INTENT(OUT) :: arr
                INTEGER :: i

                DO i=1,num_points
                        arr(i) = (i-1) * (x_max-x_min)/(num_points-1)
                        arr(i) = arr(i) + x_min
                END DO
        END SUBROUTINE

        SUBROUTINE POLAR_2_CART(r, theta, x, y)
                REAL*8, INTENT(IN) :: r, theta
                REAL*8, INTENT(OUT) :: x, y

                x = r * COS(theta)
                y = r * SIN(theta)
        END SUBROUTINE 
