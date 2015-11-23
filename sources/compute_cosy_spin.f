        subroutine compute_cosy_spin(vect,cosy_sp_mat)
*                                                                    *
**********************************************************************
 
* compute the matrix elements of the STM for each particle, using
* the input coordinates of the particle and the STM Taylor series
* expansion coefficients given by COSY:
*
*     input : vect
*             from spin.cmn   spin_ten  == array of coefficients 
*                             iligne_spi(3) == number of coefficients
*                             ip_spi == array of powers for
*                                        target quantities
*     output  in spin.cmn   cosy_sp_mat(3,3)
*             so for particular vect calculate nrow=3,ncol=3 matrix
*            where each element is sum over ii=1,iligne_spi of 
*  spin_ten(ii,ncol,nrow)*(sum of jj=2,5 vect(jj)**ip_spi(ii,jj,nrow))
c
        implicit none
c
        double precision terme, sousterme, vect(5)
        double precision cosy_sp_mat(3,3)
c
        integer ii,jj
        integer il, ipower, ncolstm
c
        INCLUDE 'spin.cmn'

* Initialisation:
        do ii=1,3
           do jj=1,3
              cosy_sp_mat(ii,jj)=0.0
        enddo
        enddo
        ncolstm = 5
* First and third row of the spin transfer matrix:
* based on the fact that they have the monomial powers
* (but with different coefficients)
        do il = 1, iligne_spi(1) 
            terme = 1.d0
            do ipower = 1, ncolstm ! x_tgt=0 thus gives no contribution
              if (ip_spi(il,ipower,1).eq.0) then
                sousterme = 1.d0
              else
                sousterme = 
     s                  vect(ipower)**ip_spi(il,ipower,1)
              endif
              terme = terme * sousterme
            enddo  ! ipower
            cosy_sp_mat(1,1) = cosy_sp_mat(1,1) + sp_ten(il,1,1)*terme
            cosy_sp_mat(2,1) = cosy_sp_mat(2,1) + sp_ten(il,2,1)*terme
            cosy_sp_mat(3,1) = cosy_sp_mat(3,1) + sp_ten(il,3,1)*terme
            cosy_sp_mat(1,3) = cosy_sp_mat(1,3) + sp_ten(il,1,3)*terme
            cosy_sp_mat(2,3) = cosy_sp_mat(2,3) + sp_ten(il,2,3)*terme
            cosy_sp_mat(3,3) = cosy_sp_mat(3,3) + sp_ten(il,3,3)*terme
        enddo    ! il
* Second row of the STM:
        do il = 1, iligne_spi(2) ! 
          if (ip_spi(il,1,2).eq.0) then
            terme = 1.d0
            do ipower = 2, ncolstm ! x_tgt=0 thus gives no contribution
              if (ip_spi(il,ipower,2).eq.0) then
                sousterme = 1.d0
              else
                sousterme = 
     s                  vect(ipower)**ip_spi(il,ipower,2)
              endif
              terme = terme * sousterme
            enddo  ! ipower
            cosy_sp_mat(1,2) = cosy_sp_mat(1,2) + sp_ten(il,1,2)*terme
            cosy_sp_mat(2,2) = cosy_sp_mat(2,2) + sp_ten(il,2,2)*terme
            cosy_sp_mat(3,2) = cosy_sp_mat(3,2) + sp_ten(il,3,2)*terme
          endif  ! power of x = 0 otherwise contribution = 0
        enddo    ! il
        return
        end
