! Harmonic-space Internal Linear Combination for full-sky HEALPix maps

program hilc

use healpix_modules

implicit none

!======================================================================================
real(DP), allocatable     :: map_in(:,:), map_out(:,:), freqmaps(:,:), cls(:,:,:,:), cov_inv(:,:), w_l(:)
complex(DPC), allocatable :: alms(:,:,:,:), alm_out(:,:,:)
integer                   :: nfreqmaps, nside, nside_in, npix_in, nmaps, ord_in, i, j, l, m, lmax
character(len=80)         :: filename, header(43)
!======================================================================================

! Read number of frequency maps
write(*,'(/,X,A)',advance='no') "Enter number of frequency maps: "
read(*,*) nfreqmaps

! Read common resolution parameter for the frequency maps
write(*,'(/,X,A)',advance='no') "Enter resulution parameter (Nside): "
read(*,*) nside

! Maximum multipole moment
lmax = 3*nside - 1

! Provide storage for frequency maps array
allocate(freqmaps(0:12*nside**2-1,nfreqmaps))

do i = 1, nfreqmaps
	! Read input file name of frequency map
	write(*,'(/,X,"File name of map ", I0,": ")',advance='no') i
	read(*,*) filename
	
	! Store frequency map in array "map_in"
	npix_in = getsize_fits(filename, nmaps=nmaps, nside=nside_in, ordering=ord_in)
	allocate(map_in(0:12*nside_in**2-1,1))
	call input_map(filename, map_in, npix_in, 1)
	
	! Change ordering to RING if necessary
	if (ord_in == 2) call convert_nest2ring(nside_in, map_in(:,1))
	
	! Change frequency maps resolution to the given resolution parameter
	if (nside_in == nside) freqmaps(:,i) = map_in(:,1)
	if (nside_in /= nside) call udgrade_ring(map_in(:,1), nside_in, freqmaps(:,i), nside)
	
	deallocate(map_in)
end do


write(*,'(/,X,A)') "Computing ILC weights"

! Provide storage for spherical harmonic coefficients
allocate(alms(nfreqmaps,1,0:lmax,0:lmax))

! Get spherical harmonic coefficients
!$OMP PARALLEL PRIVATE(i) SHARED(nfreqmaps, nside, lmax, freqmaps, alms)
!$OMP DO SCHEDULE(DYNAMIC)
do i = 1, nfreqmaps
	call map2alm(nside, lmax, lmax, freqmaps(:,i), alms(i,:,:,:))
end do
!$OMP END DO
!$OMP END PARALLEL

! Provide storage for the (cross/auto-)correlation functions "Cl"
allocate(cls(nfreqmaps,nfreqmaps,0:lmax,1))

! Compute (cross/auto-)correlation functions "Cl"
!$OMP PARALLEL PRIVATE(i, j) SHARED(nfreqmaps, lmax, alms, cls)
!$OMP DO SCHEDULE(DYNAMIC)
do j = 1, nfreqmaps
	do i = 1, nfreqmaps
		if (i /= j) call alm2cl(lmax, lmax, alms(i,:,:,:), alms(j,:,:,:), cls(i,j,:,:))
		if (i == j) call alm2cl(lmax, lmax, alms(i,:,:,:), cls(i,i,:,:))
	end do
end do
!$OMP END DO
!$OMP END PARALLEL

! Provide storage for inverse covariance matrix, weight vector, output spherical harmonic coefficients and output map
allocate(cov_inv(nfreqmaps,nfreqmaps), w_l(nfreqmaps), alm_out(1,0:lmax,0:lmax), map_out(0:12*nside**2-1,1))

!$OMP PARALLEL PRIVATE(l, cov_inv, w_l, m) SHARED(lmax, nfreqmaps, cls, alms)
!$OMP DO SCHEDULE(DYNAMIC)
do l = 0, lmax
	! Compute inverse covariance matrix
	call mtx_inv(nfreqmaps, cls(:,:,l,1), cov_inv)
	
	! Compute the ILC weights given the covariance matrix of the frequency maps
	w_l = sum(cov_inv, dim=2) / sum(cov_inv)

	! Compute CMB spherical harmonic coefficients
	do m = 0, l
		alm_out(1,l,m) = dot_product(w_l, alms(:,1,l,m))
	end do
end do
!$OMP END DO
!$OMP END PARALLEL

! Compute output map
call alm2map(nside, lmax, lmax, alm_out, map_out(:,1))

! Generate output file
write(*,'(/,X,A)',advance='no') "Enter output file name: "
read(*,*) filename
call write_minimal_header(header, 'map', nside=nside, order=1)
call output_map(map_out, header, filename)

! Free memory from allocatable objects
deallocate(alms, map_out, freqmaps, cov_inv, w_l, alm_out)

contains

! Compute the covariance matrix of "nvar" variables with "nobs" observations stacked in a matrix "A"
subroutine cov_mtx(nvar, nobs, A, C)
	integer, intent(in)     :: nvar, nobs
	real(DP), intent(inout) :: A(nvar,nobs)
	real(DP), intent(out)   :: C(nvar,nvar)
	integer                 :: i, j
	real(DP)                :: mean_i, DeltaA(nvar,nobs)
	
	! Substract the mean of each variable from each observation
	do i = 1, nvar
		mean_i = sum(A(i,:)) / nobs
		DeltaA(i,:) = A(i,:) - mean_i
	end do
	
	! Compute the covariance matrix
	do i = 1, nvar
		do j = 1, nvar
			if (i <= j) C(i,j) = dot_product(DeltaA(i,:), DeltaA(j,:)) / nobs
			if (i >  j) C(i,j) = C(j,i)
		end do
	end do
end subroutine cov_mtx

! Solve the real system of "n" symmetric linear equations in "n" unknowns in the form "A*X=B"
subroutine mtx_inv(n, A, A_inv)
	integer, intent(in)   :: n
	real(DP), intent(in)  :: A(n,n)
	real(DP), intent(out) :: A_inv(n,n)
	real(DP), allocatable :: work(:)
	real(DP)              :: get_lwork(1)
	integer               :: ipiv(n), info, lwork, i, j
	
	A_inv = A ! Initialize inverse matrix
	info  = 0 ! Status indicator ("stat/=0" indicates an error)
	
	! Calculate the optimal size of the WORK array as the first entry of the GET_LWORK array
	call dsytrf('U', n, A_inv, n, ipiv, get_lwork, -1, info)

	! Provide storage for the WORK array
	lwork = get_lwork(1)
	allocate(work(lwork))
	
	! LAPACK subroutine to compute the factorization of a real symmetric matrix A using the Bunch-Kaufman diagonal pivoting method
	call dsytrf('U', n, A_inv, n, ipiv, work, lwork, info)

	! Stop the program if necessary
	if (info /= 0) call fatal_error('Singular matrix in subroutine MTX_INV()')
	
	! Free memory from Work
	deallocate(work)
	
	! Resizing WORK array
	allocate(work(n))
	
	! LAPACK subroutine to get "A^(-1)" for a symmetric matrix "A"
	call dsytri('U', n, A_inv, n, ipiv, work, info)
	
	! Stop the program if necessary
	if (info /= 0) call fatal_error('Singular matrix in subroutine MTX_INV()')
	
	! Off-diagonal elements of the matrix in its lower triangular part
	do i = 2, n
		do j = 1, i-1
			A_inv(i,j) = A_inv(j,i)
		end do
	end do
end subroutine mtx_inv

end program hilc
