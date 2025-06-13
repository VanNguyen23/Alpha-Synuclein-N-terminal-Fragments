!  bsheet_content.f90 
!
!  FUNCTIONS:
!  Console1 - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: bsheet_content
!
!  Author:  Haoyu Wang.
!
!****************************************************************************

    module my_functions
    implicit none

    contains
    
    subroutine read_input(input_file, pep1, pep2, nc1, nc2, nb1, nb2)
    character(len=*) :: input_file
    character(len=100) :: pep1, pep2
    integer, intent(out) :: nc1, nc2, nb1, nb2

    character(len=200) :: line
    character(len=100) :: key, value, substr
    character(len=1) :: char
    integer :: i, ios, g_count, eqpos

    character(len=100), dimension(8) :: lines

    ! Read first 8 lines
    open(unit=10, file=input_file, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'Error opening file:', trim(input_file)
        stop
    end if

    do i = 1, 8
        read(10,'(A)', iostat=ios) line
        if (ios /= 0) exit
        lines(i) = adjustl(trim(line))
    end do
    close(10)

    ! read in pep1
    if (index(lines(2), '=') > 0) then
        call parse_line(lines(2), pep1)
    end if
    ! read in pep2
    if (index(lines(6), '=') > 0) then
        call parse_line(lines(6), pep2)
    end if

    ! read number of chain (nc1) for pep1
    eqpos = index(lines(4), '=')
    if (eqpos > 0) then
        substr = adjustl(trim(lines(4)(eqpos + 1:)))
        read(substr, *) nc1
    end if
    ! read number of chain (nc2) for pep2
    eqpos = index(lines(8), '=')
    if (eqpos > 0) then
        substr = adjustl(trim(lines(8)(eqpos + 1:)))
        read(substr, *) nc2
    end if

    ! Count beads: non-"G" = 4, "G" = 3
    g_count = 0
    do i = 1, len_trim(pep1)
        char = pep1(i:i)
        if (char == 'G') g_count = g_count + 1
    end do
    nb1 = 4 * len_trim(pep1) - g_count

    g_count = 0
    do i = 1, len_trim(pep2)
        char = pep2(i:i)
        if (char == 'G') g_count = g_count + 1
    end do
    nb2 = 4 * len_trim(pep2) - g_count

    end subroutine read_input

    subroutine parse_line(line, outputs)
        implicit none
        character(len=*), intent(in) :: line
        character(len=*), intent(out) :: outputs
        integer :: eq_pos

        eq_pos = index(line, '=')
        if (eq_pos > 0) then
            outputs = adjustl(trim(line(eq_pos + 1:)))
        else
            outputs = ''
        end if
    end subroutine parse_line

    
    subroutine read_energy(energy_file, ttotal, i, anneal_steps, t_series, npts)
    implicit none
    character(len=*)                  :: energy_file
    real(8)                           :: ttotal, timetemp, factor
    real(8)                           :: dummy0, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7, dummy8, dummy9
    integer                           :: i, anneal_steps, npts, unit, ios, count, j
    real(8), allocatable              :: t_series(:), t(:)
    character(len=100)                :: line, temp_file

    ! 1. Determine temperature
    if (i == 0) then
        timetemp = 0.5
    else if (i <= anneal_steps) then
        write(temp_file, '(A,I0)') 'inputs/annealtemp_', i
        open(10, file=trim(temp_file), status='old', action='read')
        read(10,*) timetemp
        close(10)
    else
        temp_file = 'inputs/simtemp'
        open(10, file=trim(temp_file), status='old', action='read')
        read(10,*) timetemp
        close(10)
    end if

    ! 2. Read energy file, we only read the reduce time
    open(10, file=trim(energy_file), status='old', action='read')
    allocate(t(10000))  ! Max lines (adjust as needed)
    count = 0
    
    do
        read(10, *, iostat=ios) dummy0, t(count+1), dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7, dummy8, dummy9
        if (ios /= 0) exit
        count = count + 1
    end do
    close(10)

    ! count is length of time array, npts is length the the array remove the first one (which is t=0)
    if (count <= 1) then
        allocate(t_series(0))
        npts = 0
        return
    end if

    allocate(t_series(count - 1))
    do j = 2, count
        t_series(j - 1) = ttotal + t(j) * 0.96d0 * 0.001d0 * 3.3d0 / sqrt(timetemp * 12.0d0)
    end do
    ttotal = t_series(count - 1)
    npts = count - 1

    end subroutine read_energy
    
    subroutine read_bptnr_multi(filename, nc1, nc2, nb1, nb2, nframes, collision_nums, bptnrs, hbmats)
    implicit none
    character(len=*) :: filename
    integer                                 :: nc1, nc2, nb1, nb2, nframes
    integer                                 :: nop1, nop2, noptotal, nchain
    integer                                 :: unit, ios, i, j, frame    
    integer                                 :: bead, partner, chain_i, chain_j
    
    integer, allocatable                    :: bptnrs(:,:,:)
    integer, allocatable                    :: hbmats(:,:,:)
    integer, allocatable                    :: hb(:)    
    integer(8), allocatable                 :: collision_nums(:)
    integer(8)                              :: coll
    
    ! Initialize
    nop1 = nc1 * nb1
    nop2 = nc2 * nb2
    noptotal = nop1 + nop2
    nchain = nc1 + nc2

    allocate(hb(noptotal))
    allocate(collision_nums(0))
    allocate(bptnrs(noptotal, 2, 0))
    allocate(hbmats(nchain, nchain, 0))

    frame = 0
    open(unit=15, file=filename, form='unformatted', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'Error opening file: ', trim(filename)
        stop
    end if

    do
        read(15, iostat=ios) coll, hb
        if (ios /= 0) exit

        frame = frame + 1

        ! Expand arrays
        call expand_arrays(frame, noptotal, nchain, collision_nums, bptnrs, hbmats)

        collision_nums(frame) = coll

        do i = 1, noptotal
            bptnrs(i, 1, frame) = i
            bptnrs(i, 2, frame) = hb(i)
        end do

        hbmats(:, :, frame) = 0
        do i = 1, noptotal
            partner = hb(i)
            if (partner /= 0) then
                if (i <= nop1) then
                    chain_i = (i - 1) / nb1 + 1
                else
                    chain_i = (i - nop1 - 1) / nb2 + nc1 + 1
                end if

                if (partner <= nop1) then
                    chain_j = (partner - 1) / nb1 + 1
                else
                    chain_j = (partner - nop1 - 1) / nb2 + nc1 + 1
                end if

                hbmats(chain_i, chain_j, frame) = hbmats(chain_i, chain_j, frame) + 1
            end if
        end do
    end do
    close(15)

    nframes = frame
    
    deallocate(hb)
    end subroutine read_bptnr_multi
    
    
    subroutine expand_arrays(frame, noptotal, nchain, collision_nums, bptnrs, hbmats)
    implicit none
    integer, intent(in) :: frame, noptotal, nchain
    integer(8), allocatable, intent(inout) :: collision_nums(:)
    integer, allocatable, intent(inout) :: bptnrs(:,:,:)
    integer, allocatable, intent(inout) :: hbmats(:,:,:)

    integer(8), allocatable :: temp_col(:)
    integer, allocatable :: temp_bptnr(:,:,:)
    integer, allocatable :: temp_hbm(:,:,:)

    allocate(temp_col(frame))
    if (frame > 1) temp_col(1:frame-1) = collision_nums
    deallocate(collision_nums)
    collision_nums = temp_col

    allocate(temp_bptnr(noptotal, 2, frame))
    if (frame > 1) temp_bptnr(:,:,1:frame-1) = bptnrs
    deallocate(bptnrs)
    bptnrs = temp_bptnr

    allocate(temp_hbm(nchain, nchain, frame))
    if (frame > 1) temp_hbm(:,:,1:frame-1) = hbmats
    deallocate(hbmats)
    hbmats = temp_hbm
    end subroutine expand_arrays

    
    
    subroutine bsheet_content(bptnr, hbmatrix, nc1, nc2, nb1, nb2, nres1, nres2, anti, para, turn, undetermined, content)
    implicit none
    ! Inputs
    integer, intent(in) :: nc1, nc2, nb1, nb2, nres1, nres2
    integer, intent(in) :: bptnr(:,:), hbmatrix(:,:)

    ! Outputs
    integer, intent(out) :: anti, para, turn, undetermined
    real(8), intent(out) :: content

    ! Locals
    integer :: nchain, total_residues
    integer, allocatable                    :: bsheet_assign(:)
    integer                                 :: i, j, m, ii, jj, jj2, ii2
    integer                                 :: nbi, nbii, nresi, nbj, nbjj, nresj
    integer, allocatable                    :: nh_i(:), nh_j(:), co_i(:), co_j(:)  
    integer, allocatable                    :: resNHi(:), resCOi(:), resNHj(:), resCOj(:) 
    integer                                 :: nresii, nresjj, nbond_ij, nbond_ji
    integer                                 :: NH_i_pos, CO_i_pos, CO_j_pos, NH_j_pos
    integer                                 :: k, x, counter, diff, nh_basg, co_basg
    real(8)                                 :: a, b
    logical                                 :: found

    anti = 0; para = 0; turn = 0; undetermined = 0
    nchain = nc1 + nc2

    total_residues = nc1 * nres1 + nc2 * nres2
    allocate(bsheet_assign(total_residues))
    bsheet_assign = 0

    do i = 1, nchain
        do j = i, nchain
            nbond_ij = 0; nbond_ji = 0                ! nbond_i: number of bonds (NH_i to CO_j), and same ......
            if (hbmatrix(i, j) > 0) then
                ! Determine nb and nres for chains i and j
                if (i <= nc1) then                      ! if current i is pep1
                    nbi = nb1                           ! number of beads in pep1
                    nresi = nres1                       ! number of residues in pep1
                    nbii = (i - 1) * nb1 + 1            ! chain i bead starting index
                    nresii = (i - 1) * nres1            ! chain i residue starting index
                else
                    nbi = nb2
                    nresi = nres2
                    nbii = (i - nc1 - 1) * nb2 + nc1 * nb1  + 1
                    nresii = (i - nc1 - 1) * nres2 + nc1 * nres1
                endif

                if (j <= nc1) then
                    nbj = nb1
                    nresj = nres1
                    nbjj = (j - 1) * nb1 + 1
                    nresjj = (j - 1) * nres1
                else
                    nbj = nb2
                    nresj = nres2
                    nbjj = (j - nc1 - 1) * nb2 + nb1 * nc1 + 1
                    nresjj = (j - nc1 - 1) * nres2 + nres1 * nc1
                endif

                
                allocate(nh_i(nresi))
                allocate(nh_j(nresj))
                allocate(co_i(nresi))
                allocate(co_j(nresj))
                nh_i = 0; co_j = 0; nh_j = 0; co_i = 0

                do m = 1, nresi
                    CO_j_pos = bptnr(nbii + nresi + m - 1, 2)
                    
                    jj  = int(real(CO_j_pos/ nb1 + 1))
                    jj2 = int(real(CO_j_pos  - nc1*nb1)/ nb2 + 1 + nc1) 
                    
                    if (( jj == j .or. jj2 == j) .and. (CO_j_pos /=0 )) then
                        nh_i(m) = nbii + nresi + m - 1             ! index of m_th NH bead in chain i
                        co_j(m) = CO_j_pos                         ! index of m_th CO bead in chain j
                        nbond_ij = nbond_ij + 1
                    endif
                enddo

                     
                do m = 1, nresj
                    CO_i_pos = bptnr(nbjj + nresj + m - 1, 2)
                    ii  = int(real(CO_i_pos/ nb1 + 1))
                    ii2 = int(real(CO_i_pos  - nc1*nb1)/ nb2 + 1 + nc1)    
                    if (( ii == i .or. ii2 == i) .and. (CO_i_pos /=0 )) then
                        nh_j(m) = nbjj + nresj + m - 1
                        co_i(m) = CO_i_pos
                        nbond_ji = nbond_ji + 1

                        
                    endif
                enddo

                nh_i = pack(nh_i, nh_i /= 0)
                nh_j = pack(nh_j, nh_j /= 0)
                co_i = pack(co_i, co_i /= 0)
                co_j = pack(co_j, co_j /= 0)
                
                
                ! Determine sheet direction from CO changes
                a = 0.0
                b = 0.0
                if (nbond_ij >= 2) then
                    do m = 2, nbond_ij
                        if (co_j(m) /= 0 .and. co_j(m - 1) /= 0) then
                            a = a + (co_j(m) - co_j(m - 1))
                        endif
                    enddo
                endif

                if (nbond_ji >= 2) then
                    do m = 2, nbond_ji
                        if (co_i(m) /= 0 .and. co_i(m - 1) /= 0) then
                            b = b + (co_i(m) - co_i(m - 1))
                        endif
                    enddo
                endif

                if (i <= nc1) then                  
                    resNHi = nh_i - nresi - (i - 1) * nbi                   ! calculate NH is which residue in chain i, array calculation
                    resCOi = co_i - 2*nresi - (i - 1) * nbi                 ! calculate CO is which residue in chain i
                else
                    resNHi = nh_i - nresi - (i - nc1 - 1) * nbi - nc1 * nb1
                    resCOi = co_i - 2*nresi - (i - nc1 - 1) * nbi - nc1 * nb1
                endif

                if (j <= nc1) then
                    resCOj = co_j - 2*nresj - (j - 1) * nbj
                    resNHj = nh_j - nresj - (j - 1) * nbj
                else
                    resCOj = co_j - 2*nresj - (j - nc1 - 1) * nbj - nc1 * nb1
                    resNHj = nh_j - nresj - (j - nc1 - 1) * nbj - nc1 * nb1
                endif
                
!!!!!!!!!!!!!!!!!!!!! LOOP 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                do m = 1, size(resNHi)
                    NH_i_pos = resNHi(m)
                    CO_j_pos = resCOj(m)

                    ! Check if CO_j_pos exists in resNHj
                    found = .false.
                    do k = 1, size(resNHj)
                        if (resNHj(k) == CO_j_pos) then
                            found = .true.
                            CO_i_pos = resCOi(k)

                            if (CO_i_pos == NH_i_pos) then                      ! antiparallel
                                if (bsheet_assign(nresii + NH_i_pos) <= 1) then
                                    bsheet_assign(nresii + NH_i_pos) = 2
                                else if (bsheet_assign(nresii + NH_i_pos) == 3) then
                                    bsheet_assign(nresii + NH_i_pos) = 5
                                endif

                                if (bsheet_assign(nresjj + CO_j_pos) <= 1) then
                                    bsheet_assign(nresjj + CO_j_pos) = 2
                                else if (bsheet_assign(nresjj + CO_j_pos) == 3) then
                                    bsheet_assign(nresjj + CO_j_pos) = 5
                                endif

                            else if (CO_i_pos == NH_i_pos - 2) then             ! parallel
                                if (bsheet_assign(nresii + NH_i_pos) <= 1) then
                                    bsheet_assign(nresii + NH_i_pos) = 3
                                else if (bsheet_assign(nresii + NH_i_pos) == 2) then
                                    bsheet_assign(nresii + NH_i_pos) = 5
                                endif

                                if(NH_i_pos - 1 > 0) then
                                    if (bsheet_assign(nresii + NH_i_pos - 1) <= 1) then
                                        bsheet_assign(nresii + NH_i_pos - 1) = 3
                                    else if (bsheet_assign(nresii + NH_i_pos - 1) == 2) then
                                        bsheet_assign(nresii + NH_i_pos - 1) = 5
                                    endif
                                endif

                                if (bsheet_assign(nresjj + CO_j_pos) <= 1) then
                                    bsheet_assign(nresjj + CO_j_pos) = 3
                                else if (bsheet_assign(nresjj + CO_j_pos) == 2) then
                                    bsheet_assign(nresjj + CO_j_pos) = 5
                                endif

                                if (bsheet_assign(nresii + CO_i_pos) <= 1) then
                                    bsheet_assign(nresii + CO_i_pos) = 3
                                else if (bsheet_assign(nresii + CO_i_pos) == 2) then
                                    bsheet_assign(nresii + CO_i_pos) = 5
                                end if
                            end if
                            
                            exit
                        end if
                    end do
                
                    if (.not. found) then
                        ! CO_j_pos + 2 in resNHj
                        found = .false.
                        do k = 1, size(resNHj)
                            if (resNHj(k) == CO_j_pos + 2) then
                                found = .true.
                                CO_i_pos = resCOi(k)

                                if (CO_i_pos == NH_i_pos) then             ! parallel
                                    if (bsheet_assign(nresii + NH_i_pos) <= 1) then
                                        bsheet_assign(nresii + NH_i_pos) = 3
                                    else if (bsheet_assign(nresii + NH_i_pos) == 2) then
                                        bsheet_assign(nresii + NH_i_pos) = 5
                                    end if
                
                                    if (bsheet_assign(nresjj + CO_j_pos) <= 1) then
                                        bsheet_assign(nresjj + CO_j_pos) = 3
                                    else if (bsheet_assign(nresjj + CO_j_pos) == 2) then
                                        bsheet_assign(nresjj + CO_j_pos) = 5
                                    end if
                
                                    if (bsheet_assign(nresjj + CO_j_pos + 1) <= 1) then
                                        bsheet_assign(nresjj + CO_j_pos + 1) = 3
                                    else if (bsheet_assign(nresjj + CO_j_pos + 1) == 2) then
                                        bsheet_assign(nresjj + CO_j_pos + 1) = 5
                                    end if

                                    if (bsheet_assign(nresjj + CO_j_pos + 2) <= 1) then
                                        bsheet_assign(nresjj + CO_j_pos + 2) = 3
                                    else if (bsheet_assign(nresjj + CO_j_pos + 2) == 2) then
                                        bsheet_assign(nresjj + CO_j_pos + 2) = 5
                                    end if
                                end if
                                exit
                            end if
                        end do

                        if (.not. found) then
                            ! no pair — single H-bond
                            if (bsheet_assign(nresii + NH_i_pos) <= 1) then
                                bsheet_assign(nresii + NH_i_pos) = 1
                            end if
                            
                            if (bsheet_assign(nresjj + CO_j_pos) <= 1) then
                                bsheet_assign(nresjj + CO_j_pos) = 1
                            end if
                        end if
                    end if
                            
                    found = .false.
                    do k = 1, size(resCOi)
                        if (resCOi(k) == NH_i_pos - 2) then
                            found = .true.
                            NH_j_pos = resCOi(k)  ! same as CO_i_p2_pos
                            exit
                        end if
                    end do

                    if (found) then
                        if (NH_j_pos == CO_j_pos + 2) then
                            if (bsheet_assign(nresii + NH_i_pos) <= 1) then
                                bsheet_assign(nresii + NH_i_pos) = 2
                            else if (bsheet_assign(nresii + NH_i_pos) == 3) then
                                bsheet_assign(nresii + NH_i_pos) = 5
                            end if

                            if(NH_i_pos - 1 > 0) then
                                if (bsheet_assign(nresii + NH_i_pos - 1) <= 1) then
                                    bsheet_assign(nresii + NH_i_pos - 1) = 2
                                else if (bsheet_assign(nresii + NH_i_pos - 1) == 3) then
                                    bsheet_assign(nresii + NH_i_pos - 1) = 5
                                endif
                            endif

                            if (NH_i_pos - 2 > 0) then
                                if (bsheet_assign(nresii + NH_i_pos - 2) <= 1) then
                                    bsheet_assign(nresii + NH_i_pos - 2) = 2
                                else if (bsheet_assign(nresii + NH_i_pos - 2) == 3) then
                                    bsheet_assign(nresii + NH_i_pos - 2) = 5
                                endif
                            endif

                            if (bsheet_assign(nresjj + CO_j_pos) <= 1) then
                                bsheet_assign(nresjj + CO_j_pos) = 2
                            else if (bsheet_assign(nresjj + CO_j_pos) == 3) then
                                bsheet_assign(nresjj + CO_j_pos) = 5
                            endif

                            if (bsheet_assign(nresjj + CO_j_pos + 1) <= 1) then
                                bsheet_assign(nresjj + CO_j_pos + 1) = 2
                            else if (bsheet_assign(nresjj + CO_j_pos + 1) == 3) then
                                bsheet_assign(nresjj + CO_j_pos + 1) = 5
                            end if

                            if (bsheet_assign(nresjj + CO_j_pos + 2) <= 1) then
                                bsheet_assign(nresjj + CO_j_pos + 2) = 2
                            else if (bsheet_assign(nresjj + CO_j_pos + 2) == 3) then
                                bsheet_assign(nresjj + CO_j_pos + 2) = 5
                            end if

                        end if
                    end if     
                end do
!!!!!!!!!!!!!!!!!!!!! LOOP 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                
                do m = 1, size(resNHj)
                    NH_j_pos = resNHj(m)
                    CO_i_pos = resCOi(m)

                    ! Check if CO_j_pos exists in resNHj
                    found = .false.
                    do k = 1, size(resNHi)
                        if (resNHi(k) == CO_i_pos) then
                            found = .true.
                            CO_j_pos = resCOj(k)

                            if (CO_j_pos == NH_j_pos) then                      ! antiparallel
                                if (bsheet_assign(nresjj + NH_j_pos) <= 1) then
                                    bsheet_assign(nresjj + NH_j_pos) = 2
                                else if (bsheet_assign(nresjj + NH_j_pos) == 3) then
                                    bsheet_assign(nresjj + NH_j_pos) = 5
                                end if

                                if (bsheet_assign(nresii + CO_i_pos) <= 1) then
                                    bsheet_assign(nresii + CO_i_pos) = 2
                                else if (bsheet_assign(nresii + CO_i_pos) == 3) then
                                    bsheet_assign(nresii + CO_i_pos) = 5
                                end if
!!!!
                            else if (CO_i_pos == NH_i_pos - 2) then             ! parallel
                                if (bsheet_assign(nresjj + NH_j_pos) <= 1) then
                                    bsheet_assign(nresjj + NH_j_pos) = 3
                                else if (bsheet_assign(nresjj + NH_j_pos) == 2) then
                                    bsheet_assign(nresjj + NH_j_pos) = 5
                                end if

                                if (NH_j_pos - 1 > 0) then
                                    if (bsheet_assign(nresjj + NH_j_pos - 1) <= 1) then
                                        bsheet_assign(nresjj + NH_j_pos - 1) = 3
                                    else if (bsheet_assign(nresjj + NH_j_pos - 1) == 2 ) then
                                        bsheet_assign(nresjj + NH_j_pos - 1) = 5
                                    end if
                                endif

                                if (bsheet_assign(nresii + CO_i_pos) <= 1) then
                                    bsheet_assign(nresii + CO_i_pos) = 3
                                else if (bsheet_assign(nresii + CO_i_pos) == 2) then
                                    bsheet_assign(nresii + CO_i_pos) = 5
                                end if

                                if (bsheet_assign(nresjj + CO_j_pos) <= 1) then
                                    bsheet_assign(nresjj + CO_j_pos) = 3
                                else if (bsheet_assign(nresjj + CO_j_pos) == 2) then
                                    bsheet_assign(nresjj + CO_j_pos) = 5
                                end if
                            end if
                            
                            exit
                        end if
                    end do
!!!!!                
                    if (.not. found) then
                        found = .false.
                        do k = 1, size(resNHi)
                            if (resNHi(k) == CO_i_pos + 2) then
                                found = .true.
                                CO_j_pos = resCOj(k)

                                if (CO_j_pos == NH_j_pos) then             ! parallel
                                    if (bsheet_assign(nresjj + NH_j_pos) <= 1) then
                                        bsheet_assign(nresjj + NH_j_pos) = 3
                                    else if (bsheet_assign(nresjj + NH_j_pos) == 2) then
                                        bsheet_assign(nresjj + NH_j_pos) = 5
                                    end if
                
                                    if (bsheet_assign(nresii + CO_i_pos) <= 1) then
                                        bsheet_assign(nresii + CO_i_pos) = 3
                                    else if (bsheet_assign(nresii + CO_i_pos) == 2) then
                                        bsheet_assign(nresii + CO_i_pos) = 5
                                    end if
                
                                    if (bsheet_assign(nresii + CO_i_pos + 1) <= 1) then
                                        bsheet_assign(nresii + CO_i_pos + 1) = 3
                                    else if (bsheet_assign(nresii + CO_i_pos + 1) == 2) then
                                        bsheet_assign(nresii + CO_i_pos + 1) = 5
                                    end if

                                    if (bsheet_assign(nresii + CO_i_pos + 2) <= 1) then
                                        bsheet_assign(nresii + CO_i_pos + 2) = 3
                                    else if (bsheet_assign(nresii + CO_i_pos + 2) == 2) then
                                        bsheet_assign(nresii + CO_i_pos + 2) = 5
                                    end if
                                end if
                                exit
                            end if
                        end do

                        if (.not. found) then
                            ! no pair — single H-bond
                            if (bsheet_assign(nresjj + NH_j_pos) <= 1) then
                                bsheet_assign(nresjj + NH_j_pos) = 1
                            end if
                            
                            if (bsheet_assign(nresii + CO_i_pos) <= 1) then
                                bsheet_assign(nresii + CO_i_pos) = 1
                            end if
                        end if
                    end if
                            
                    found = .false.
                    do k = 1, size(resCOj)
                        if (resCOj(k) == NH_j_pos - 2) then
                            found = .true.
                            NH_i_pos = resCOj(k)  ! same as CO_i_p2_pos
                            exit
                        end if
                    end do

                    if (found) then
                        if (NH_i_pos == CO_i_pos + 2) then
                            if (bsheet_assign(nresjj + NH_j_pos) <= 1) then
                                bsheet_assign(nresjj + NH_j_pos) = 2
                            else if (bsheet_assign(nresjj + NH_j_pos) == 3) then
                                bsheet_assign(nresjj + NH_j_pos) = 5
                            end if
                            
                            if(NH_j_pos - 1 > 0) then
                                if (bsheet_assign(nresjj + NH_j_pos - 1) <= 1) then
                                    bsheet_assign(nresjj + NH_j_pos - 1) = 2
                                else if (bsheet_assign(nresjj + NH_j_pos - 1) == 3) then
                                    bsheet_assign(nresjj + NH_j_pos - 1) = 5
                                endif
                            endif
                            
                            if (NH_j_pos - 2 > 0) then
                                if (bsheet_assign(nresjj + NH_j_pos - 2) <= 1) then
                                    bsheet_assign(nresjj + NH_j_pos - 2) = 2
                                else if (bsheet_assign(nresjj + NH_j_pos - 2) == 3) then
                                    bsheet_assign(nresjj + NH_j_pos - 2) = 5
                                endif
                            endif
                            
                            if (bsheet_assign(nresii + CO_i_pos) <= 1) then
                                bsheet_assign(nresii + CO_i_pos) = 2
                            else if (bsheet_assign(nresii + CO_i_pos) == 3) then
                                bsheet_assign(nresii + CO_i_pos) = 5
                            endif

                            if (bsheet_assign(nresii + CO_i_pos + 1) <= 1) then
                                bsheet_assign(nresii + CO_i_pos + 1) = 2
                            else if (bsheet_assign(nresii + CO_i_pos + 1) == 3) then
                                bsheet_assign(nresii + CO_i_pos + 1) = 5
                            end if

                            if (bsheet_assign(nresii + CO_i_pos + 2) <= 1) then
                                bsheet_assign(nresii + CO_i_pos + 2) = 2
                            else if (bsheet_assign(nresii + CO_i_pos + 2) == 3) then
                                bsheet_assign(nresii + CO_i_pos + 2) = 5
                            end if

                        end if
                    end if     
                end do                
                
                
                if ((a>0 .and. b>=0) .or. (a>=0 .and. b>0)) then        ! the two chains are parallel
                    do x=1, nresi                                       ! two chains are equal length, otherwise, need two judgement
                        if ((bsheet_assign(nresii + x) == 3) .or. (bsheet_assign(nresii + x) == 5)) para = para + 1
                    enddo
                    
                    do x=1, nresj
                        if ((bsheet_assign(nresjj + x) == 3) .or. (bsheet_assign(nresjj + x) == 5)) para = para + 1
                    enddo
                else if (( a<0 .and. b<=0) .or. (a<=0 .and. b<0)) then
                    do x=1, nresj
                        if ((bsheet_assign(nresjj + x) == 2) .or. (bsheet_assign(nresjj + x) == 5)) anti = anti + 1
                    enddo
                    
                    do x=1, nresi
                        if ((bsheet_assign(nresii + x) == 2) .or. (bsheet_assign(nresii + x) == 5)) anti = anti + 1
                    enddo                    
                else
                    do x=1, nresi
                        if ((bsheet_assign(nresii + x) == 2)) then
                            anti = anti + 1
                        elseif ((bsheet_assign(nresii + x) == 3)) then
                            para = para + 1
                        endif
                    enddo                    
                    
                    do x=1, nresj
                        if ((bsheet_assign(nresjj + x) == 2)) then
                            anti = anti + 1
                        elseif ((bsheet_assign(nresjj + x) == 3)) then
                            para = para + 1
                        endif
                    enddo
                    
                    !write(*,*) "error in determining direction in pair:", i, j, a, b
                    !write(*,*) "this can happen when a peptide has a turn"
                    !write(*,*) nh, co, nh2, co2
                    undetermined = undetermined + 1     
                endif 

                if (i == j .and. a <= 0) then
                    do m = 1, size(resNHi)
                        NH_i_pos = resNHi(m)
                        CO_j_pos = resCOj(m)
                        diff = abs(NH_i_pos - CO_j_pos)
                        found = .false.
                        CO_i_pos = -1
                        do k = 1, size(resNHi)
                            if (resNHi(k) == CO_j_pos) then
                                found = .true.
                                CO_i_pos = resCOj(k)
                                exit
                            endif
                        enddo
                        
                        if (found .and. CO_i_pos == NH_i_pos) then
                            if (diff >= 3 .and. diff <= 5) then
                                if ((NH_i_pos-2 >= 1) .and. (CO_j_pos+2 <= nresi)) then
                                    nh_basg = bsheet_assign(nresii + NH_i_pos-2)
                                    co_basg = bsheet_assign(nresii + CO_j_pos+2)
                                    if (((nh_basg == 2) .or. (nh_basg == 5)) .and. ((co_basg == 2) .or. (co_basg == 5))) then
                                        turn = turn + 1
                                    endif
                                endif
                            endif
                            
                        else if (.not. found) then 
                            if (diff >= 3 .and. diff <= 5) then
                                if ((NH_i_pos-2 >= 1) .and. (CO_j_pos+2 <= nresi)) then
                                    nh_basg = bsheet_assign(nresii + NH_i_pos-2)
                                    co_basg = bsheet_assign(nresii + CO_j_pos+2)
                                    if (((nh_basg == 2) .or. (nh_basg == 5)) .and. ((co_basg == 2) .or. (co_basg == 5))) then
                                        turn = turn + 2
                                    endif
                                endif
                            endif 
                        endif
                    enddo
                    
                endif
                
                deallocate(nh_i)
                deallocate(nh_j)
                deallocate(co_i)
                deallocate(co_j)
                
            endif
        enddo
    enddo

    counter = 0
    do k=1, total_residues
        if (bsheet_assign(k) >= 2) counter = counter + 1
    enddo
    content = real(counter)/real(total_residues)
    
    turn = turn / 2
    deallocate(bsheet_assign)
    end subroutine bsheet_content
    
    end module my_functions
    
    program bsheet_content_main
    
    use my_functions
    
    implicit none

    character(len=100)                      :: pep1, pep2
    character(len=100)                      :: energy_file, bptnr_file
    integer                                 :: anneal_steps, sim_steps, start_step, end_step
    integer                                 :: nc1, nc2, nb1, nb2, i, j, k, npts, nframes, n1, n2, nres1, nres2
    integer                                 :: anti, para, undetermined, turn
    integer, allocatable                    :: bptnrs(:,:,:), hbmats(:,:,:), hb(:), anti_ratio(:), para_ratio(:)
    
    integer(8), allocatable                 :: collision_nums(:), contents(:)
    real(8)                                 :: ttotal, content
    real(8), allocatable                    :: t_series(:), t(:), temp(:)
    
    
    anneal_steps = 9  !! annealing steps. change value here
    sim_steps = 300    !! equil steps. change value here
    start_step = 0   !! change value here
    
    end_step = anneal_steps + sim_steps 
    
    call read_input("input.txt", pep1, pep2, nc1, nc2, nb1, nb2)
    open(10,file='bsheet_content.txt', status='replace', action='write')
    close(10)
    nres1 = len_trim(pep1)
    nres2 = len_trim(pep2)
    
    
    !!!!!!!!!!!!!!!!! get start time !!!!!!!!!!!!!!!!!!!!!!!!!!
    ttotal = 0   ! the cumulative time after reading each energy file
    if (start_step /= 0 .and. start_step /= 1) then
        do i = 0, start_step - 1
            write(energy_file, '(A,I4.4,A)') 'results/run', i, '.energy'
            call read_energy(energy_file, ttotal, i, anneal_steps, t_series, npts)

            if (npts > 0) then
                ttotal =  t_series(npts)
            end if

            deallocate(t_series)
        end do
    end if
    
    k = 0
    do j = start_step, end_step
        write(bptnr_file, '(A,I4.4,A)') 'results/run', j, '.bptnr'
        write(energy_file, '(A,I4.4,A)') 'results/run', j, '.energy'
        
        call read_energy(energy_file, ttotal, j, anneal_steps, t_series, npts)
        
        if (.not. allocated(t)) then
            allocate(t(size(t_series)))
            t = t_series
        else ! combine t_series (time in one simulaiton) to the whole time period
            n1 = size(t)
            n2 = size(t_series)
    
            allocate(temp(n1 + n2))
            if (n1 > 0) temp(1:n1) = t(1:n1)
            if (n2 > 0) temp(n1+1:n1+n2) = t_series(1:n2)
    
            deallocate(t)
            
            t = temp
            deallocate(temp)
        endif
        
        
        
        call read_bptnr_multi(bptnr_file, nc1, nc2, nb1, nb2, nframes, collision_nums, bptnrs, hbmats)
    

        open(10,file='bsheet_content.txt', access="append")
        do i = 1, nframes
            
            call bsheet_content(bptnrs(:,:,i), hbmats(:,:,i), nc1, nc2, nb1, nb2, nres1, nres2, anti, para, turn, undetermined, content)
            write(10, "(F10.5, 3I10, F15.5)") t_series(i), anti, para, turn, content

        enddo
        close(10)
        
        deallocate(t_series)
        deallocate(collision_nums)
        deallocate(bptnrs)
        deallocate(hbmats)
    enddo
    
    end program bsheet_content_main

