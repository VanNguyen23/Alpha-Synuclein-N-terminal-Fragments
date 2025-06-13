program vmddssp
	implicit none
	integer i, chains, res, totalres,dnwti,j, count, pos(5), countel, counter,k, hairpins, beta 
	integer, allocatable :: resid(:),numchain(:),hairpinchains(:)
	character, allocatable :: dssp(:)
	character colorcode*1, empty*8
	
	open(unit=1,file='DSSP-P3Next-V40A-Run7',status='old')
	do i = 1,9
		read(1,*)
	enddo
	chains = 24
	res = 31
	totalres = chains*res
	allocate(resid(totalres),dssp(totalres),numchain(totalres),hairpinchains(chains))
	do i = 1,totalres
		read(1,*) resid(i),colorcode,empty,dnwti,dssp(i)
		
		numchain(i) = int((i-1)/res)+1
		!print*, i, numchain(i)
	enddo
	print*, ''
	hairpins = 0
	do i = 1, chains
		count = 0
		pos = 0
		countel = 0
		counter = 0
		do j = (i-1)*res+6, i*res-6
			!print*, i, j
			if ((dssp(j).eq.'T').and.(dssp(j).eq.dssp(j-1))) then
				count = count+1
				pos(count) = j
			endif
		enddo
		!print*, i, count
		if ((count.ge.1).and.(count.le.5)) then
			
			do k = 1, 5
				if ((numchain(pos(1)-1-k).eq.i).and.(dssp(pos(1)-1-k).eq.'E')) then
					countel = countel+1
				endif
				if ((numchain(pos(count)+k).eq.i).and.(dssp(pos(count)+k).eq.'E')) then
					counter = counter+1
				endif
			enddo
			!print*, i, countel, counter, pos(count)
			if ((countel.ge.4).and.(counter.ge.4)) then
				hairpins = hairpins+1
				hairpinchains(hairpins) = i
			endif
		endif
	enddo
	beta = 0
	do i = 1, chains
		count = 0
		do j = (i-1)*res+1, i*res
			if (dssp(j).eq.'E') then
				count = count+1
			endif
		enddo
		if (count.gt.3) then
			beta = beta+1
		endif
	enddo
	print*, 'Number of beta-hairpins in the beta-sheet: ',hairpins
	print*, 'Peptide Chains that formed beta-hairpins: ' 
	print*, pack(hairpinchains,hairpinchains/=0)
	print*, 'Number of beta-strains in the beta-sheet: ', beta
	close(1)
end program