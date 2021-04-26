        program opt
        implicit real*8(a-h,o-z)
        parameter (ndim=27)
        parameter (natoms=9)
        parameter (ftol=1e-50)
        parameter (nmono=3)
        dimension psips(ndim),psips_opt(ndim),coord(3,natoms),y(ndim),
     1  grad(ndim),coord2(3,natoms),x(natoms,3),x_opt(natoms,3)
        open(unit=8,file='pud.dat',status='old')
        open(unit=9,file='pud-opt-f50-2.dat',status='unknown')
        open(unit=10,file='pud-opt2.xyz',status='unknown')
C		Given a geometry, finds the minimum on a potential energy
C		surface using the conjugate gradient method. 
C		Does not have to be the global minimum, can be local
C		Inputs:
C		Geometry of a molecule in units of angstrom
C		Outputs: 
C		Coordinates of the minimum energy structure.
        read(8,*) ((coord(j,i),j=1,3),i=1,natoms) 
        ip = 0
        do i =1,natoms
            do j = 1,3
                ip = ip + 1
                y(ip) = coord(j,i)/0.52917721067
            enddo
        enddo
        do j = 1,ndim
            psips(j) = y(j)
        enddo
        do j = 1,ndim
            psips(j) = psips(j)*0.52917721067
        enddo
        call calcpot(nmono,v,psips(:))
        v =v/627.509474
        do j = 1,ndim
            psips(j) = psips(j)/0.52917721067
        enddo
        do j = 1,ndim
            psips_opt(j) = psips(j)
        enddo
        call cpu_time(start)
        avg_iter = 0
C		Details of functions here can be found in 
C		Numerical Recipes in Fortran 77 Section 10.6
        call FRPRMN(psips_opt(:),ndim,ftol,iter,fret)
        print *, iter,fret*219474.6,v*219474.6
        avg_iter = avg_iter + dfloat(iter)
        print *, avg_iter,'average amount of iterations'
        call cpu_time(finish)
        total = finish-start
        print *, total,'seconds'
        ip = 0
        do i =1 ,natoms
            do j = 1,3
                ip = ip + 1
                coord2(j,i) = psips_opt(ip)
            enddo
        enddo
        write(9,*) ((coord2(j,i)*0.52917721067,j=1,3),i=1,natoms)
        ip = 0
        do j = 1,natoms
            do k = 1,3
                ip = ip + 1
                x(j,k) = psips(ip)
                x_opt(j,k) = psips_opt(ip)
            enddo
        enddo
        write (10,*) natoms
        write (10,*) 'before optimzation',v*219474.6
        do j =1,natoms
            if ((j.eq.1).or.(j.eq.4).or.(j.eq.7)) then
                write(10,*) 'O',(x(j,k)*0.529177,k=1,3)
            else
                write(10,*) 'H',(x(j,k)*0.529177,k=1,3)
            endif
        enddo
        write (10,*) ' '
        write (10,*) natoms
        write (10,*) 'after optimzation',fret*219474.6
        do j =1,natoms
            if ((j.eq.1).or.(j.eq.4).or.(j.eq.7)) then
                write(10,*) 'O',(x_opt(j,k)*0.529177,k=1,3)
            else
                write(10,*) 'H',(x_opt(j,k)*0.529177,k=1,3)
            endif
        enddo
        write (10,*) ' '
        call calc_grad(psips_opt(:),grad(:))
C		Checks that geometry is a minimum based on gradient is 0.
        do i = 1,ndim
            print *, grad(i),i
        enddo
        end program
