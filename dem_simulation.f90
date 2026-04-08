program dem_simulation
    use omp_lib
    implicit none
    ! =====================
    ! PARAMETERS
    ! =====================
    integer, parameter :: TEST_FREEFALL=1, TEST_CONSTVEL=2, TEST_BOUNCE=3,MULTIPLE=4
    integer :: test_case
    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: dt = 1.0d-4
    real(dp), parameter :: t_end = 5.0d0
    real(dp), parameter :: kn = 1.0d5
    real(dp), parameter :: gamma_n = 10.0d0
    integer, parameter :: N = 200
    real(dp), parameter :: Lx = 1.0d0, Ly = 1.0d0, Lz = 1.0d0
    real(8) :: t_start, t_fin
    ! =====================
    ! VARIABLES
    ! =====================
    
    real(dp) :: x(N,3), v(N,3), f(N,3)
    real(dp) :: m(N), R(N)
    real(dp) :: g(3)
    real(dp) :: t
    real(dp) :: z0, z_analytical, error
    integer :: i
    integer :: step, nsteps
    test_case = MULTIPLE  ! change to 2 or 3
    nsteps = int(t_end / dt)
    ! =====================
    ! INITIALIZATION
    ! =====================
    select case(test_case)

    case(TEST_FREEFALL)
        print *, "Running Free Fall Test"
        x(1,:) = (/0.5_dp, 0.5_dp, 0.8_dp/)
        v(1,:) = 0.0_dp
        m(1) = 1.0_dp
        R(1) = 0.02_dp
        g = (/0.0_dp, 0.0_dp, -9.81_dp/)
        open(20, file="freefall.csv", status="replace")

    case(TEST_CONSTVEL)
        print *, "Running Constant Velocity Test"
        x(1,:) = (/0.5_dp, 0.5_dp, 0.5_dp/)
        v(1,:) = (/1.0_dp, 0.0_dp, 0.0_dp/)
        m(1) = 1.0_dp
        R(1) = 0.02_dp
        g = (/0.0_dp, 0.0_dp, 0.0_dp/)
        open(20, file="constvel.csv", status="replace")

    case(TEST_BOUNCE)
        print *, "Running Bounce Test"
        x(1,:) = (/0.5_dp, 0.5_dp, 0.5_dp/)
        v(1,:) = 0.0_dp
        m(1) = 1.0_dp
        R(1) = 0.02_dp
        g = (/0.0_dp, 0.0_dp, -9.81_dp/)
        open(20, file="bounce.csv", status="replace")

    case(MULTIPLE)
        print *, "N=",N
        m = 1.0_dp
        R = 0.02_dp
        call initialize_particles(x,v,m,R)
        g = (/0.0_dp, 0.0_dp, -9.81_dp/)

    end select
    z0 = x(1,3)
    ! =====================
    ! TIME LOOP
    ! =====================
    call cpu_time(t_start)
    do step = 1, nsteps
        
        call zero_forces(f)

        call add_gravity(f, m, g)

        call particle_contacts(x, v, R, f)

        if (test_case == TEST_BOUNCE .OR. test_case==MULTIPLE) then
            call wall_contacts(x, v, R, f)
        end if

        call integrate(x, v, f, m)
        t = step * dt
        z_analytical = z0 - 0.5d0 * 9.81d0 * t**2
        error = abs(x(1,3) - z_analytical)
        select case(test_case)

        case(TEST_FREEFALL)
            write(20,*) t, x(1,3)

        case(TEST_CONSTVEL)
            write(20,*) t, x(1,1), v(1,1)

        case(TEST_BOUNCE)
            write(20,*) t, x(1,3), v(1,3)
        case(MULTIPLE)
            if (mod(step, 100) == 0) then
            call write_vtk(x, step)
            end if
        end select

    end do

    print *, "Simulation complete."
    close(20)
    call cpu_time(t_fin)
    print *, "Total runtime (s):", t_fin - t_start
    print *, "Threads:", omp_get_max_threads()
contains

subroutine initialize_particles(x, v, m, R)
    implicit none
    real(dp), intent(in) :: m(:),R(:)
    real(8), intent(out) :: x(:,:), v(:,:)

    integer :: i, j, max_attempts, attempt
    real(8) :: rand(3), dist, rij(3)
    logical :: overlap

    call random_seed()

    max_attempts = 10000

    do i = 1, size(m)

        attempt = 0

        do
            attempt = attempt + 1
            if (attempt > max_attempts) then
                print *, "ERROR: Could not place particle ", i
                stop
            end if

            ! Generate random position inside box (with margin)
            call random_number(rand)
            x(i,1) = R(i) + (Lx - 2.0d0*R(i)) * rand(1)
            x(i,2) = R(i) + (Ly - 2.0d0*R(i)) * rand(2)
            x(i,3) = R(i) + (Lz - 2.0d0*R(i)) * rand(3)

            overlap = .false.

            ! Check overlap with previous particles
            do j = 1, i-1
                rij = x(i,:) - x(j,:)
                dist = sqrt(sum(rij**2))

                if (dist < (R(i) + R(j))) then
                    overlap = .true.
                    exit
                end if
            end do

            if (.not. overlap) exit

        end do

        ! Initialize velocity and properties
        v(i,:) = 0.0d0  

    end do

end subroutine

subroutine zero_forces(f)
    implicit none
    real(kind=8), intent(out) :: f(:,:)
    f = 0.0d0
end subroutine

subroutine add_gravity(f, m, g)
    implicit none
    real(kind=8), intent(inout) :: f(:,:)
    real(kind=8), intent(in) :: m(:), g(3)
    integer :: i

    do i = 1, size(m)
        f(i,:) = f(i,:) + m(i) * g
    end do
end subroutine

subroutine particle_contacts(x, v, R, f)
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp), intent(in) :: x(:,:), v(:,:), R(:)
    real(dp), intent(inout) :: f(:,:)

    integer :: i, j
    real(dp) :: rij(3), nij(3), vij(3)
    real(dp) :: dij, delta, vn, Fn

    !$omp parallel do private(i,j,rij,nij,vij,dij,delta,vn,Fn) shared(x,v,R,f)
    do i = 1, size(R)
        do j = i+1, size(R)

            rij = x(j,:) - x(i,:)
            dij = sqrt(sum(rij**2))

            if (dij > 0.0d0) then
                nij = rij / dij
            else
                cycle
            end if

            delta = R(i) + R(j) - dij

            if (delta > 0.0d0) then
                vij = v(j,:) - v(i,:)
                vn = dot_product(vij, nij)

                Fn = max(0.0d0, kn*delta - gamma_n*vn)

                ! --- ATOMIC UPDATES ---
                !$omp atomic
                f(i,1) = f(i,1) - Fn * nij(1)
                !$omp atomic
                f(i,2) = f(i,2) - Fn * nij(2)
                !$omp atomic
                f(i,3) = f(i,3) - Fn * nij(3)

                !$omp atomic
                f(j,1) = f(j,1) + Fn * nij(1)
                !$omp atomic
                f(j,2) = f(j,2) + Fn * nij(2)
                !$omp atomic
                f(j,3) = f(j,3) + Fn * nij(3)

            end if

        end do
    end do
    !$omp end parallel do
end subroutine

subroutine wall_contacts(x, v, R, f)
    implicit none
    real(kind=8), intent(in) :: x(:,:), v(:,:), R(:)
    real(kind=8), intent(inout) :: f(:,:)

    real(kind=8), parameter :: kn = 1.0d5, gamma_n = 10.0d0
    real(kind=8), parameter :: Lx=1.0d0, Ly=1.0d0, Lz=1.0d0

    integer :: i
    real(kind=8) :: delta, vn

    do i = 1, size(R)

        ! =====================
        ! Z walls
        ! =====================

        ! Bottom wall (z = 0)
        delta = R(i) - x(i,3)
        if (delta > 0.0d0) then
            vn = v(i,3)
            f(i,3) = f(i,3) + max(0.0d0, kn*delta - gamma_n*vn)
        end if

        ! Top wall (z = Lz)
        delta = R(i) - (Lz - x(i,3))
        if (delta > 0.0d0) then
            vn = -v(i,3)
            f(i,3) = f(i,3) - max(0.0d0, kn*delta - gamma_n*vn)
        end if


        ! =====================
        ! X walls
        ! =====================

        ! Left wall (x = 0)
        delta = R(i) - x(i,1)
        if (delta > 0.0d0) then
            vn = v(i,1)
            f(i,1) = f(i,1) + max(0.0d0, kn*delta - gamma_n*vn)
        end if

        ! Right wall (x = Lx)
        delta = R(i) - (Lx - x(i,1))
        if (delta > 0.0d0) then
            vn = -v(i,1)
            f(i,1) = f(i,1) - max(0.0d0, kn*delta - gamma_n*vn)
        end if


        ! =====================
        ! Y walls
        ! =====================

        ! Front wall (y = 0)
        delta = R(i) - x(i,2)
        if (delta > 0.0d0) then
            vn = v(i,2)
            f(i,2) = f(i,2) + max(0.0d0, kn*delta - gamma_n*vn)
        end if

        ! Back wall (y = Ly)
        delta = R(i) - (Ly - x(i,2))
        if (delta > 0.0d0) then
            vn = -v(i,2)
            f(i,2) = f(i,2) - max(0.0d0, kn*delta - gamma_n*vn)
        end if
    end do

end subroutine

subroutine integrate(x, v, f, m)
    implicit none
    real(kind=8), intent(inout) :: x(:,:), v(:,:)
    real(kind=8), intent(in) :: f(:,:), m(:)

    real(kind=8), parameter :: dt = 1.0d-4
    integer :: i

    do i = 1, size(m)
        v(i,:) = v(i,:) + (f(i,:) / m(i)) * dt
        x(i,:) = x(i,:) + v(i,:) * dt
    end do

end subroutine

subroutine write_vtk(x, step)
    implicit none
    integer, intent(in) :: step
    real(8), intent(in) :: x(:,:)

    integer :: i, N
    character(len=100) :: filename

    N = size(x,1)

    write(filename, '(A,I5.5,A)') "particles_", step, ".vtk"

    open(unit=10, file=filename, status='replace')

    write(10,'(A)') '# vtk DataFile Version 3.0'
    write(10,'(A)') 'DEM particles'
    write(10,'(A)') 'ASCII'
    write(10,'(A)') 'DATASET POLYDATA'
    write(10,'(A,I0,A)') 'POINTS ', N, ' float'

    do i = 1, N
        write(10,'(3F12.6)') x(i,1), x(i,2), x(i,3)
    end do

    ! write(10,*) "VERTICES", N, 2*N
    ! do i = 1, N
    !     write(10,*) 1, i-1
    ! end do

    close(10)
end subroutine

end program dem_simulation