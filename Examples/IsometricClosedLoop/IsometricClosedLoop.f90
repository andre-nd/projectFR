program IsometricClosedLoop    
    use MotorUnitPoolClass
    use NeuralTractClass
    use InterneuronPoolClass
    use SynapticNoiseClass
    use ConfigurationClass
    use ogpf 
    use randomSeedInitialize
    use SynapsesFactoryModule
    use jointAnkleForceTaskClass
    use AfferentPoolClass
    implicit none 
    !integer, parameter :: wp = kind(1.0d0)
    type(Configuration) :: conf
    real(wp), parameter :: pi = 4 * atan(1.0_wp)    
    real(wp) :: dt
    real(wp) :: tf
    integer :: timeLength
    integer :: i, j
    real(wp), dimension(:), allocatable :: t, IaFR
    real(wp) :: tic, toc
    type(gpf) :: gp
    real(wp) :: FR
    integer :: GammaOrder 
    character(len = 80) :: pool, muscle
    character(len = 80) :: filename = 'confIsometricClosedLoop.rmto'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools   
    type(AfferentPool), dimension(:), allocatable:: afferentPools    
    type(jointAnkleForceTask) :: ankle
    real(wp) :: angle
    
    
    call init_random_seed()

    conf = Configuration(filename)
    allocate(afferentPools(1))
    pool = 'Ia'
    muscle = 'SOL'
    afferentPools(1) = AfferentPool(conf, pool, muscle)
    allocate(neuralTractPools(1))
    pool = 'CMExt'
    neuralTractPools(1) = NeuralTract(conf, pool)
    allocate(motorUnitPools(1))
    pool = 'SOL'
    motorUnitPools(1) = MotorUnitPool(conf, pool)    
    ankle = jointAnkleForceTask(conf, motorUnitPools)
    allocate(interneuronPools(0))

    synapticNoisePools = synapseFactory(conf, neuralTractPools, &
                                        motorUnitPools, &
                                        interneuronPools, &
                                        afferentPools)
    
    tf = conf%simDuration_ms
    dt = conf%timeStep_ms
    timeLength = int(tf/dt)

    allocate(t(timeLength))
    allocate(IaFR(timeLength))
    
    t = [(dt*(i-1), i=1, timeLength)]
    
    FR = 1000.0/12.0
    GammaOrder = 10
    
    do j = 1, size(neuralTractPools)
        call neuralTractPools(j)%reset()
    end do
    do j = 1, size(motorUnitPools)
        call motorUnitPools(j)%reset()
        call afferentPools(j)%reset()
    end do

    call cpu_time(tic)
    do i = 2, size(t)
        angle = 0.1*sin(2*pi*t(i)/1000.0)    
        call ankle%atualizeAnkle(t(i), angle)      
        do j = 1, size(neuralTractPools)
            call neuralTractPools(j)%atualizePool(t(i), FR, GammaOrder)
        end do
        do j = 1, size(motorUnitPools)
            call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 32.0_wp, 32.0_wp)
            call afferentPools(j)%atualizeAfferentPool(t(i), motorUnitPools(j)%spindle%IaFR_Hz)            
        end do
        call ankle%computeTorque(t(i))
        IaFR(i) = motorUnitPools(1)%spindle%IaFR_Hz
    end do    
    call cpu_time(toc)

    print '(F15.6, A)', toc - tic, ' seconds'

    
    
    
    call neuralTractPools(1)%listSpikes()
    call motorUnitPools(1)%listSpikes()
    call afferentPools(1)%listSpikes()
    call motorUnitPools(1)%getMotorUnitPoolEMG()
    
    call gp%title('MN spike instants at the soma')
    call gp%xlabel('t (s))')
    call gp%ylabel('Motoneuron index')
    call gp%plot(motorUnitPools(1)%poolSomaSpikes(:,1), motorUnitPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

    call gp%title('MN spike instants at the terminal')
    call gp%xlabel('t (s))')
    call gp%ylabel('Motoneuron index')
    call gp%plot(motorUnitPools(1)%poolTerminalSpikes(:,1), &
    motorUnitPools(1)%poolTerminalSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

    call gp%title('Descending command spike instants at the terminal')
    call gp%xlabel('t (s))')
    call gp%ylabel('Descending command index')
    call gp%plot(neuralTractPools(1)%poolTerminalSpikes(:,1), &
    neuralTractPools(1)%poolTerminalSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

    call gp%title('afferents spike instants at the terminal')
    call gp%xlabel('t (s))')
    call gp%ylabel('Ia index')
    call gp%plot(afferentPools(1)%poolTerminalSpikes(:,1), &
    afferentPools(1)%poolTerminalSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')


    call gp%title('Muscle force')
    call gp%xlabel('t (ms))')
    call gp%ylabel('Force (N)')
    call gp%plot(t, motorUnitPools(1)%HillMuscle%force, 'with line lw 2 lc rgb "#0008B0"')

    call gp%title('Muscle length')
    call gp%xlabel('t (ms))')
    call gp%ylabel('length (m)')
    call gp%plot(t, motorUnitPools(1)%HillMuscle%length_m, 'with line lw 2 lc rgb "#0008B0"')

    call gp%title('Muscle velocity')
    call gp%xlabel('t (ms))')
    call gp%ylabel('velocity (m/ms)')
    call gp%plot(t, motorUnitPools(1)%HillMuscle%velocity_m_ms, 'with line lw 2 lc rgb "#0008B0"')

    call gp%title('Ankle angle')
    call gp%xlabel('t (ms))')
    call gp%ylabel('angle (degree)')
    call gp%plot(t, ankle%ankleAngle_rad, 'with line lw 2 lc rgb "#0008B0"')
    
    call gp%title('Ia Firing rate')
    call gp%xlabel('t (ms))')
    call gp%ylabel('Ia FR (Hz)')
    call gp%plot(t, IaFR, 'with line lw 2 lc rgb "#0008B0"')
end program IsometricClosedLoop