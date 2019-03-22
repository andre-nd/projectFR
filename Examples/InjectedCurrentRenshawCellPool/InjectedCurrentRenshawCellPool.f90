program InjectedCurrentRenshawCellPool
    use MotorUnitPoolClass
    use NeuralTractClass
    use InterneuronPoolClass
    use ConfigurationClass
    use ogpf 
    use randomSeedInitialize
    use SynapsesFactoryModule
    use SynapticNoiseClass
    use AfferentPoolClass
    implicit none 
    !integer, parameter :: wp = kind(1.0d0)
    type(Configuration) :: conf
    real(wp), parameter :: PI = 4 * atan(1.0_wp)    
    real(wp) :: dt
    real(wp) :: tf
    integer :: timeLength
    integer :: i, j
    real(wp), dimension(:), allocatable :: t, MNv_mV
    real(wp) :: tic, toc
    type(gpf) :: gp
    real(wp) :: FR
    integer :: GammaOrder 
    character(len = 80) :: pool, group
    character(len = 80) :: filename = 'confInjectedCurrentRenshawCellPool.rmto'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools    
    type(AfferentPool), dimension(:), allocatable:: afferentPools    

    call init_random_seed()

    conf = Configuration(filename)
    allocate(neuralTractPools(0))
    allocate(motorUnitPools(0))
    pool = 'RC'
    group = 'ext'
    allocate(interneuronPools(1))
    interneuronPools(1) = InterneuronPool(conf, pool, group)
    allocate(afferentPools(0))
    synapticNoisePools = synapseFactory(conf, neuralTractPools, &
                motorUnitPools, interneuronPools, afferentPools)
    
    conf = Configuration(filename)
    
    tf = conf%simDuration_ms
    dt = conf%timeStep_ms
    timeLength = int(tf/dt)
    
    allocate(t(timeLength))
    allocate(MNv_mV(timeLength))

    
      
    t = [(dt*(i-1), i=1, timeLength)]
    
    do j = 1, size(synapticNoisePools)
        call synapticNoisePools(j)%reset()
    end do
    
    do j = 1, size(interneuronPools)
        call interneuronPools(j)%reset()
    end do

    call cpu_time(tic)
    do i = 1, size(t)        
        do j = 1, size(synapticNoisePools)
            call synapticNoisePools(j)%atualizePool(t(i))
        end do
        interneuronPools(1)%iInjected(:) = 5.0
        do j = 1, size(interneuronPools)
            call interneuronPools(j)%atualizeInterneuronPool(t(i))
        end do
    end do
    call cpu_time(toc)

    print '(F15.6, A)', toc - tic, ' seconds'

    
    
    
    
    

    
    call interneuronPools(1)%listSpikes()
    

    call gp%title('Interneuron spike instants at the terminal')
    call gp%xlabel('t (s))')
    call gp%ylabel('interneuron index')
    call gp%plot(interneuronPools(1)%poolSomaSpikes(:,1), &
    interneuronPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

    call synapticNoisePools(1)%listSpikes()
    call gp%title('synaptic noise spike instants at the terminal')
    call gp%xlabel('t (s))')
    call gp%ylabel('noise index')
    call gp%plot(synapticNoisePools(1)%poolTerminalSpikes(:,1), &
    synapticNoisePools(1)%poolTerminalSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')
    
end program InjectedCurrentRenshawCellPool