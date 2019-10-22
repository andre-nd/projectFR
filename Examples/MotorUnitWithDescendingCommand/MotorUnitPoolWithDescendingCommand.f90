program MotorUnitPoolWithDescendingCommand
    use MotorUnitPoolClass
    use NeuralTractClass
    use InterneuronPoolClass
    use SynapticNoiseClass
    use ConfigurationClass
    use ogpf 
    use randomSeedInitialize
    use MuscleNoHillClass
    use CharacterArrayClass
    use CharacterMatrixClass
    use QueueClass
    use SynapsesFactoryModule
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
    real(wp), dimension(:,:), allocatable :: unitsEMG
    real(wp) :: tic, toc
    type(gpf) :: gp
    real(wp) :: FR
    integer :: GammaOrder 
    character(len = 80) :: pool = 'CMExt'
    character(len = 80) :: filename = 'confMNPoolWithDescendingCommand.rmto'
    type(CharacterArray) :: dynamics  
    Character(len = 80) :: synapseDyn
    type(Queue) :: spikesQueue
    integer :: newItem
    logical, dimension(10) :: logList
    integer, dimension(10) :: logCount
    type(CharacterMatrix) :: Synapses
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools     
    type(AfferentPool), dimension(:), allocatable:: afferentPools     
    character(len = 80) :: muindex
    call init_random_seed()

    conf = Configuration(filename)
    allocate(neuralTractPools(1))
    neuralTractPools(1) = NeuralTract(conf, pool)
    pool = 'SOL'
    allocate(motorUnitPools(1))
    motorUnitPools(1) = MotorUnitPool(conf, pool)    
    allocate(interneuronPools(0))
    allocate(afferentPools(0))
    synapticNoisePools = synapseFactory(conf, neuralTractPools, &
                                        motorUnitPools, interneuronPools, &
                                        afferentPools)
    
    
    
    tf = conf%simDuration_ms
    dt = conf%timeStep_ms
    timeLength = int(tf/dt)
    
    allocate(t(timeLength))
    allocate(MNv_mV(timeLength))

    do j = 1, size(neuralTractPools)
        call neuralTractPools(j)%reset()
    end do

    call motorUnitPools(1)%reset()
    
      
    t = [(dt*(i-1), i=1, timeLength)]
    
    FR = 60
    GammaOrder = 10

    call cpu_time(tic)
    do i = 1, size(t)        
        do j = 1, size(neuralTractPools)
            call neuralTractPools(j)%atualizePool(t(i), FR, GammaOrder)            
        end do
        do j = 1, size(motorUnitPools)
            call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 30.0_wp, 30.0_wp)            
            MNv_mV(i) = motorUnitPools(j)%v_mV(30)        
        end do
    end do
    call cpu_time(toc)

    print '(F15.6, A)', toc - tic, ' seconds'    
    
    call neuralTractPools(1)%listSpikes()
    call motorUnitPools(1)%listSpikes()
    call motorUnitPools(1)%getMotorUnitPoolEMG()
    allocate(unitsEMG(20, size(t)))
    do j = 1, size(t)
        do i = 1, 20
            unitsEMG(i,j) = motorUnitPools(1)%unit(i)%getEMG(t(j))
        end do
    end do
    
    call gp%title('EMG')
    call gp%xlabel('t (ms))')
    call gp%ylabel('EMG (mV)')
    call gp%plot(t, motorUnitPools(1)%emg, 'with line lw 2 lc rgb "#0008B0"')  

    do i = 1, 20
        write(muindex, '(I3)')i
        call gp%title('EMG' // muindex )
        call gp%xlabel('t (ms))')
        call gp%ylabel('EMG (mV)')
        call gp%plot(t, unitsEMG(i,:), 'with line lw 2 lc rgb "#0008B0"')  
    end do
    
    ! call gp%title('MN spike instants at the soma')
    ! call gp%xlabel('t (s))')
    ! call gp%ylabel('Motoneuron index')
    ! call gp%plot(motorUnitPools(1)%poolSomaSpikes(:,1), motorUnitPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

    ! call gp%title('MN spike instants at the terminal')
    ! call gp%xlabel('t (s))')
    ! call gp%ylabel('Motoneuron index')
    ! call gp%plot(motorUnitPools(1)%poolTerminalSpikes(:,1), &
    ! motorUnitPools(1)%poolTerminalSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

    ! call gp%title('MN spike instants at the terminal')
    ! call gp%xlabel('t (s))')
    ! call gp%ylabel('Motoneuron index')
    ! call gp%plot(neuralTractPools(1)%poolTerminalSpikes(:,1), &
    ! neuralTractPools(1)%poolTerminalSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

    ! call gp%title('Muscle force')
    ! call gp%xlabel('t (ms))')
    ! call gp%ylabel('Force (N)')
    ! call gp%plot(t, motorUnitPools(1)%NoHillMuscle%force, 'with line lw 2 lc rgb "#0008B0"')
    
    
end program MotorUnitPoolWithDescendingCommand