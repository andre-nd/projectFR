program HReflexMWave
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
    type(Configuration), target :: conf
    real(wp), parameter :: pi = 4 * atan(1.0_wp)    
    real(wp) :: dt
    real(wp) :: tf
    integer :: timeLength
    integer :: i, j, k
    real(wp), dimension(:), allocatable :: t
    real(wp) :: tic, toc
    type(gpf) :: gp
    character(len = 80) :: pool, muscle, paramTag
    character(len = 80) :: filename = 'confH-ReflexM-WaveCurve.rmto'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools   
    type(AfferentPool), dimension(:), allocatable:: afferentPools    
    type(jointAnkleForceTask) :: ankle
    real(wp) :: angle
    real(wp) , dimension(200):: Amp, phase
    integer :: Nsim
    real(wp) :: FirstStim, LastStim
    real(wp), dimension(:), allocatable :: Mp, Hp, Stim
    real(wp), dimension(:,:), allocatable :: emg
    character(len = 80) :: value1, value2
    real(wp) :: spindleFR
    character(len = 3) :: fileNumber

    Nsim = 35
    FirstStim = 9.0
    LastStim = 25
    
    call init_random_seed()

    conf = Configuration(filename)

    tf = conf%simDuration_ms
    dt = conf%timeStep_ms
    timeLength = int(tf/dt)

    allocate(t(timeLength))
    
    t = [(dt*(i-1), i=1, timeLength)]
    
    allocate(emg(timeLength, Nsim))
    allocate(Mp(Nsim))
    allocate(Hp(Nsim))
    allocate(Stim(Nsim))   
    
    

    allocate(afferentPools(2))
    pool = 'Ia'
    muscle = 'SOL'
    afferentPools(1) = AfferentPool(conf, pool, muscle)

    pool = 'Ia'
    muscle = 'LG'
    afferentPools(2) = AfferentPool(conf, pool, muscle)

    allocate(neuralTractPools(0))
    
    allocate(motorUnitPools(1))
    pool = 'SOL'
    motorUnitPools(1) = MotorUnitPool(conf, pool)    
    

    ankle = jointAnkleForceTask(conf, motorUnitPools)
    allocate(interneuronPools(0))

    synapticNoisePools = synapseFactory(conf, neuralTractPools, &
                                        motorUnitPools, &
                                        interneuronPools, &
                                        afferentPools)
    

    do j = 1, Nsim
        Stim(j) = FirstStim + j * (LastStim - FirstStim) / (Nsim - 1.0)
        paramTag = 'stimIntensity_PTN'
        write(value1, '(F15.6)')Stim(j)
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1,value2)
        
        do i = 1, size(motorUnitPools(1)%unit)
            call motorUnitPools(1)%unit(i)%createStimulus()
        end do
        
        do i = 1, size(afferentPools(1)%unit)
            call afferentPools(1)%unit(i)%createStimulus()
        end do
        do i = 1, size(afferentPools(2)%unit)
            call afferentPools(2)%unit(i)%createStimulus()
        end do

        call motorUnitPools(1)%reset()
        call afferentPools(1)%reset()
        call afferentPools(2)%reset()
        spindleFR = 0.0
        call cpu_time(tic)
        do i = 1, timeLength
            call motorUnitPools(1)%atualizeMotorUnitPool(t(i), 32.0_wp, 32.0_wp)
            call afferentPools(1)%atualizeAfferentPool(t(i),spindleFR)
            call afferentPools(2)%atualizeAfferentPool(t(i),spindleFR)
        end do
        call cpu_time(toc)
        print '(F15.6, A)', toc - tic, ' seconds'
        call motorUnitPools(1)%getMotorUnitPoolEMG()
        Mp(j) = maxval(motorUnitPools(1)%emg(1:500)) - minval(motorUnitPools(1)%emg(1:500))
        Hp(j) = maxval(motorUnitPools(1)%emg(501:)) - minval(motorUnitPools(1)%emg(501:))
        print '(F15.6)', Mp(j)
        print '(F15.6)', Hp(j)
        emg(:,j) = motorUnitPools(1)%emg
        call motorUnitPools(1)%listSpikes()
        call afferentPools(1)%listSpikes()
        call afferentPools(2)%listSpikes()
        write(fileNumber, '(I2)')j
        filename = 'spikes'// fileNumber //'.txt'
        open(1, file=filename, status = 'replace') 
        
        do i = 1, size(motorUnitPools(1)%poolTerminalSpikes, 1)            
            write(1, '(F15.6, 1X, F15.1)') motorUnitPools(1)%poolTerminalSpikes(i,1),&
                        motorUnitPools(1)%poolTerminalSpikes(i,2)
        end do
        close(1)
        
        call motorUnitPools(1)%reset()
        call afferentPools(1)%reset()
        call afferentPools(2)%reset()
    end do
   
    call gp%title('EMG')
    call gp%xlabel('t (ms))')
    call gp%ylabel('emg (mV)')
    call gp%plot(t, emg(:,10),  'with line lw 2 lc rgb "#0008B0"')
        
    call gp%title('Hp')
    call gp%xlabel('Stimulus intensity (mA)')
    call gp%ylabel('EMG peak (mV)')
    call gp%plot(Stim, Hp, 'with line lw 2 lc rgb "#0008B0"')
    
    call gp%title('Mp')
    call gp%xlabel('Stimulus intensity (mA)')
    call gp%ylabel('EMG peak (mV)')
    call gp%plot(Stim, Mp, 'with line lw 2 lc rgb "#0008B0"')

    call gp%title('Mp-Hp')
    call gp%xlabel('Stimulus intensity (mA)')
    call gp%ylabel('EMG peak (mV)')
    call gp%plot(x1=Stim, y1=Mp, ls1='with line lw 2 lc rgb "blue"', &
                 x2 = Stim, y2 = Hp, ls2='with line lw 2 lc rgb "red"')

    print '(F15.6)', maxval(Hp)/maxval(Mp)   

    
end program HReflexMWave