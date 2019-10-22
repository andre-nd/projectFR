program ForceVariability
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
    include 'mkl.fi'
    !integer, parameter :: wp = kind(1.0d0)
    type(Configuration) :: conf
    real(wp), parameter :: pi = 4 * atan(1.0_wp)    
    real(wp) :: dt
    real(wp) :: tf
    integer :: timeLength
    integer :: i, j
    real(wp), dimension(:), allocatable :: t
    real(wp) :: tic, toc
    type(gpf) :: gp
    real(wp), dimension(9) :: FR
    integer, dimension(9) :: GammaOrder 
    real(wp), dimension(5) :: condVelS, condVelFR, condVelFF1, condVelFF2
    character(len = 80) :: pool, muscle, paramTag
    character(len = 80) :: filename = 'confForceVariability.rmto'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools   
    type(AfferentPool), dimension(:), allocatable:: afferentPools    
    type(jointAnkleForceTask) :: ankle
    real(wp) :: angle
    real(wp) , dimension(200):: Amp, phase
    real(wp), dimension(:,:,:), allocatable :: forceSOL, forceMG, forceLG, torque
    integer :: m, l, k
    character(len = 80) :: value1, value2, fileNumber, velNumber
    
    call init_random_seed()

    conf = Configuration(filename)

    tf = conf%simDuration_ms
    dt = conf%timeStep_ms
    timeLength = nint(tf/dt)

    allocate(t(timeLength))
    
    t = [(dt*(i-1), i=1, timeLength)]

    FR = [95.0, 110.0, 125.0, 140.0, &
          155.0, 170.0, 185.0, 200.0, 250.0]

    GammaOrder = [7, 5, 4,4,4,3,2,2,1]

    condVelS = [((88-22)*(i-1)/4+22, i=1, 5)]
    condVelFR = [((102-25.5)*(i-1)/4+25.5, i=1, 5)]
    condVelFF1 = [((104-26)*(i-1)/4+26, i=1, 5)]
    condVelFF2 = [((106-26.5)*(i-1)/4+26.5, i=1, 5)]

    
    allocate(forceSOL(timeLength, size(FR), size(condVelS)))
    allocate(forceMG(timeLength, size(FR), size(condVelS)))
    allocate(forceLG(timeLength, size(FR), size(condVelS)))
    allocate(torque(timeLength, size(FR), size(condVelS)))

    allocate(afferentPools(0))
    allocate(neuralTractPools(1))
    allocate(motorUnitPools(3))
    allocate(interneuronPools(0))
    
    


    do m = 1, size(condVelS)
        paramTag = 'axonDelayCondVel:MG-S' 
        write(value1, '(F15.6)')condVelS(m)
        write(value2, '(F15.6)')condVelFR(m)
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'axonDelayCondVel:SOL-S'

        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'axonDelayCondVel:LG-S'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'axonDelayCondVel:MG-FR'
        
        write(value1, '(F15.6)')condVelFR(m)
        write(value2, '(F15.6)')condVelFF1(m)
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'axonDelayCondVel:SOL-FR'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'axonDelayCondVel:LG-FR'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'axonDelayCondVel:MG-FF'
        
        write(value1, '(F15.6)')condVelFF1(m)
        write(value2, '(F15.6)')condVelFF2(m)
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'axonDelayCondVel:SOL-FF'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'axonDelayCondVel:LG-FF'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        
        
        pool = 'CMExt'
        neuralTractPools(1) = NeuralTract(conf, pool)
        pool = 'SOL'
        motorUnitPools(1) = MotorUnitPool(conf, pool)
        pool = 'MG'
        motorUnitPools(2) = MotorUnitPool(conf, pool)    
        pool = 'LG'
        motorUnitPools(3) = MotorUnitPool(conf, pool)
        ankle = jointAnkleForceTask(conf, motorUnitPools)
        synapticNoisePools = synapseFactory(conf, neuralTractPools, &
                                            motorUnitPools, &
                                            interneuronPools, &
                                            afferentPools)
        
        do l = 1,  size(FR)
            tic =  dsecnd()
            do j = 1, size(neuralTractPools)
                call neuralTractPools(j)%reset()
            end do            
            do j = 1, 3
                call motorUnitPools(j)%reset()
            end do            
            call ankle%reset()
            do i = 1, size(t)
                angle = 0.0
                call ankle%atualizeAnkle(t(i), angle)  
                if (mod(t(i),1000.0_wp)==0) then
                    print *, t(i)
                end if
                do j = 1, size(neuralTractPools)
                    call neuralTractPools(j)%atualizePool(t(i), FR(l), GammaOrder(l))
                end do
                do j = 1, 3
                    call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 0.0_wp, 0.0_wp)
                end do
                call ankle%computeTorque(t(i))                
            end do    
            toc =  dsecnd()

            print '(F15.6, A)', toc - tic, ' seconds'
            
            do j = 1, 3
                call motorUnitPools(j)%listSpikes()
            end do

            call neuralTractPools(1)%listSpikes()
            
            
            
            
            do i = 1, timeLength
                forceSOL(i,l, m) = motorUnitPools(1)%NoHillMuscle%force(i)
                forceMG(i,l, m) = motorUnitPools(2)%NoHillMuscle%force(i)
                forceLG(i,l, m) = motorUnitPools(3)%NoHillMuscle%force(i)
                torque(i,l, m) = ankle%ankleTorque_Nm(i)
            end do

            write(fileNumber, '(F15.2)')FR(l)
            call gp%title('Ankle torque' // fileNumber)
            call gp%xlabel('t (ms))')
            call gp%ylabel('torque (N.m)')
            call gp%plot(t, ankle%ankleTorque_Nm, 'with line lw 2 lc rgb "#0008B0"')            

            call gp%title('force SOLeus' // fileNumber)
            call gp%xlabel('t (ms))')
            call gp%ylabel('force (N)')
            call gp%plot(t, motorUnitPools(1)%NoHillMuscle%force, 'with line lw 2 lc rgb "#0008B0"')

            call gp%title('force MG' // fileNumber)
            call gp%xlabel('t (ms))')
            call gp%ylabel('force (N)')
            call gp%plot(t, motorUnitPools(2)%NoHillMuscle%force, 'with line lw 2 lc rgb "#0008B0"')

            call gp%title('force LG' // fileNumber)
            call gp%xlabel('t (ms))')
            call gp%ylabel('force (N)')
            call gp%plot(t, motorUnitPools(3)%NoHillMuscle%force, 'with line lw 2 lc rgb "#0008B0"')
            
            call gp%title('MN spike instants at the soma')
            call gp%xlabel('t (s))')
            call gp%ylabel('Motoneuron index')
            call gp%plot(motorUnitPools(1)%poolSomaSpikes(:,1), &
            motorUnitPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

            call gp%title('MG MN spike instants at the soma')
            call gp%xlabel('t (s))')
            call gp%ylabel('Motoneuron index')
            call gp%plot(motorUnitPools(2)%poolSomaSpikes(:,1), &
            motorUnitPools(2)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

            call gp%title('LG MN spike instants at the soma')
            call gp%xlabel('t (s))')
            call gp%ylabel('Motoneuron index')
            call gp%plot(motorUnitPools(3)%poolSomaSpikes(:,1), &
            motorUnitPools(3)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

            call gp%title('NT spike instants at the soma')
            call gp%xlabel('t (s))')
            call gp%ylabel('NT index')
            call gp%plot(neuralTractPools(1)%poolTerminalSpikes(:,1), &
            neuralTractPools(1)%poolTerminalSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')
            
            ! write(velNumber, '(I2)')
            ! filename = 'spikesSOL_FR'// fileNumber //'vel' // velNumber // '.txt'
            ! open(1, file=filename, status = 'replace')
            ! do i = 1, size(t)
            !     write(1, '(F15.6, 1X)') (motorUnitPools(1)%poolTerminalSpikes, j=1,size(condVelS))
            ! end do
            ! close(1)

            ! filename = 'spikesMG_FR'// fileNumber //'vel' // velNumber // '.txt'
            ! open(1, file=filename, status = 'replace')
            ! do i = 1, size(t)        
            !     write(1, '(F15.6, 1X)') (motorUnitPools(2)%poolTerminalSpikes, j=1,size(condVelS))
            ! end do
            ! close(1)

            ! filename = 'spikesLG_FR'// fileNumber //'vel' // velNumber // '.txt'
            ! open(1, file=filename, status = 'replace')
            ! do i = 1, size(t)
            !     write(1, '(F15.6, 1X)') (motorUnitPools(3)%poolTerminalSpikes, j=1,size(condVelS))
            ! end do
            ! close(1)
            
            
            
        end do
        deallocate(synapticNoisePools)
        if (allocated(motorUnitPools)) deallocate(motorUnitPools)
        allocate(motorUnitPools(3))
        if (allocated(neuralTractPools)) deallocate(neuralTractPools)
        allocate(neuralTractPools(1))
        
    end do
    do l = 1, size(FR)
        write(fileNumber, '(F4.0)')FR(l)
        filename = 'forceSOL_FR'// trim(fileNumber) //'linear.txt'
        open(1, file=filename, status = 'replace') 
        do i = 1, size(t)        
            write(1, '(F15.6, 1X,F15.6,1X,F15.6,1X,F15.6,1X,F15.6)') (forceSOL(i,l,j), j=1,size(condVelS))
        end do
        close(1)

        filename = 'forceMG_FR'// trim(fileNumber) //'linear.txt'
        open(1, file=filename, status = 'replace') 
        do i = 1, size(t)        
            write(1, '(F15.6, 1X,F15.6,1X,F15.6,1X,F15.6,1X,F15.6)') (forceMG(i,l,j), j=1,size(condVelS))
        end do
        close(1)

        filename = 'forceLG_FR'// trim(fileNumber) //'linear.txt'
        open(1, file=filename, status = 'replace') 
        do i = 1, size(t)        
            write(1, '(F15.6, 1X,F15.6,1X,F15.6,1X,F15.6,1X,F15.6)') (forceLG(i,l,j), j=1,size(condVelS))
        end do
        close(1)

        filename = 'torque_FR'// trim(fileNumber) //'linear.txt'
        open(1, file=filename, status = 'replace') 
        do i = 1, size(t)        
            write(1, '(F15.6, 1X,F15.6,1X,F15.6,1X,F15.6,1X,F15.6)') (torque(i,l,j), j=1,size(condVelS))
        end do
        close(1)
    end do
    
    
    
    
end program ForceVariability