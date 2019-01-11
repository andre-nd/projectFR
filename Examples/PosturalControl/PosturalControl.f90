program PosturalControl
    use MotorUnitPoolClass
    use NeuralTractClass
    use InterneuronPoolClass
    use SynapticNoiseClass
    use ConfigurationClass
    use ogpf 
    use randomSeedInitialize
    use SynapsesFactoryModule
    use postureTaskClass
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
    real(wp), dimension(:), allocatable :: t, gII_V_mV, IaIn_V_mV,  IbIn_V_mV, IbFR
    real(wp) :: tic, toc
    type(gpf) :: gp
    real(wp) :: FR
    integer :: GammaOrder 
    character(len = 80) :: pool, muscle, group
    character(len = 95) :: filename = 'confPosturalControl.rmto'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools   
    type(AfferentPool), dimension(:), allocatable:: afferentPools    
    type(postureTask) :: body
    real(wp) :: angle
    real(wp) , dimension(200):: Amp, phase
    logical :: continueFlag
    real(wp) :: dynGamma, statGamma
    integer :: subject, trial
    character*198 buf
    integer :: location
    character(len = 80) :: paramTag


    call init_random_seed()
    print *, filename
    call mkl_get_version_string(buf)
    write(*,'(a)') buf

    filename = 'confPosturalControl.rmto'
    conf = Configuration(filename)

    tf = conf%simDuration_ms
    dt = conf%timeStep_ms
    timeLength = int(tf/dt)

    allocate(t(timeLength))
    
    t = [(dt*(i-1), i=1, timeLength)]
    allocate(gII_V_mV(timeLength))
    allocate(IaIn_V_mV(timeLength))
    allocate(IbIn_V_mV(timeLength))
    allocate(IbFR(timeLength))

    FR = 77
    GammaOrder = 25
    
    filename = 'GammaFallTime3.txt'
    open(1, file=filename, status = 'replace')

    if (allocated(afferentPools)) deallocate(afferentPools)
    if (allocated(motorUnitPools)) deallocate(motorUnitPools)
    if (allocated(neuralTractPools)) deallocate(neuralTractPools)
    if (allocated(synapticNoisePools)) deallocate(synapticNoisePools)
    if (allocated(interneuronPools)) deallocate(interneuronPools)

    subject = 1
    do while (subject <= 10)
        
        allocate(afferentPools(12))

        paramTag = 'Con:Ia-SOL>SOL-S@dendrite|excitatory'
        call conf%showConfigurationParameter(paramTag)

        paramTag = 'Con:Ia-LG>SOL-S@dendrite|excitatory'
        call conf%showConfigurationParameter(paramTag)

        pool = 'Ia'
        muscle = 'SOL'
        afferentPools(1) = AfferentPool(conf, pool, muscle)

        pool = 'Ia'
        muscle = 'MG'
        afferentPools(2) = AfferentPool(conf, pool, muscle)

        pool = 'Ia'
        muscle = 'LG'
        afferentPools(3) = AfferentPool(conf, pool, muscle)

        pool = 'Ia'
        muscle = 'TA'
        afferentPools(4) = AfferentPool(conf, pool, muscle)

        pool = 'II'
        muscle = 'SOL'
        afferentPools(5) = AfferentPool(conf, pool, muscle)

        pool = 'II'
        muscle = 'MG'
        afferentPools(6) = AfferentPool(conf, pool, muscle)

        pool = 'II'
        muscle = 'LG'
        afferentPools(7) = AfferentPool(conf, pool, muscle)

        pool = 'II'
        muscle = 'TA'
        afferentPools(8) = AfferentPool(conf, pool, muscle)

        pool = 'Ib'
        muscle = 'SOL'
        afferentPools(9) = AfferentPool(conf, pool, muscle)

        pool = 'Ib'
        muscle = 'MG'
        afferentPools(10) = AfferentPool(conf, pool, muscle)

        pool = 'Ib'
        muscle = 'LG'
        afferentPools(11) = AfferentPool(conf, pool, muscle)

        pool = 'Ib'
        muscle = 'TA'
        afferentPools(12) = AfferentPool(conf, pool, muscle)

        allocate(neuralTractPools(1))
        pool = 'CMExt'
        neuralTractPools(1) = NeuralTract(conf, pool)
        
        allocate(motorUnitPools(4))
        pool = 'SOL'
        motorUnitPools(1) = MotorUnitPool(conf, pool)

        pool = 'MG'
        motorUnitPools(2) = MotorUnitPool(conf, pool)    

        pool = 'LG'
        motorUnitPools(3) = MotorUnitPool(conf, pool)    

        pool = 'TA'
        motorUnitPools(4) = MotorUnitPool(conf, pool)    

        body = postureTask(conf, motorUnitPools)
        
        allocate(interneuronPools(3))
        pool = 'gII'
        group = 'ext'
        interneuronPools(1) = InterneuronPool(conf, pool, group)

        pool = 'IaIn'
        group = 'ext'
        interneuronPools(2) = InterneuronPool(conf, pool, group)

        pool = 'IbIn'
        group = 'ext'
        interneuronPools(3) = InterneuronPool(conf, pool, group)

        synapticNoisePools = synapseFactory(conf, neuralTractPools, &
                                            motorUnitPools, &
                                            interneuronPools, &
                                            afferentPools)
        

        ! print '(a, F15.6)', 'length', interneuronPools(1)%unit(1)%Compartments(1)%length_mum
        ! print '(a, F15.6)', 'gleak', interneuronPools(1)%unit(1)%Compartments(1)%gLeak_muS
        ! print '(a, F15.6)', 'diameter', interneuronPools(1)%unit(1)%Compartments(1)%diameter_mum
        ! print '(a, F15.6)', 'capacitance', interneuronPools(1)%unit(1)%Compartments(1)%capacitance_nF
        
        ! print '(a, F15.6)', 'length', interneuronPools(2)%unit(1)%Compartments(1)%length_mum
        ! print '(a, F15.6)', 'gleak', interneuronPools(2)%unit(1)%Compartments(1)%gLeak_muS
        ! print '(a, F15.6)', 'diameter', interneuronPools(2)%unit(1)%Compartments(1)%diameter_mum
        ! print '(a, F15.6)', 'capacitance', interneuronPools(2)%unit(1)%Compartments(1)%capacitance_nF

        ! print '(a, F15.6)', 'length', interneuronPools(3)%unit(1)%Compartments(1)%length_mum
        ! print '(a, F15.6)', 'gleak', interneuronPools(3)%unit(1)%Compartments(1)%gLeak_muS
        ! print '(a, F15.6)', 'diameter', interneuronPools(3)%unit(1)%Compartments(1)%diameter_mum
        ! print '(a, F15.6)', 'capacitance', interneuronPools(3)%unit(1)%Compartments(1)%capacitance_nF
        call random_number(dynGamma)
        dynGamma = 30.5!29.5 + 8*dynGamma
        call random_number(statGamma)
        statGamma = 42.0!34.9 + 8*statGamma
        print '(A, F15.6)', 'Dynamic Gamma = ',  dynGamma
        print '(A, F15.6)', 'Static Gamma = ', statGamma
        continueFlag = .true.
        trial = 1
        do while (continueFlag .and. trial <= 2)
            ! print *, 'IaIn synpases ' , size(interneuronPools(2)%unit(1)%Compartments(1)%SynapsesIn(1)%gmax_muS),&
            !          size(interneuronPools(2)%unit(50)%Compartments(1)%SynapsesIn(1)%gmax_muS),&
            !          size(interneuronPools(2)%unit(100)%Compartments(1)%SynapsesIn(1)%gmax_muS) 
            ! print *, 'SOL inhibitory synapses ', size(motorUnitPools(1)%unit(1)%Compartments(2)%SynapsesIn(2)%gmax_muS),&
            !          size(motorUnitPools(1)%unit(200)%Compartments(2)%SynapsesIn(2)%gmax_muS),&
            !          size(motorUnitPools(1)%unit(400)%Compartments(2)%SynapsesIn(2)%gmax_muS)
            ! print *, 'MG excitatory synapses ', size(motorUnitPools(2)%unit(1)%Compartments(1)%SynapsesIn(1)%gmax_muS),&
            !          size(motorUnitPools(2)%unit(200)%Compartments(1)%SynapsesIn(1)%gmax_muS),&
            !          size(motorUnitPools(2)%unit(400)%Compartments(1)%SynapsesIn(1)%gmax_muS)
            ! print *, 'SOL excitatory synapses ',size(motorUnitPools(1)%unit(1)%Compartments(1)%SynapsesIn(1)%gmax_muS),&
            !          size(motorUnitPools(1)%unit(200)%Compartments(1)%SynapsesIn(1)%gmax_muS),&
            !          size(motorUnitPools(1)%unit(400)%Compartments(1)%SynapsesIn(1)%gmax_muS)        
            ! print *, 'Ia MG making synapses ', size(afferentPools(2)%unit(1)%transmitSpikesThroughSynapses),&
            !         size(afferentPools(2)%unit(50)%transmitSpikesThroughSynapses),&
            !         size(afferentPools(2)%unit(100)%transmitSpikesThroughSynapses) 
            continueFlag = .false.
            tic =  dsecnd()
            i = 1            
            do while (i <= size(t))
                do j = 1, size(neuralTractPools)
                    ! if (t(i) > 1000 .and. t(i) < 4000) then 
                    !     FR = 26.3 + 10.0 * (body%ankleTorque_Nm(nint(t(i)/conf%timeStep_ms)))
                    ! end if
                    ! FR = 62
                    !if (t(i) > 1000 .and. t(i) < 4001) print *, 'FR = ', FR
                    call neuralTractPools(j)%atualizePool(t(i), FR, GammaOrder)
                end do
                do j = 1, 3
                    call motorUnitPools(j)%atualizeMotorUnitPool(t(i), dynGamma, statGamma)
                end do
                call motorUnitPools(4)%atualizeMotorUnitPool(t(i), 0.0_wp, 0.0_wp)
                call interneuronPools(1)%atualizeInterneuronPool(t(i))
                call interneuronPools(2)%atualizeInterneuronPool(t(i))
                call interneuronPools(3)%atualizeInterneuronPool(t(i))
                gII_V_mV(i) = interneuronPools(1)%v_mV(1)
                !IaIn_V_mV(i) = interneuronPools(2)%v_mV(1)
                IbIn_V_mV(i) = interneuronPools(3)%v_mV(1)
                call body%atualizeBody(t(i))
                do j = 1, 4
                    call afferentPools(j)%atualizeAfferentPool(t(i), motorUnitPools(j)%spindle%IaFR_Hz)
                    call afferentPools(j+4)%atualizeAfferentPool(t(i), motorUnitPools(j)%spindle%IIFR_Hz)
                    call afferentPools(j+8)%atualizeAfferentPool(t(i), &
                                           motorUnitPools(j)%GTO%IbFR_Hz)
                end do
                IbFR(i) = motorUnitPools(2)%GTO%IbFR_Hz
                if (body%ankleAngle_rad(i) > 20*pi/180 .or. body%ankleAngle_rad(i) < -15*pi/180) then
                    print '(A, 1X, F15.6, 1X, A)', 'Body fell after ', t(i)/1000.0, ' seconds'
                    continueFlag = .true.
                    write(1, '(F15.6, 1X, F15.6,1X,F15.6, 1X, I2, 1X, I2)') ([t(i), dynGamma, statGamma]), ([subject, trial])
                    i = size(t)
                    call random_number(dynGamma)
                    dynGamma = 30.5!29.5 + 8*dynGamma
                    call random_number(statGamma)
                    statGamma = 42.0!34.9 + 8*statGamma
                    print '(A, F15.6)', 'Dynamic Gamma = ', dynGamma
                    print '(A, F15.6)', 'Static Gamma = ', statGamma
                end if
                i = i + 1
            end do
            toc =  dsecnd()

            print '(F15.6, A)', toc - tic, ' seconds'
            
            ! call neuralTractPools(1)%listSpikes()
            call motorUnitPools(2)%listSpikes()
            call afferentPools(2)%listSpikes()
            call afferentPools(4)%listSpikes()
            call afferentPools(6)%listSpikes()
            call afferentPools(10)%listSpikes()
            ! !call motorUnitPools(1)%getMotorUnitPoolEMG()
            
            ! call gp%title('MN spike instants at the soma')
            ! call gp%xlabel('t (s))')
            ! call gp%ylabel('Motoneuron index')
            ! call gp%plot(motorUnitPools(2)%poolSomaSpikes(:,1), &
            ! motorUnitPools(2)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

            call gp%title('MG AF Ia spike instants at the terminal')
            call gp%xlabel('t (s))')
            call gp%ylabel('Afferent index')
            call gp%plot(afferentPools(2)%poolTerminalSpikes(:,1), &
            afferentPools(2)%poolTerminalSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

            ! call gp%title('AF Ib MG spike instants at the terminal')
            ! call gp%xlabel('t (s))')
            ! call gp%ylabel('Afferent index')
            ! call gp%plot(afferentPools(10)%poolTerminalSpikes(:,1), &
            ! afferentPools(10)%poolTerminalSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

            ! call gp%title('TA AF Ia spike instants at the terminal')
            ! call gp%xlabel('t (s))')
            ! call gp%ylabel('Afferent index')
            ! call gp%plot(afferentPools(4)%poolTerminalSpikes(:,1), &
            ! afferentPools(4)%poolTerminalSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

            ! call gp%title('AF II spike instants at the soma')
            ! call gp%xlabel('t (s))')
            ! call gp%ylabel('Afferent II index')
            ! call gp%plot(afferentPools(6)%poolTerminalSpikes(:,1), &
            ! afferentPools(6)%poolTerminalSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')
            
            ! call gp%title('Membrane Potential gII')
            ! call gp%xlabel('t (ms))')
            ! call gp%ylabel('V (mV)')
            ! call gp%plot(t, gII_V_mV, 'with line lw 2 lc rgb "#0008B0"')

            ! call gp%title('Membrane Potential IaIn')
            ! call gp%xlabel('t (ms))')
            ! call gp%ylabel('V (mV)')
            ! call gp%plot(t, IaIn_V_mV, 'with line lw 2 lc rgb "#0008B0"')

            ! call gp%title('FR Ib mg')
            ! call gp%xlabel('t (ms))')
            ! call gp%ylabel('FR (Hz)')
            ! call gp%plot(t, IbFR, 'with line lw 2 lc rgb "#0008B0"')

            ! call gp%title('Membrane Potential IbIn')
            ! call gp%xlabel('t (ms))')
            ! call gp%ylabel('V (mV)')
            ! call gp%plot(t, IbIn_V_mV, 'with line lw 2 lc rgb "#0008B0"')

            ! call gp%title('MN spike instants at the terminal')
            ! call gp%xlabel('t (s))')
            ! call gp%ylabel('Motoneuron index')
            ! call gp%plot(motorUnitPools(2)%poolTerminalSpikes(:,1), &
            ! motorUnitPools(2)%poolTerminalSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

            ! call gp%title('SOL force')
            ! call gp%xlabel('t (ms))')
            ! call gp%ylabel('Force (N)')
            ! call gp%plot(t, motorUnitPools(1)%HillMuscle%force, 'with line lw 2 lc rgb "#0008B0"')

            ! call gp%title('SOL LENGTH')
            ! call gp%xlabel('t (ms))')
            ! call gp%ylabel('Force (N)')
            ! call gp%plot(t, motorUnitPools(1)%HillMuscle%length_m, 'with line lw 2 lc rgb "#0008B0"')

            ! call gp%title('MG force')
            ! call gp%xlabel('t (ms))')
            ! call gp%ylabel('torque (N.m)')
            ! call gp%plot(t, motorUnitPools(2)%HillMuscle%tendonForce_N, 'with line lw 2 lc rgb "#0008B0"')

            ! call gp%title('MG LENGTH')
            ! call gp%xlabel('t (ms))')
            ! call gp%ylabel('Force (N)')
            ! call gp%plot(t, motorUnitPools(2)%HillMuscle%length_m, 'with line lw 2 lc rgb "#0008B0"')

            ! call gp%title('LG LENGTH')
            ! call gp%xlabel('t (ms))')
            ! call gp%ylabel('Force (N)')
            ! call gp%plot(t, motorUnitPools(3)%HillMuscle%length_m, 'with line lw 2 lc rgb "#0008B0"')

            ! call gp%title('TA LENGTH')
            ! call gp%xlabel('t (ms))')
            ! call gp%ylabel('Force (N)')
            ! call gp%plot(t, motorUnitPools(4)%HillMuscle%length_m, 'with line lw 2 lc rgb "#0008B0"')

            ! call gp%title('Ankle torque')
            ! call gp%xlabel('t (ms))')
            ! call gp%ylabel('torque (N.m)')
            ! call gp%plot(t, body%ankleTorque_Nm, 'with line lw 2 lc rgb "#0008B0"')

            call gp%title('Ankle angle')
            call gp%xlabel('t (ms)')
            call gp%ylabel('angle (degree)')
            call gp%plot(t, body%ankleAngle_rad*180.0/pi, 'with line lw 2 lc rgb "#0008B0"')

            do j = 1, size(neuralTractPools)
                call neuralTractPools(j)%reset()
            end do

            do j = 1, 4
                call motorUnitPools(j)%reset()
            end do            

            do j = 1, 12
                call afferentPools(j)%reset()
            end do     

            call interneuronPools(1)%reset()       
            call interneuronPools(2)%reset()    
            call interneuronPools(3)%reset()       
            call body%reset()
            IbFR(:) = 0.0
            gII_V_mV(:) = 0.0
            IaIn_V_mV(:) = 0.0
            IbIn_V_mV(:) = 0.0
            i = 1
            trial = trial + 1
        end do
        
        write(1, '(F15.6, 1X, F15.6, 1X, F15.6, 1X, I2, 1X, I2)') ([t(i), dynGamma, statGamma]), ([subject, trial])
        subject = subject + 1
        deallocate(afferentPools)
        deallocate(motorUnitPools)
        deallocate(neuralTractPools)
        deallocate(synapticNoisePools)
        deallocate(interneuronPools)
    end do
    close(1)
    
end program PosturalControl