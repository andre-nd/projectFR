! '''
!     Neuromuscular simulator in Fortran.
!     Copyright (C) 2019  Renato Naville Watanabe
!                         Pablo Alejandro    
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     any later version.

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.

!     Contact: renato.watanabe@ufabc.edu.br
! '''

program AntidromicStimulationofMNandRC
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
    use AfferentPoolClass
    use SynapsesFactoryModule
    
    implicit none 
    !integer, parameter :: wp = kind(1.0d0)
    type(Configuration) :: conf
    real(wp), parameter :: pi = 4 * atan(1.0_wp)    
    real(wp) :: dt
    real(wp) :: tf
    integer :: timeLength
    integer :: i, j
    real(wp), dimension(:), allocatable :: t, MNv_mV, RCv_mV
    real(wp) :: tic, toc
    type(gpf) :: gp
    real(wp) :: FR
    integer :: GammaOrder 
    character(len = 80) :: pool, group
    character(len = 80) :: filename = 'confAntidromicStimulationofMNandRC.rmto'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools     
    type(AfferentPool), dimension(:), allocatable:: afferentPools     

    call init_random_seed()

    conf = Configuration(filename)
    allocate(neuralTractPools(0))
    pool = 'SOL'
    allocate(motorUnitPools(1))
    motorUnitPools(1) = MotorUnitPool(conf, pool)    
    pool = 'RC'
    group = 'ext'    
    allocate(interneuronPools(1))
    interneuronPools(1) = InterneuronPool(conf, pool, group)
    allocate(afferentPools(0))
    synapticNoisePools = synapseFactory(conf, neuralTractPools, &
                                        motorUnitPools, &
                                        interneuronPools, &
                                        afferentPools)
    
    conf = Configuration(filename)
    
    tf = conf%simDuration_ms
    dt = conf%timeStep_ms
    timeLength = nint(tf/dt)
    
    allocate(t(timeLength))
    allocate(MNv_mV(timeLength))
    allocate(RCv_mV(timeLength))

    
      
    t = [(dt*(i-1), i=1, timeLength)]
    
    do j = 1, size(motorUnitPools)
        call motorUnitPools(j)%reset()
    end do

    do j = 1, size(interneuronPools)
        call interneuronPools(j)%reset()
    end do

    do j = 1, size(synapticNoisePools)
        call synapticNoisePools(j)%reset()
    end do

      

     

    call cpu_time(tic)
    do i = 1, size(t)        
        do j = 1, size(synapticNoisePools)
            call synapticNoisePools(j)%atualizePool(t(i))
        end do
        do j = 1, size(interneuronPools)
            call interneuronPools(j)%atualizeInterneuronPool(t(i))
            RCv_mV(i) = interneuronPools(j)%v_mV(1)        
        end do
        do j = 1, size(motorUnitPools)
            call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 32.0_wp, 32.0_wp)
            MNv_mV(i) = motorUnitPools(j)%v_mV(2)      
        end do
    end do

    call cpu_time(toc)

    print '(F15.6, A)', toc - tic, ' seconds'
    
    call motorUnitPools(1)%listSpikes()
    call interneuronPools(1)%listSpikes()
    
    call gp%title('PTN stimulus')
    call gp%xlabel('t (ms))')
    call gp%ylabel('Stimulus (mA)')
    call gp%plot(t, motorUnitPools(1)%unit(1)%nerveStimulus_mA, 'with line lw 2 lc rgb "#0008B0"')  
    
    call gp%title('Membrane potential of the soma of the MN #1')
    call gp%xlabel('t (ms))')
    call gp%ylabel('Descending command index')
    call gp%plot(t, MNv_mV, 'with line lw 2 lc rgb "#0008B0"')  

    call gp%title('Membrane potential of the soma of the RC #1')
    call gp%xlabel('t (ms))')
    call gp%ylabel('Descending command index')
    call gp%plot(t, RCv_mV, 'with line lw 2 lc rgb "#0008B0"')  


    call gp%title('MN spike instants at the soma')
    call gp%xlabel('t (s))')
    call gp%ylabel('Motoneuron index')
    call gp%plot(motorUnitPools(1)%poolSomaSpikes(:,1), &
    motorUnitPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

    call gp%title('RC spike instants at the soma')
    call gp%xlabel('t (s))')
    call gp%ylabel('Motoneuron index')
    call gp%plot(interneuronPools(1)%poolSomaSpikes(:,1), &
    interneuronPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')


    call gp%title('Muscle force')
    call gp%xlabel('t (ms))')
    call gp%ylabel('Force (N)')
    call gp%plot(t, motorUnitPools(1)%NoHillMuscle%force, 'with line lw 2 lc rgb "#0008B0"')
    
    
end program AntidromicStimulationofMNandRC