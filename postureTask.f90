! '''
!     Neuromuscular simulator in Fortran.
!     Copyright (C) 2019  Renato Naville Watanabe

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

module postureTaskClass
    use ConfigurationClass
    use MotorUnitPoolClass
    use MusclePointerClass
    implicit none
    private
    integer, parameter :: wp = kind(1.0d0)
    real(wp), parameter :: pi = 4 * atan(1.0_wp)    
    real(wp), parameter :: g = 9.80655
    public :: postureTask
    

    type postureTask
        type(Configuration), pointer :: conf
        type(MusclePointer), dimension(:), allocatable :: muscles
        real(wp), dimension(:), allocatable :: ankleAngle_rad, ankleTorque_Nm, ankleOmega_rad_s
        real(wp) :: momentOfInertia, mass, height
    
        contains
            procedure :: atualizeBody
            procedure :: computeTorque
            procedure :: reset

    end type postureTask

    interface postureTask   
        module procedure init_postureTask
    end interface postureTask

    contains

        type(postureTask) function   init_postureTask(conf, pools)
            class(Configuration), intent(in), target :: conf
            class(MotorUnitPool), intent(in) :: pools(:)
            integer :: numberOfPools, i, timeLength

            init_postureTask%conf => conf
            
            numberOfPools = size(pools)
            allocate(init_postureTask%muscles(numberOfPools))

            do i = 1, numberOfPools
                if (pools(i)%pool == 'SOL' .or. pools(i)%pool == 'MG' .or. pools(i)%pool == 'LG' .or. pools(i)%pool == 'TA') then
                    init_postureTask%muscles(i) = MusclePointer()
                    call init_postureTask%muscles(i)%assignMuscle(pools(i))
                end if
            end do
            ! ##
            timeLength = nint(conf%simDuration_ms/conf%timeStep_ms)
            allocate(init_postureTask%ankleAngle_rad(timeLength))
            init_postureTask%ankleAngle_rad(:) = 5.0*pi/180.0
            allocate(init_postureTask%ankleTorque_Nm(timeLength))
            init_postureTask%ankleTorque_Nm(:) = 0.0
            allocate(init_postureTask%ankleOmega_rad_s(timeLength))
            init_postureTask%ankleOmega_rad_s(:) = 0.0

            init_postureTask%mass = 60.0
            init_postureTask%height = 0.85

            init_postureTask%momentOfInertia = 4.0*(init_postureTask%mass*init_postureTask%height**2)/3.0

            print '(A)', 'Posture Task built'

        end function 

        subroutine atualizeBody(self, t)
            ! '''
            ! Update the ankle joint.
            ! Atualizes the musculotendon length and the moment-arm of each muscle.
            ! Updates the angle and angular velocity by numerically solving 
            ! the differential equations with the Euler method.
            ! 
            ! Inputs:
            !
            ! * t: real 
            ! 
            
            class(postureTask), intent(inout) :: self
            real(wp), intent(in) ::t
            integer :: i, timeIndex
            real(wp) :: dthetadt, domegadt, angle


            timeIndex = nint(t/self%conf%timeStep_ms)+1
            
            angle = self%ankleAngle_rad(timeIndex)*180/pi
            if (self%muscles(1)%muscle%hillModel == 'No') then
                do i = 1, size(self%muscles)
                    
                    call self%muscles(i)%muscle%NoHillMuscle%atualizeMusculoTendonLength(angle)
                    call self%muscles(i)%muscle%NoHillMuscle%atualizeMomentArm(angle)
                end do
            else 
                do i = 1, size(self%muscles)
                    call self%muscles(i)%muscle%HillMuscle%atualizeMusculoTendonLength(angle)
                    call self%muscles(i)%muscle%HillMuscle%atualizeMomentArm(angle)
                end do
            end if
            
            call self%computeTorque(t)

            if (t > 1000.0) then
                dthetadt = self%ankleOmega_rad_s(timeIndex)
                domegadt = self%ankleTorque_Nm(timeIndex)/self%momentOfInertia
                self%ankleOmega_rad_s(timeIndex + 1) = self%ankleOmega_rad_s(timeIndex) + &
                                                       self%conf%timeStep_ms*domegadt/1000.0
                
                self%ankleAngle_rad(timeIndex + 1) = self%ankleAngle_rad(timeIndex) + &
                                                     self%conf%timeStep_ms*dthetadt/1000.0
            else                
                self%ankleOmega_rad_s(timeIndex + 1) = 0.0
                self%ankleAngle_rad(timeIndex + 1) = 5.0*pi/180.0
            end if 
            
        end subroutine


        subroutine computeTorque(self, t)
            ! '''
            ! '''
            class(postureTask), intent(inout) :: self
            real(wp), intent(in) ::t
            real(wp) :: muscularTorque, velocity, acceleration
            integer :: i
            integer :: timeIndex

            timeIndex = nint(t/self%conf%timeStep_ms) + 1
                
            muscularTorque = 0.0

            if (self%muscles(1)%muscle%hillModel == 'No') then
                do i = 1, size(self%muscles)
                    muscularTorque = muscularTorque + &
                    self%muscles(i)%muscle%NoHillMuscle%force(timeIndex) * &
                    self%muscles(i)%muscle%NoHillMuscle%momentArm_m(timeIndex)
                end do
            else
                do i = 1, size(self%muscles)
                    muscularTorque = muscularTorque + &
                    self%muscles(i)%muscle%HillMuscle%force(timeIndex) * &
                    self%muscles(i)%muscle%HillMuscle%momentArm_m(timeIndex)
                end do
            end if
            
            if (mod(t, 1000.0_wp) == 0) then
                print *, 'Sol Force ', self%muscles(1)%muscle%HillMuscle%force(timeIndex)
                print *, 'MG Force ', self%muscles(2)%muscle%HillMuscle%force(timeIndex)
                print *, 'LG Force ', self%muscles(3)%muscle%HillMuscle%force(timeIndex)
                print *, 'TA Force ', self%muscles(4)%muscle%HillMuscle%force(timeIndex)
                print *, 'Sol Activation ', self%muscles(1)%muscle%HillMuscle%activationTypeI(timeIndex)
                print *, 'MG Activation ', self%muscles(2)%muscle%HillMuscle%activationTypeI(timeIndex)
                print *, 'LG Activation ', self%muscles(3)%muscle%HillMuscle%activationTypeI(timeIndex)
                print *, 'TA Activation ', self%muscles(4)%muscle%HillMuscle%activationTypeI(timeIndex)
                print *, 'Sol length ', self%muscles(1)%muscle%HillMuscle%length_m(timeIndex)
                print *, 'MG length ', self%muscles(2)%muscle%HillMuscle%length_m(timeIndex)
                print *, 'LG length ', self%muscles(3)%muscle%HillMuscle%length_m(timeIndex)
                print *, 'TA length ', self%muscles(4)%muscle%HillMuscle%length_m(timeIndex)
                print *, 'muscle ', muscularTorque
                print *, 'passive', -0.65*self%mass*g*self%height*(self%ankleAngle_rad(timeIndex)-0*pi/180)
                print *, 'gravity', self%mass*g*self%height*sin(self%ankleAngle_rad(timeIndex))
                print *, 'viscosity', - 5.81*self%ankleOmega_rad_s(timeIndex)
            end if

            self%ankleTorque_Nm(timeIndex) = muscularTorque - 5.81*self%ankleOmega_rad_s(timeIndex) - &
                                             0.65*self%mass*g*self%height*(self%ankleAngle_rad(timeIndex)-0*pi/180) + &
                                             self%mass*g*self%height*sin(self%ankleAngle_rad(timeIndex))
        
        end subroutine

        subroutine reset(self)
            ! '''
            ! '''
            class(postureTask), intent(inout) :: self
            
            self%ankleAngle_rad(:) = 5.0*pi/180.0
            self%ankleTorque_Nm(:) = 0.0
            self%ankleOmega_rad_s(:) = 0.0
        end subroutine
end module postureTaskClass        
    