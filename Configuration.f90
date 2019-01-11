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

! '''
! \mainpage ReMoto in Fortran

! This program is a neuronal simulation system, intended for studying spinal cord neuronal 
! networks responsible for muscle control. These networks are affected by descending drive, 
! afferent drive, and electrical nerve stimulation. The simulator may be used to investigate
! phenomena at several levels of organization, e.g., at the neuronal membrane level or at 
! the whole muscle behavior level (e.g., muscle force generation). This versatility is due 
! to the fact that each element (neurons, synapses, muscle fibers) has its own specific 
! mathematical model, usually involving the action of voltage- or neurotransmitter-dependent
! ionic channels. The simulator should be helpful in activities such as interpretation of
! results obtained from neurophysiological experiments in humans or mammals, proposal of 
! hypothesis or testing models or theories on neuronal dynamics or neuronal network processing,
! validation of experimental protocols, and teaching neurophysiology.

! The elements that take part in the system belong to the following classes: motoneurons, 
! muscle fibers (electrical activity and force generation), Renshaw cells, Ia inhibitory 
! interneurons, Ib inhibitory interneurons, Ia and Ib afferents. The neurons are interconnected
! by chemical synapses, which can be exhibit depression or facilitation.

! The system simulates the following nuclei involved in flexion and extension of the human or
! cat ankle: Medial Gastrocnemius (MG), Lateral Gastrocnemius (LG), Soleus (SOL), and Tibialis
! Anterior (TA).

! A web-based version can be found in [remoto.leb.usp.br](http://remoto.leb.usp.br/remoto/index.html).
! The version to which this documentation  refers is from a Python program that can be found in
! [github.com/rnwatanabe/projectFR](https://github.com/rnwatanabe/projectFR).

! '''

module ConfigurationClass
    ! '''
    ! Class that builds an object of Configuration, based on a configuration file.
    ! '''
    use CharacterArrayClass
    use CharacterMatrixClass
    implicit none
    private
    integer, parameter :: wp = kind( 1.0d0 )
    public :: Configuration

    type Configuration
        character(len = 80) :: filename
        real(wp) :: timeStep_ms, simDuration_ms, skinThickness_mm
        real(wp) :: timeStepByTwo_ms, timeStepBySix_ms
        real(wp) :: EMGAttenuation_mm1, EMGWidening_mm1, EMGNoiseEMG
        character(len = 80) :: MUParameterDistribution
        type(CharacterMatrix) :: confMatrix
        
        contains
            procedure :: parameterSet
            procedure :: determineSynapses
            procedure :: changeConfigurationParameter
            procedure :: showConfigurationParameter

    end type Configuration

    interface Configuration
        module procedure init_Configuration
    end interface

    contains

            type(Configuration) function init_Configuration(filename)
                ! '''
                ! Constructor.
                
                ! Builds the Configuration object. A Configuration object is responsible to set the variables
                ! that are used in the whole system, such as timeStep and simDuration.
                
                ! - Inputs:
                !     + **filename**: name of the file with the parameter values. The extension  of the file should be .rmto.
                
                ! '''
                character(*), intent(in) :: filename
                integer :: ierr, il, j, stop1, i
                character(len = 80) :: line
                character(len = 80) :: param1, param2, param3
                type(CharacterArray) :: newLine

                init_Configuration%filename = filename
                init_Configuration%confMatrix = CharacterMatrix()
                open(1,file = init_Configuration%filename, status='old',iostat=ierr)
                
                do while (ierr.eq.0)
                    read(1, '(A)', iostat=ierr) line
                    il=len_trim(line)
                    j = 1
                    do i = 1, il
                        if (line(i:i) == ',') then
                            if (j.eq.1) then 
                                param1 = line(1:i-1)
                                j = j + 1
                                stop1 = i
                            else if (j.eq.2) then 
                                param2 = line(stop1+1:i-1)
                                param3 = line(i+1:il)
                            end if 
                        end if            
                    end do
                    !## Time step of the numerical solution of the differential equation.
                    if (j == 2) then
                        if (param1=='timeStep') then
                            read(param2(1:len_trim(param2)), *)init_Configuration%timeStep_ms
                        end if
                        !## Total length of the simulation in ms.    
                        if (param1=='simDuration') then 
                            read(param2(1:len_trim(param2)), *)init_Configuration%simDuration_ms
                        end if
                        !## skin thickness, in mm.
                        if (param1=='skinThickness') then 
                            read(param2(1:len_trim(param2)), *)init_Configuration%skinThickness_mm
                        end if
                        !## EMG attenuation factor, in 1/mm.  
                        if (param1=='EMGAttenuationFactor') then 
                            read(param2(1:len_trim(param2)), *)init_Configuration%EMGAttenuation_mm1
                        end if
                        !## EMG widening factor, in 1/mm.
                        if (param1=='EMGWideningFactor') then
                            read(param2(1:len_trim(param2)), *)init_Configuration%EMGWidening_mm1
                        end if
                        !## EMG widening factor.
                        if (param1=='EMGNoiseEMG') then
                            read(param2(1:len_trim(param2)), *)init_Configuration%EMGNoiseEMG
                        end if
                        !## Distribution of the parameters along the motor units.                
                        if (param1=='MUParameterDistribution') then 
                            init_Configuration%MUParameterDistribution = param2
                        end if
                        !## The variable  timeStep divided by two, for computational efficiency.
                        init_Configuration%timeStepByTwo_ms = init_Configuration%timeStep_ms / 2.0 
                        !## The variable  timeStep divided by six, for computational efficiency.
                        init_Configuration%timeStepBySix_ms = init_Configuration%timeStep_ms / 6.0
                        newLine = CharacterArray()
                        call newLine%AddToList(param1)
                        call newLine%AddToList(param2)
                        call newLine%AddToList(param3)
                        call init_Configuration%confMatrix%append(newLine)
                    end if
                end do
                
                close(unit = 1)
                
            
            end function init_Configuration

            character(len=80) function parameterSet(self,  paramTag, pool, index) result(requestedParamater)
                ! '''
                ! Function that returns the value of wished parameter specified in the paramTag variable.
                ! In the case of min/max parameters, the value returned is the specific to the index of the unit that called the
                ! function. 


                ! - Inputs: 

                !     + **paramTag**: string with the name of the wished parameter as in the first column of the rmto file.

                !     + **pool**: pool from which the unit that will receive the parameter value belongs. For example SOL. 
                !     It is used only in the parameters that have a range.

                !     + **index**: index of the unit. It is is an integer.

                ! - Outputs:
                    
                !     + required parameter value
                ! '''
                class(Configuration), intent(inout) :: self
                character(*), intent(in) :: pool
                integer, intent(in) :: index
                integer :: ierr, il, j, stop1, i, k
                character(len = 80) :: line
                character(len = 80) :: param1, param2, param3
                real(wp) :: param2Real, param3Real
                integer :: MUnumber_S, MUnumber_FR, MUnumber_FF, Nnumber
                real(wp), dimension(:), allocatable :: paramVec_S, paramVec_FR, paramVec_FF, paramVec
                character(len=50), intent(in) ::paramTag
                logical :: distribute, wholePool
                real(wp), dimension(:), allocatable :: indexUnits
                
                distribute = .true.
                wholePool = .false.
                MUnumber_S = 0
                MUnumber_FR = 0
                MUnumber_FF = 0
                Nnumber = 0              
                
                do k = 1, size(self%confMatrix%item)
                    param1 = self%confMatrix%item(k)%item(1)%string
                    param2 = self%confMatrix%item(k)%item(2)%string
                    param3 = self%confMatrix%item(k)%item(3)%string
                    
                    
                    if (pool=='SOL'.or.pool=='MG'.or.pool=='LG'.or.pool=='TA') then
                        if (param1.eq.('MUnumber_' // trim(pool) // '-S')) then
                            read(param2(1:len_trim(param2)), *)MUnumber_S
                        else if (param1.eq.('MUnumber_' // trim(pool) // '-FR')) then
                            read(param2(1:len_trim(param2)), *)MUnumber_FR
                        else if (param1.eq.('MUnumber_' // trim(pool) // '-FF')) then
                            read(param2(1:len_trim(param2)), *)MUnumber_FF
                        end if
                        Nnumber = MUnumber_S + MUnumber_FR + MUnumber_FF 
                    else 
                        if (trim(param1).eq.('Number_' // trim(pool))) then 
                            read(param2(1:len_trim(param2)), *)Nnumber
                        end if
                    end if
                end do
                
                allocate(paramVec(Nnumber)) 
                if (allocated(paramVec_S)) deallocate(paramVec_S)
                if (allocated(paramVec_FR)) deallocate(paramVec_FR)
                if (allocated(paramVec_FF)) deallocate(paramVec_FF)
                
                do k = 1, size(self%confMatrix%item)
                    param1 = self%confMatrix%item(k)%item(1)%string
                    param2 = self%confMatrix%item(k)%item(2)%string
                    param3 = self%confMatrix%item(k)%item(3)%string
                    
                    
                    if (trim(param1)==trim(paramTag)) then 
                        requestedParamater = param2
                        distribute = .false.
                    else if (trim(self%MUParameterDistribution)=='linear') then     
                        if (trim(param1).eq.(trim(paramTag) // ':' // trim(pool) // '-S')) then
                            if (MUnumber_S>0) allocate(paramVec_S(MUnumber_S))
                            read(param2,*)param2Real
                            read(param3,*)param3Real
                            paramVec_S = [((param3Real-param2Real)/(MUnumber_S+1)*(i-1)+param2Real, i=1, MUnumber_S)]                                    
                            paramVec(1:MUnumber_S) = paramVec_S                            
                        else if (trim(param1).eq.(trim(paramTag) // ':' // trim(pool) // '-FR')) then
                            if (MUnumber_FR>0) allocate(paramVec_FR(MUnumber_FR))  
                            read(param2,*)param2Real
                            read(param3,*)param3Real
                            paramVec_FR = [((param3Real-param2Real)/(MUnumber_FR+1)*(i-1)+param2Real, i=1, MUnumber_FR)]                                    
                        else if (trim(param1).eq.(trim(paramTag) // ':' // trim(pool) // '-FF')) then
                            if (MUnumber_FF>0) allocate(paramVec_FF(MUnumber_FF))
                            read(param2,*)param2Real
                            read(param3,*)param3Real
                            paramVec_FF = [((param3Real-param2Real)/(MUnumber_FF+1)*(i-1)+param2Real, i=1, MUnumber_FF)]                                    
                        else if (trim(param1).eq.(trim(paramTag) // ':' // trim(pool) // '-')) then
                            read(param2,*)param2Real
                            read(param3,*)param3Real
                            paramVec = [((param3Real-param2Real)/(Nnumber+1)*(i-1) + param2Real, i=1, Nnumber)]                                     
                            wholePool = .true.
                        end if
                    else if (self%MUParameterDistribution == 'exponential') then                         
                        indexUnits = [(i-1, i = 1, Nnumber)]
                        if (trim(param1).eq.(trim(paramTag) // ':' // trim(pool) // '-S')) then
                            allocate(paramVec_S(2))
                            read(param2,*)param2Real
                            read(param3,*)param3Real
                            paramVec_S = [param2Real, param3Real]
                        else if (trim(param1).eq.(trim(paramTag) // ':' // trim(pool) // '-FR')) then
                            allocate(paramVec_FR(2))
                            read(param2,*)param2Real
                            read(param3,*)param3Real
                            paramVec_FR = [param2Real, param3Real]
                        else if (trim(param1).eq.(trim(paramTag) // ':' // trim(pool) // '-FF')) then
                            allocate(paramVec_FF(2))
                            read(param2,*)param2Real
                            read(param3,*)param3Real
                            paramVec_FF = [param2Real, param3Real]
                        else if (trim(param1).eq.(trim(paramTag) // ':' // trim(pool) // '-')) then
                            read(param2,*)param2Real
                            read(param3,*)param3Real
                            distribute = .false.
                            if (abs(param2Real)>1e-10) then
                                paramVec = param2Real*exp(1.0/Nnumber*log(param3Real/param2Real) * indexUnits)
                            else 
                                paramVec = exp(1.0/Nnumber*log(param3Real + 1.0) * indexUnits) - 1.0
                            end if
                            write(requestedParamater, '(F15.6)')paramVec(index)
                        end if
                    end if
                end do
                
                if (trim(self%MUParameterDistribution).eq.'linear'.and.distribute) then
                    if (MUnumber_FR > 0 .and..not.wholePool) then
                        paramVec(MUnumber_S+1:MUnumber_S+MUnumber_FR) = paramVec_FR
                    end if
                    if (MUnumber_FF > 0 .and..not.wholePool) then 
                        paramVec(MUnumber_S+MUnumber_FR+1:MUnumber_S+MUnumber_FR+MUnumber_FF) = paramVec_FF
                    end if
                    write(requestedParamater, '(F15.6)')paramVec(index)
                else if (self%MUParameterDistribution == 'exponential' .and.distribute) then
                    if (allocated(paramVec_S)) then                        
                        if (paramTag == 'twitchPeak' .or. paramTag == 'bSatSOCDS') then
                            paramVec = paramVec_S(1)*exp(1.0/Nnumber*log(paramVec_FF(2)/paramVec_S(1)) * indexUnits)   
                        else
                            paramVec = ((paramVec_S(1) - (paramVec_S(2)+paramVec_FR(1))/2.0) * exp(-5.0*indexUnits/MUnumber_S)&
                                + ((paramVec_S(2)+paramVec_FR(1))/2.0 - paramVec_FF(2)) &
                                * (1 - exp(1.0/MUnumber_FF*log(((paramVec_FR(2)+paramVec_FF(1))/2.0 - &
                                (paramVec_S(2) + paramVec_FR(1))/2.0)/(paramVec_FF(2)- &
                                (paramVec_S(2)+paramVec_FR(1))/2.0)) * (Nnumber - indexUnits)))&
                                + paramVec_FF(2)) 
                        end if
                        write(requestedParamater, '(F15.6)')paramVec(index)
                    end if
                end if

                if (allocated(paramVec)) deallocate(paramVec)
                if (allocated(paramVec_S)) deallocate(paramVec_S)
                if (allocated(paramVec_FF)) deallocate(paramVec_FF)
                if (allocated(paramVec_FR)) deallocate(paramVec_FR)
            end function parameterSet
            
            type(CharacterMatrix) function determineSynapses(self, neuralSource) result(Synapses)
                ! '''
                ! Function used to determine all the synapses that a given pool makes. It is used in the SynapsesFactory class.
                
                ! - Inputs:
                !     + **neuralSource** - string with the pool name from which is desired to know what synapses it will make.

                ! - Outputs:
                !     + array of strings with all the synapses target that the neuralSource will make.
                ! '''
                class(Configuration), intent(inout) :: self
                character(len=80), intent(in) :: neuralSource
                character(len=80) :: line, param1, param2, param3, paramTag, param
                integer :: ierr, il, j, i, stop1, pos, posUnitKind, posComp, posKind, k
                real(wp) :: paramReal
                type(CharacterArray) :: newSynapse 
                

                Synapses = CharacterMatrix()
                
                paramTag = 'Con:' // trim(neuralSource)
                do k = 1, size(self%confMatrix%item)
                    param1 = self%confMatrix%item(k)%item(1)%string
                    param2 = self%confMatrix%item(k)%item(2)%string
                    param3 = self%confMatrix%item(k)%item(3)%string

                    il = len_trim(param1)
                    pos = 0
                    do i = 1, il
                        if (param1(1:i).eq.paramTag) then
                            pos = i
                            read(param2, *)paramReal
                        end if
                    end do
                                       
                    if ((pos > 0).and.(paramReal > 0.0)) then
                        posUnitKind = 0
                        do i = pos+2, il
                            if (param1(i:i).eq.'-') posUnitKind = i-1
                        end do
                        posComp = 0
                        do i = posUnitKind+1, il
                            if (param1(i:i).eq.'@') posComp = i-1
                        end do
                        posKind = 0
                        do i = posComp, il
                            if (param1(i:i).eq.'|') posKind = i-1
                        end do
                        newSynapse = CharacterArray()
                        param = param1(pos+2:posUnitKind)
                        call newSynapse%AddToList(param)
                        param = param1(posUnitKind+2:posComp)
                        call newSynapse%AddToList(param)
                        param = param1(posComp+2:posKind)
                        call newSynapse%AddToList(param)
                        param = param1(posKind+2:il)
                        call newSynapse%AddToList(param)

                        call Synapses%append(newSynapse)
                    end if               
                end do

                if (.not.allocated(Synapses%item)) allocate(Synapses%item(0))               
                
            end function
            
            
            subroutine changeConfigurationParameter(self, paramTag, value1, value2)
                ! '''
                ! '''
                class(Configuration), intent(inout) :: self
                character(len = 80), intent(in) :: paramTag
                character(len = 80), intent(in) :: value1, value2 
                integer :: i
                logical :: found
                
                found = .false.
                do i = 1, size(self%confMatrix%item)
                    if (self%confMatrix%item(i)%item(1)%string == paramTag) then
                        self%confMatrix%item(i)%item(2)%string = trim(value1)
                        self%confMatrix%item(i)%item(3)%string = trim(value2)
                        found = .true.
                    end if
                end do
                
                ! In case the parameter did not match any tag
                if (.not.found) then
                    print *, "Following parameter tag was not found on configuration file:"
                    print *, paramTag
                    stop 1
                end if
            end subroutine

            subroutine showConfigurationParameter(self, paramTag)
                ! '''
                ! '''
                class(Configuration), intent(inout) :: self
                character(len = 80), intent(in) :: paramTag
                
                integer :: i
                logical :: found
                
                found = .false.
                do i = 1, size(self%confMatrix%item)
                    if (self%confMatrix%item(i)%item(1)%string == paramTag) then
                        print *, paramTag, self%confMatrix%item(i)%item(2)%string, self%confMatrix%item(i)%item(3)%string
                        found = .true.
                    end if
                end do
                
                ! In case the parameter did not match any tag
                if (.not.found) then
                    print *, "Following parameter tag was not found on configuration file:"
                    print *, paramTag
                    stop 1
                end if
            end subroutine

end module

        
        
   
     
    
    
    
        