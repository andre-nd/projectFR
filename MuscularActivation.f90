! '''
!     Neuromuscular simulator in Fortran.
!     Copyright (C) 2018  Renato Naville Watanabe

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

module MuscularActivationClass
    use ConfigurationClass
    use DynamicalArrays
    use MotorUnitClass
    use mkl_spblas
    implicit none
    private
    integer, parameter :: wp = kind( 1.0d0 )
    real(wp), parameter :: PI = 4 * atan(1.0_wp)    
    public :: MuscularActivation

    type MuscularActivation
        type(Configuration), pointer :: conf
        character(len = 6) :: pool
        integer :: MUnumber
        character(len = 80) :: activationModel
        real(wp), dimension(:,:), allocatable :: ActMatrix
        real(wp), dimension(:), allocatable :: an, activation_nonSat, bSat, activation_Sat
        real(wp), dimension(:), allocatable :: diracDeltaValue
        type(MotorUnit), pointer:: unit(:)
        integer :: spIndexing, spRows, spCols, spNumberOfElements, spOperation
        integer, dimension(:), allocatable :: spRowStart, spRowEnd, spColIdx
        real(wp), dimension(:), allocatable :: spValues
        type(sparse_matrix_t) :: ActMatrixSp
        type(matrix_descr) :: spDescr
        real(wp) :: spAlpha, spBeta
        
        contains
            procedure :: atualizeActivationSignal
            procedure :: reset
    end type MuscularActivation
    
    interface MuscularActivation
        module procedure init_MuscularActivation
    end interface MuscularActivation

    contains

        type(MuscularActivation) function init_MuscularActivation(conf, pool, MUnumber, unit)
            class(Configuration), intent(in), target :: conf        
            character(len = 6), intent(in) :: pool
            integer, intent(in) :: MUnumber
            class(MotorUnit), dimension(MUnumber),  intent(in), target:: unit
            character(len = 80) :: paramTag
            integer :: i, j, stat
             

            init_MuscularActivation%conf => conf
            init_MuscularActivation%pool = pool
            init_MuscularActivation%MUnumber = MUnumber
            !allocate(init_MuscularActivation%unit(init_MuscularActivation%MUnumber))
            nullify(init_MuscularActivation%unit)
            init_MuscularActivation%unit => unit(:)
            
            !## Model of the activation signal. For now, it can be *SOCDS* (second order critically damped system).
            paramTag = 'activationModel'
            init_MuscularActivation%activationModel = init_MuscularActivation%conf%parameterSet(paramTag, pool, 0)

            if (trim(init_MuscularActivation%activationModel).eq.'SOCDS') then
                ! ## Matrix that multiplied by the vector formed as the formula below gives the activation
                ! ## signal at instant \f$n\f$:
                ! ## \f{equation}{
                ! ##    \resizebox{0.95\hsize}{!}{$Av(n) = \left[\begin{array}{ccccccccccc}a_1(n-1)&a_1(n-2)&e_1(n-1)&...&a_i(n-i)&a_i(n-2)&e_i(n-1)&...&a__{N_{MU}}(n-1)&a__{N_{MU}}(n-2)&e_{N_{MU}}(n-1)\end{array}\right]^T$}                    
                ! ## \f}
                ! ## where \f$a_i(n)\f$ is the activation signal of the motor unit \f$i\f$, \f$e_i(n)\f$ is
                ! ## 1/T (inverse of simulation time step, Dirac's delta approximation) if the motor unit \f$i\f$,
                ! ## fired at instant \f$n\f$. The vector \f$Av\f$ is updated every step at the function
                ! ## atualizeActivationSignal.
                ! ## The activation matrix itself is formed as:
                ! ## \f{equation}{
                ! ##      \resizebox{0.95\hsize}{!}{$\scriptstyle
                ! ##      A = \left[\begin{array}{ccccccccccc}\scriptscriptstyle  2\exp\left(-\frac{T}{T_{c_1}}\right)&\scriptscriptstyle -\exp\left(-2\frac{T}{T_{c_1}}\right)&\scriptscriptstyle  \frac{T^2}{T_{c_1}}\exp\left(1-\frac{T}{T_{c_1}} \right)&\scriptscriptstyle 0&\scriptscriptstyle ...&\scriptscriptstyle  0&\scriptscriptstyle  0& \scriptscriptstyle 0&\scriptscriptstyle 0&\scriptscriptstyle 0&\scriptscriptstyle 0\\
                ! ##                \scriptscriptstyle 0&\scriptscriptstyle 0&\scriptscriptstyle 0&\scriptscriptstyle \ddots&\scriptscriptstyle ...&&&&&\scriptscriptstyle ...&\scriptscriptstyle 0\\
                ! ##                \scriptscriptstyle 0&\scriptscriptstyle ...&&\scriptscriptstyle 0&\scriptscriptstyle 2\exp\left(-\frac{T}{T_{c_i}}\right)&\scriptscriptstyle -\exp\left(-2\frac{T}{T_{c_i}}\right)&\scriptscriptstyle \frac{T^2}{T_{c_i}}\exp\left(1-\frac{T}{T_{c_i}} \right)&\scriptscriptstyle 0&&&\scriptscriptstyle 0\\
                ! ##                \scriptscriptstyle0&\scriptscriptstyle0&\scriptscriptstyle...&&&\scriptscriptstyle0&\scriptscriptstyle 0&\scriptscriptstyle\ddots&\scriptscriptstyle0&\scriptscriptstyle0\\
                ! ##                \scriptscriptstyle0&\scriptscriptstyle0&\scriptscriptstyle0&\scriptscriptstyle...&&&&\scriptscriptstyle0&\scriptscriptstyle 2\exp\left(-\frac{T}{T_{c_{N_{MU}}}}\right)&\scriptscriptstyle -\exp\left(-2\frac{T}{T_{c_{N_{MU}}}}\right)&\scriptscriptstyle \frac{T^2}{T_{c_{{MU}}}}\exp\left(1-\frac{T}{T_{c_{N_{MU}}}} \right)\end{array}\right]$}
                ! ## \f} 
                ! ## The nonsaturated activation signal \f$a\f$ of all the motor units is obtained with:
                ! ## \f{equation}{
                ! ##   a = A.Av 
                ! ## \f}
                ! ## where each elemement o \f$a\f$ is the activation signal of a motor unit.
                allocate(init_MuscularActivation%ActMatrix(init_MuscularActivation%MUnumber, 3*init_MuscularActivation%MUnumber))
                init_MuscularActivation%ActMatrix(:,:) = 0.0

                do i = 1, init_MuscularActivation%MUnumber
                    init_MuscularActivation%ActMatrix(i,3*(i-1)+1:3*(i-1)+3) = &
                                    [2*exp(-init_MuscularActivation%conf%timeStep_ms/&
                                     init_MuscularActivation%unit(i)%TwitchTc_ms),&
                                     -exp(-2*init_MuscularActivation%conf%timeStep_ms/&
                                     init_MuscularActivation%unit(i)%TwitchTc_ms),& 
                                     (init_MuscularActivation%conf%timeStep_ms**2)/&
                                      init_MuscularActivation%unit(i)%TwitchTc_ms*&
                                      exp(1.0-init_MuscularActivation%conf%timeStep_ms/&
                                      init_MuscularActivation%unit(i)%TwitchTc_ms)]
                end do
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                init_MuscularActivation%spIndexing = 1
                init_MuscularActivation%spRows = size(init_MuscularActivation%ActMatrix, 1)
                init_MuscularActivation%spCols = size(init_MuscularActivation%ActMatrix, 2)
                init_MuscularActivation%spNumberOfElements = 0
                allocate(init_MuscularActivation%spRowStart(init_MuscularActivation%spRows))
                allocate(init_MuscularActivation%spRowEnd(init_MuscularActivation%spRows))
                init_MuscularActivation%spRowStart(:) = 0
                do i = 1, init_MuscularActivation%spRows
                    do j = 1, init_MuscularActivation%spCols
                        if (abs(init_MuscularActivation%ActMatrix(i,j))>1e-10) then
                            init_MuscularActivation%spNumberOfElements = init_MuscularActivation%spNumberOfElements + 1
                            if (init_MuscularActivation%spRowStart(i) == 0) then
                                init_MuscularActivation%spRowStart(i) = init_MuscularActivation%spNumberOfElements
                            end if                    
                            call AddToList(init_MuscularActivation%spValues, init_MuscularActivation%ActMatrix(i,j))
                            call integerAddToList(init_MuscularActivation%spColIdx, j)
                        end if                
                    end do
                    init_MuscularActivation%spRowEnd(i) = init_MuscularActivation%spNumberOfElements + 1
                end do

                ! Create a Sparse Matrix for performance purposes (init_MotorUnitPool%GSp)
                stat = mkl_sparse_d_create_csr(init_MuscularActivation%ActMatrixSp, &
                                               init_MuscularActivation%spIndexing, &
                                               init_MuscularActivation%spRows, &
                                               init_MuscularActivation%spCols, &
                                               init_MuscularActivation%spRowStart, &
                                               init_MuscularActivation%spRowEnd, &
                                               init_MuscularActivation%spColIdx, &
                                               init_MuscularActivation%spValues)

                ! Values for the matrix-vector operation matInt = GV
                init_MuscularActivation%spDescr%type = SPARSE_MATRIX_TYPE_GENERAL
                init_MuscularActivation%spAlpha = 1.0 
                init_MuscularActivation%spOperation = SPARSE_OPERATION_NON_TRANSPOSE 
                init_MuscularActivation%spBeta = 0.0    
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    
                ! ## Is a vector formed as:
                ! ## \f{equation}{
                ! ##    \resizebox{0.95\hsize}{!}{$Av(n) = \left[\begin{array}{ccccccccccc}a_1(n-1)&a_1(n-2)&e_1(n-1)&...&a_i(n-i)&a_i(n-2)&e_i(n-1)&...&a__{N_{MU}}(n-1)&a__{N_{MU}}(n-2)&e_{N_{MU}}(n-1)\end{array}\right]^T$}                    
                ! ## \f}
                ! ## It is multiplied by the matriz actMatrix to obtain the activation signal 
                ! ## (see actMatrix explanation)
                allocate(init_MuscularActivation%an(3*init_MuscularActivation%MUnumber))
                init_MuscularActivation%an(:) = 0.0
            end if

            ! ## The non-saturated activation signal of all motor units (see actMatrix explanation).
            allocate(init_MuscularActivation%activation_nonSat(init_MuscularActivation%MUnumber))    
            init_MuscularActivation%activation_nonSat(:) = 0.0
            ! ## The parameter \f$b\f$ (see twitchSaturation function explanation) of 
            ! ## each motor unit.
            allocate(init_MuscularActivation%bSat(init_MuscularActivation%MUnumber))
            do i = 1, init_MuscularActivation%MUnumber     
                init_MuscularActivation%bSat(i) = init_MuscularActivation%unit(i)%bSat
            end do
                       
            ! ## The non-saturated activation signal of all motor units (see actMatrix explanation).
            allocate(init_MuscularActivation%activation_Sat(init_MuscularActivation%MUnumber))
            init_MuscularActivation%activation_Sat(:) = 0.0
            ! ## Dirac's delta approximation amplitude value. Is the inverse
            ! ## of the simulation time step (\f$1/T\f$). 
            init_MuscularActivation%diracDeltaValue = - init_MuscularActivation%bSat / init_MuscularActivation%conf%timeStep_ms

            !self.MUindices = np.arange(0, self.MUnumber)
        end function

        subroutine atualizeActivationSignal(self, t)
            ! '''
            ! Update the activation signal of the motor units.
            !
            ! - Inputs:
            !     + **t**: current instant, in ms.        
            ! '''
            class(MuscularActivation), intent(inout) :: self
            real(wp), intent(in) :: t
            !class(MotorUnit), intent(in), dimension(self%MUnumber)  :: unit
            integer :: i, sizeTrain, stat            

            do i = 1, self%MUnumber
                self%an(3*(i-1)+2) = self%an(3*(i-1)+1)
                self%an(3*(i-1)+1) = self%activation_nonSat(i)
                if (allocated(self%unit(i)%terminalSpikeTrain)) then                    
                    sizeTrain = size(self%unit(i)%terminalSpikeTrain)                    
                    if (abs(t - self%unit(i)%terminalSpikeTrain(sizeTrain)) < 1e-6) then
                        self%an(3*(i-1)+3) = self%diracDeltaValue(i)
                    else
                        self%an(3*(i-1)+3) = 0.0
                    end if   
                else
                    self%an(3*(i-1)+3) = 0.0
                end if
            end do  

             ! Create a Sparse Matrix for performance purposes (init_MotorUnitPool%GSp)
            stat = mkl_sparse_d_create_csr(self%ActMatrixSp, &
                                           self%spIndexing, &
                                           self%spRows, &
                                           self%spCols, &
                                           self%spRowStart, &
                                           self%spRowEnd, &
                                           self%spColIdx, &
                                           self%spValues)
          

            stat = mkl_sparse_d_mv(self%spOperation, &
                                   self%spAlpha, &
                                   self%ActMatrixSp, &
                                   self%spDescr, &
                                   self%an, &
                                   self%spBeta, &
                                   self%activation_nonSat)

            ! self%activation_nonSat = matmul(self%ActMatrix, self%an)
            
            ! \f{equation}{
            ! a_{sat} = \frac{1-e^{-b.a_{nSat}}}{1+e^{-b.a_{nSat}}}
            ! \f}
            self%activation_Sat(:) = 2.0 / (1.0 + exp(self%activation_nonSat)) - 1.0
            stat = mkl_sparse_destroy(self%ActMatrixSp)
        end subroutine

        subroutine reset(self)
            ! '''

            ! '''
            class(MuscularActivation), intent(inout) :: self
            
            self%an(:) = 0.0
            self%activation_nonSat(:) = 0.0
            self%activation_Sat(:) = 0.0
        end subroutine

end module MuscularActivationClass

    
    
