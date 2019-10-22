gfortran-8 -c ../../queue.f90
gfortran-8 -c ../../ogpf.f90
gfortran-8 -c ../../DynamicalArrays.f90
gfortran-8 -c ../../String.f90
gfortran-8 -c ../../CharacterArray.f90
gfortran-8 -c ../../CharacterMatrix.f90
gfortran-8 -c ../../randomSeedInitialize.f90
gfortran-8 -c ../../randomGen.f90
gfortran-8 -c ../../Configuration.f90
gfortran-8 -c ../../PulseConductanceState.f90
gfortran-8 -c ../../ChannelConductance.f90
gfortran-8 -c ../../Synapse.f90
gfortran-8 -c ../../SynapsePointer.f90
gfortran-8 -c ../../Compartment.f90
gfortran-8 -c ../../AxonDelay.f90
gfortran-8 -c ../../MotorUnit.f90
gfortran-8 -c $MKLROOT/include/mkl_spblas.f90
gfortran-8 -c ../../MuscularActivation.f90
gfortran-8 -c ../../MuscleNoHill.f90
gfortran-8 -c ../../MuscleHill.f90
gfortran-8 -c ../../MuscleSpindle.f90
gfortran-8 -c ../../GolgiTendonOrgan.f90
gfortran-8 -c ../../MotorUnitPool.f90
gfortran-8 -c ../../PointProcessGenerator.f90
gfortran-8 -c ../../NeuralTractUnit.f90
gfortran-8 -c ../../NeuralTract.f90
gfortran-8 -c ../../Interneuron.f90
gfortran-8 -c ../../InterneuronPool.f90
gfortran-8 -c ../../SynapticNoise.f90
gfortran-8 -c ../../AfferentUnit.f90
gfortran-8 -c ../../AfferentPool.f90
gfortran-8 -c ../../SynapsesFactory.f90
gfortran-8 -c ../../MusclePointer.f90
gfortran-8 -c ../../jointAnkleForceTask.f90
gfortran-8 -c ../../postureTask.f90

gfortran-8 MotorUnitPoolWithDescendingCommand.f90 -o MotorUnitPoolWithDescendingCommand -O3 -g -fbacktrace ../../GolgiTendonOrgan.f90 ../../MusclePointer.f90 ../../jointAnkleForceTask.f90 ../../AfferentUnit.f90 ../../AfferentPool.f90 ../../String.f90 ../../SynapticNoise.f90 ../../MuscleSpindle.f90 ../../MuscleHill.f90 ../../Interneuron.f90 ../../InterneuronPool.f90 ../../SynapsePointer.f90 ../../SynapsesFactory.f90 ../../Synapse.f90 ../../CharacterMatrix.f90 ../../PointProcessGenerator.f90 ../../NeuralTractUnit.f90 ../../NeuralTract.f90  ../../queue.f90 ../../CharacterArray.f90 ../../MuscleNoHill.f90 ../../MuscularActivation.f90 ../../MotorUnitPool.f90 ../../AxonDelay.f90 ../../Compartment.f90 ../../MotorUnit.f90   ../../DynamicalArrays.f90 ../../randomSeedInitialize.f90  ../../randomGen.f90 ../../ogpf.f90 ../../Configuration.f90 ../../PulseConductanceState.f90 ../../ChannelConductance.f90 ../../postureTask.f90  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl -m64 -I$MKLROOT/include


