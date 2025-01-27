ifort -c ../../queue.f90
ifort -c ../../ogpf.f90
ifort -c ../../DynamicalArrays.f90
ifort -c ../../CharacterArray.f90
ifort -c ../../CharacterMatrix.f90
ifort -c ../../randomSeedInitialize.f90
ifort -c ../../Configuration.f90
ifort -c ../../PulseConductanceState.f90
ifort -c ../../ChannelConductance.f90
ifort -c ../../Synapse.f90
ifort -c ../../SynapsePointer.f90
ifort -c ../../Compartment.f90
ifort -c ../../AxonDelay.f90
ifort -c ../../MotorUnit.f90
ifort -c ../../MuscularActivation.f90
ifort -c ../../MuscleNoHill.f90
ifort -c ../../MuscleHill.f90
ifort -c ../../MuscleSpindle.f90
ifort -c ../../MotorUnitPool.f90
ifort -c ../../PointProcessGenerator.f90
ifort -c ../../NeuralTractUnit.f90
ifort -c ../../NeuralTract.f90
ifort -c ../../Interneuron.f90
ifort -c ../../InterneuronPool.f90
ifort -c ../../SynapticNoise.f90
ifort -c ../../SynapsesFactory.f90

ifort NeuralTractExample.f90 -o NeuralTractExample -O3 ../../SynapticNoise.f90 ../../MuscleSpindle.f90 ../../MuscleHill.f90 ../../Interneuron.f90 ../../InterneuronPool.f90 ../../SynapsePointer.f90 ../../SynapsesFactory.f90 ../../Synapse.f90 ../../CharacterMatrix.f90 ../../PointProcessGenerator.f90 ../../NeuralTractUnit.f90 ../../NeuralTract.f90  ../../queue.f90 ../../CharacterArray.f90 ../../MuscleNoHill.f90 ../../MuscularActivation.f90 ../../MotorUnitPool.f90 ../../AxonDelay.f90 ../../Compartment.f90 ../../MotorUnit.f90   ../../DynamicalArrays.f90 ../../randomSeedInitialize.f90 ../../ogpf.f90 ../../Configuration.f90 ../../PulseConductanceState.f90 ../../ChannelConductance.f90 


