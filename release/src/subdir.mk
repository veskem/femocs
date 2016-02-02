################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/AtomReader.cpp \
../src/Femocs.cpp \
../src/Medium.cpp \
../src/Mesh.cpp \
../src/Mesher.cpp \
../src/Surface.cpp \
../src/SurfaceExtractor.cpp \
../src/Vacuum.cpp 

OBJS += \
./src/AtomReader.o \
./src/Femocs.o \
./src/Medium.o \
./src/Mesh.o \
./src/Mesher.o \
./src/Surface.o \
./src/SurfaceExtractor.o \
./src/Vacuum.o 

CPP_DEPS += \
./src/AtomReader.d \
./src/Femocs.d \
./src/Medium.d \
./src/Mesh.d \
./src/Mesher.d \
./src/Surface.d \
./src/SurfaceExtractor.d \
./src/Vacuum.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__cplusplus=201103L -O3 -g3 -Wall -c -fmessage-length=0 -std=c++11 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


