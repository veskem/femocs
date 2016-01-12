################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/AtomReader.cpp \
../src/Femocs.cpp \
../src/Medium.cpp \
../src/Mesher.cpp \
../src/SurfaceExtractor.cpp 

OBJS += \
./src/AtomReader.o \
./src/Femocs.o \
./src/Medium.o \
./src/Mesher.o \
./src/SurfaceExtractor.o 

CPP_DEPS += \
./src/AtomReader.d \
./src/Femocs.d \
./src/Medium.d \
./src/Mesher.d \
./src/SurfaceExtractor.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__cplusplus=201103L -O3 -g3 -Wall -c -fmessage-length=0 -std=c++11 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


