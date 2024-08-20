# Use your own image from Docker Hub as the base image
FROM grayfoxdevhaus/libfec:latest AS build
#FROM grayfoxdevhaus/libfec:local AS build

# Install necessary packages
RUN apk update && apk add --no-cache \
    build-base \
    gcc \
    g++ \
    libstdc++ \
    cmake  \
    bash

# Set the working directory
WORKDIR /app

# Copy the source code into the container
COPY src/ref_test.cpp .
COPY src/codec.cpp .
COPY include/rs /usr/local/include

# Compile the C++ code
# Non-Debug
RUN g++ -v -O3 -std=c++23 -I/usr/local/include -L/usr/local/lib -o ref_test ref_test.cpp -lfec
# Debug
#RUN g++ -v -g -O0 -march=native -std=c++23 -I/usr/local/include -L/usr/local/lib -o ref_test ref_test.cpp -lfec

# Override any existing entrypoint and set command to run the executable
ENTRYPOINT []
CMD ["./ref_test"]
#CMD ["/bin/bash"]
