# Infection Analysis - MPI

## Introduction
Infection Analysis - MPI is a parallel computing simulation that models the spread of a virus in a moving population. Using MPI (Message Passing Interface), the program simulates how individuals interact, infect each other, and recover over time. The simulation operates on a distributed system to efficiently handle large populations and compute infection dynamics.

### Features
- Simulates virus spread based on proximity and exposure time.
- Models individual movement with velocity and direction changes.
- Tracks susceptible, infected, and immune individuals over time.
- Divides the environment into countries for regional infection statistics.
- Uses parallel processing for improved performance on distributed systems.

## Installation
This project requires MPI for execution. The following steps outline how to install and run the simulation.

### Prerequisites
- MPI 
- Tested against OpenMPI (recommended implementation)

### Steps for Setup
1. Install an MPI implementation on your system.
2. Configure environment variables:
   ```sh
   export TMPDIR=/tmp
   export PATH="$HOME/opt/usr/local/bin:$PATH"
   ```
3. Compile the program:
   ```sh
   mpicc -o virus mpi_proj.c utils
   ```
4. Run the simulation:
   ```sh
   mpirun -np <PROCESSES> ./virus <N> <I> <W> <L> <w> <l> <v> <d> <t>
   ```
   Replace the parameters with:
   - `PROCESSES`: Number of parallel processes
   - `N`: Number of individuals
   - `I`: Initial infected individuals
   - `W`: Grid width
   - `L`: Grid length
   - `w`: Country width
   - `l`: Country length
   - `v`: Maximum velocity
   - `d`: Infection distance
   - `t`: Time step for the simulation

## Approach and Assumptions
- The simulation uses MPI for parallel execution across multiple processes.
- Individuals move with random velocities and can change direction upon hitting boundaries.
- The infection model considers exposure time and distance for transmission.
- The environment is modeled as a grid with macro blocks representing countries.

## Algorithm
Each process manages a portion of individuals and follows these steps per time step:
1. Move individuals and update their infection status.
2. Gather infected individuals' positions and broadcast them.
3. Check proximity-based infection rules.
4. Update statuses and output daily statistics.

## Testing and Performance
The program was tested on various hardware configurations with different numbers of processes. Key findings include:
- More processes generally improve performance.
- Performance gains are higher for larger populations.
- Increasing the number of cores speeds up simulations, but gains diminish after a certain threshold.

## Challenges
- Handling large-scale simulations with memory constraints.
- Optimizing performance bottlenecks in data communication between processes.
- Balancing load across processes when individual counts are not evenly divisible.

## Conclusions
The MPI-based infection simulation provides an efficient way to model virus spread in a population. While not perfectly scalable, it shows significant improvements when leveraging distributed computing. Future optimizations could focus on reducing communication overhead and improving load balancing.

## Contributors
- Bruno Morelli
- Cristian Lo Muto
- Vincenzo Martelli

## License
This project is open-source and distributed under the MIT License.

