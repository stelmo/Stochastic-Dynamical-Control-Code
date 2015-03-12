## Instructions for running the code in this folder

Point Julia's working directory to this folder using `cd("~\CSTR_Model")`. Then:

1. Execute `include("Reactor_test.jl")` to run the integration (Runge-Kutta 4th order) test.
2. Execute `include("Reactor_Qual.jl")` to view qualitative analysis plots.
3. Execute `include("Reactor.jl")` to run some simulations on the reactor. Edit this script to change the reactor parameters.
4. Execute `include("Reactor_Lin.jl")` to run the linearisation procedure around an arbitrary number of linearisation points.
