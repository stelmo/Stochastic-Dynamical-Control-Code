## Instructions for running the code in this folder

Point Julia's working directory to this folder using `cd("~\Linear_Latent_Dynamical_Models")`. Then:

1. Execute `include("LLDS_tests.jl")` to run the tests.
2. Execute `include("Linear_Reactor.jl")` to run filtering and prediction on the CSTR reactor model. The nonlinear model is used to generate observations while the linear model is used to perform inference. Edit this script if you would like to change parameters.

### External Contributions

The code found in the module `Confidence.jl` is based entirely on the Matlab error ellipse code found [here](http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/).
